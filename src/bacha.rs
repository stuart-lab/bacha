use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use flate2::read::MultiGzDecoder;
use flate2::Compression;
use gzp::{deflate::Gzip, par::compress::ParCompressBuilder};
use lexical_core::parse;
use log::error;
use log::info;
use ndarray::Array2;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rust_lapper::{Interval, Lapper};
use rustc_hash::FxHashMap;
use smallvec::SmallVec;
use std::{
    collections::HashSet,
    error::Error,
    fs::File,
    io,
    io::BufRead,
    io::BufReader,
    io::Read,
    io::Write,
    path::Path,
    path::PathBuf,
    sync::atomic::{AtomicBool, AtomicU64, Ordering},
    sync::Mutex,
    thread,
    time::Instant,
};
use xz2::read::XzDecoder;

// A parsed fragment routed from a parser worker to the shard owner.
// One Vec<ForwardRecord> per chunk-flush keeps channel traffic to ~O(N)
// sends per chunk rather than O(N_lines).
#[derive(Clone, Copy)]
struct ForwardRecord {
    local_idx: u32,
    chrom_id: u8,
    start: u32,
    end: u32,
}

// Open the fragment file, auto-detecting BGZF. BGZF (Blocked GZIP, used by
// CellRanger and the htslib ecosystem) is a sequence of independent gzip
// blocks and can be inflated in parallel; we delegate to noodles-bgzf's
// MultithreadedReader when `bgzf_workers >= 1`. Plain gzip (single deflate
// stream) cannot be decoded in parallel and falls back to MultiGzDecoder.
//
// Detection follows the BGZF spec: gzip magic 1F 8B, FEXTRA flag bit set in
// the header flags byte, and the well-known "BC" extra subfield ID at offset
// 12. Anything else is treated as plain gzip.
fn is_bgzf_file(path: &Path) -> io::Result<bool> {
    let mut probe = File::open(path)?;
    let mut header = [0u8; 16];
    let n = probe.read(&mut header)?;
    Ok(n >= 16
        && header[0] == 0x1f
        && header[1] == 0x8b
        && (header[3] & 0x04) != 0
        && header[12] == b'B'
        && header[13] == b'C')
}

const DEFAULT_BGZF_WORKERS: usize = 2;

struct ThreadBudget {
    bgzf: usize,
    counting: usize,
}

fn thread_budget(num_threads: usize, is_bgzf: bool) -> ThreadBudget {
    let total = num_threads.max(1);
    // At -t 1 the single-threaded path runs, which has no BGZF pool at all.
    if !is_bgzf || total == 1 {
        return ThreadBudget {
            bgzf: 0,
            counting: total,
        };
    }
    let requested = std::env::var("BACHA_BGZF_WORKERS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(DEFAULT_BGZF_WORKERS)
        .max(1);
    let bgzf = requested.min(total - 1);
    ThreadBudget {
        bgzf,
        counting: total - bgzf,
    }
}

// Decode threads for the index file. Index loading finishes before counting
// starts, so it has the budget to itself and takes the same share, bounded by
// what was requested.
fn bgzf_worker_count(num_threads: usize) -> usize {
    match std::env::var("BACHA_BGZF_WORKERS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
    {
        Some(n) => n.max(1),
        None => DEFAULT_BGZF_WORKERS.min(num_threads.max(1)),
    }
}

fn open_fragment_reader(path: &Path, bgzf_workers: usize) -> io::Result<Box<dyn BufRead + Send>> {
    let file = File::open(path)?;
    if is_bgzf_file(path)? {
        let workers =
            std::num::NonZeroUsize::new(bgzf_workers.max(1)).unwrap_or(std::num::NonZeroUsize::MIN);
        eprintln!(
            "Detected BGZF fragment file; using parallel BGZF decoder ({} workers)",
            workers
        );
        let reader = noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, file);
        Ok(Box::new(BufReader::with_capacity(8 * 1024 * 1024, reader)))
    } else {
        eprintln!("Fragment file is not BGZF; using single-threaded gzip decoder");
        Ok(Box::new(BufReader::with_capacity(
            8 * 1024 * 1024,
            MultiGzDecoder::new(file),
        )))
    }
}

// The 28-byte empty BGZF block that terminates a well-formed BGZF file (see
// the BGZF section of the SAM specification). `bgzip` and every htslib-based
// writer append it when the file is closed, so its absence means the file is
// truncated or still being written.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

// Refuse a fragment file that is visibly incomplete. Called once up front so
// the check applies to every counting path (single-threaded and sharded alike)
// and costs one small read rather than being entangled with decoder setup.
fn check_fragment_file_complete(path: &Path) -> io::Result<()> {
    if is_bgzf_file(path)? {
        check_bgzf_eof(path)?;
    }
    Ok(())
}

fn check_bgzf_eof(path: &Path) -> io::Result<()> {
    use std::io::{Seek, SeekFrom};

    let mut file = File::open(path)?;
    let len = file.seek(SeekFrom::End(0))?;
    if len >= BGZF_EOF.len() as u64 {
        file.seek(SeekFrom::End(-(BGZF_EOF.len() as i64)))?;
        let mut tail = [0u8; BGZF_EOF.len()];
        file.read_exact(&mut tail)?;
        if tail == BGZF_EOF {
            return Ok(());
        }
    }
    Err(io::Error::new(
        io::ErrorKind::UnexpectedEof,
        format!(
            "{} is BGZF but does not end with the BGZF end-of-file marker, so it is \
             truncated or still being written.",
            path.display()
        ),
    ))
}

// Report lines that were read but could not be used. Both are warnings rather
// than errors: a handful of malformed lines is a data-quality issue for the
// user to judge
fn report_skipped_lines(malformed: u64, total_lines: u64, unterminated_final_line: bool) {
    if malformed > 0 {
        eprintln!(
            "WARNING: skipped {} unparseable line(s) out of {} ({:.2}%): fewer than 4 \
             tab-separated fields, a non-UTF-8 chromosome or barcode, or a non-numeric \
             start/end coordinate.",
            malformed,
            total_lines,
            100.0 * malformed as f64 / total_lines.max(1) as f64,
        );
    }
    if unterminated_final_line {
        eprintln!(
            "WARNING: the fragment file does not end with a newline. Its final line was \
             counted, but this can also mean the file is truncated; verify the fragment \
             count against the source file."
        );
    }
}

// flat structure to store indices efficiently
pub struct FlatIndVec {
    data: Vec<u32>,    // All indices concatenated
    offsets: Vec<u32>, // Start position of each interval's indices
}

impl FlatIndVec {
    #[inline(always)]
    fn get(&self, index: usize) -> Option<&[u32]> {
        if index >= self.offsets.len() - 1 {
            return None;
        }
        let start = self.offsets[index] as usize;
        let end = self.offsets[index + 1] as usize;
        Some(&self.data[start..end])
    }
}

// Add one count to every feature column an index interval feeds. A site
// interval feeds exactly one column; a control interval feeds one column per
// feature that selected it.
//
// The feature index carried by the index file is 1-based. `build_index_result`
// validates at load time that every site index lies in 1..=n_tfs and every
// control index in n_tfs+1..=2*n_tfs, so `- 1` always yields a valid column of a
// row of length `n_features` (= 2 * n_tfs) and no per-hit bounds test is needed.
#[inline(always)]
fn bump_interval(
    interval_idx: usize,
    row: &mut [u32],
    site_indices: &[u32],
    control_indices: &FlatIndVec,
    num_site_intervals: usize,
) {
    if interval_idx < num_site_intervals {
        let col = site_indices[interval_idx] as usize - 1;
        row[col] += 1;
    } else {
        let control_idx = interval_idx - num_site_intervals;
        if let Some(indices) = control_indices.get(control_idx) {
            for &raw in indices {
                let col = raw as usize - 1;
                row[col] += 1;
            }
        }
    }
}

// Count one fragment into `row`, the contiguous count slice for its cell. The
// caller obtains that slice (e.g. via `Array2::row_mut(cell_idx).into_slice()`
// or `ArrayViewMut2::row_mut(...)`), which lets this routine serve the
// single-threaded path, the sharded workers, and the forwarded-record drainer
// without duplicating the inner loop.
//
// `ins1` and `ins2` are the fragment's two 0-based Tn5 insertion coordinates.
//
// Paired insertion counting: an interval containing either insertion is
// counted once
#[inline(always)]
fn process_fragment(
    ins1: u32,
    ins2: u32,
    row: &mut [u32],
    lapper: &Lapper<u32, u32>,
    site_indices: &[u32],
    control_indices: &FlatIndVec,
    num_site_intervals: usize,
    n_features: usize,
) {
    debug_assert!(row.len() >= n_features);

    for interval in lapper.find(ins1, ins1 + 1) {
        bump_interval(
            interval.val as usize,
            row,
            site_indices,
            control_indices,
            num_site_intervals,
        );
    }
    for interval in lapper.find(ins2, ins2 + 1) {
        // Already counted via ins1.
        if interval.start <= ins1 && ins1 < interval.stop {
            continue;
        }
        bump_interval(
            interval.val as usize,
            row,
            site_indices,
            control_indices,
            num_site_intervals,
        );
    }
}

pub fn bacha(
    frag_file: &Path,
    index_file: &Path,
    cell_file: &Path,
    output_path: &Path,
    num_threads: usize,
    lambda: f32,
    pseudocount: f32,
    check_bounds: bool,
) -> Result<(), Box<dyn Error>> {
    info!("Received fragment file: {:?}", frag_file);
    info!("Received cell file: {:?}", cell_file);
    info!("Received index file: {:?}", index_file);
    info!("Received output file: {:?}", output_path);

    precheck(output_path, num_threads, lambda, pseudocount)?;

    let (cells, cell_names) = read_cells_file(cell_file).map_err(|e| {
        error!("Failed to read cell file: {}", e);
        e
    })?;
    let sources = vec![FragmentSource {
        sample: frag_file
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| frag_file.display().to_string()),
        path: frag_file.to_path_buf(),
        cells,
        cell_names,
    }];

    bacha_sources(
        &sources,
        index_file,
        output_path,
        num_threads,
        lambda,
        pseudocount,
        check_bounds,
    )
}

/// Score several fragment files listed in a manifest into one output matrix.
///
/// The manifest is a tab-separated file of `sample`, `fragments`, `cells` rows
/// (see `read_manifest`). Files are counted one at a time, each using all
/// `num_threads` workers, into a single matrix covering every cell in every
/// file; the index is loaded once for the whole run.
pub fn bacha_manifest(
    manifest_file: &Path,
    index_file: &Path,
    output_path: &Path,
    num_threads: usize,
    lambda: f32,
    pseudocount: f32,
    check_bounds: bool,
) -> Result<(), Box<dyn Error>> {
    info!("Received manifest file: {:?}", manifest_file);
    info!("Received index file: {:?}", index_file);
    info!("Received output file: {:?}", output_path);

    precheck(output_path, num_threads, lambda, pseudocount)?;

    let rows = read_manifest(manifest_file)?;
    eprintln!("Manifest lists {} fragment file(s)", rows.len());
    let sources = build_sources(rows)?;

    bacha_sources(
        &sources,
        index_file,
        output_path,
        num_threads,
        lambda,
        pseudocount,
        check_bounds,
    )
}

// Cheap argument validation, run before any input file is opened so that a bad
// thread count or an output path that already exists fails immediately rather
// than after minutes of counting.
fn precheck(
    output_path: &Path,
    num_threads: usize,
    lambda: f32,
    pseudocount: f32,
) -> io::Result<()> {
    validate_thread_count(num_threads)?;
    validate_scoring_parameters(lambda, pseudocount)?;

    // Refuse to clobber an existing output file. Reported as an error rather
    // than exiting the process directly, so callers embedding this function
    // can handle it and `main` formats it like every other failure.
    if output_path.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!("Output file already exists: {}", output_path.display()),
        ));
    }
    Ok(())
}

fn bacha_sources(
    sources: &[FragmentSource],
    index_file: &Path,
    output_path: &Path,
    num_threads: usize,
    lambda: f32,
    pseudocount: f32,
    check_bounds: bool,
) -> Result<(), Box<dyn Error>> {
    // Run everything inside a rayon pool sized to `num_threads`
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;

    // Returns io::Result, not Box<dyn Error>: rayon requires the closure's
    // return type to be Send, and Box<dyn Error> is not.
    let t0 = Instant::now();
    pool.install(|| {
        run_count(
            sources,
            index_file,
            output_path,
            num_threads,
            lambda,
            pseudocount,
            check_bounds,
        )
    })?;
    if std::env::var("BACHA_PROFILE").is_ok() {
        eprintln!("[PERF] total={:.2}s", t0.elapsed().as_secs_f64());
    }
    Ok(())
}

// Deterministic genome/assembly mismatch error. A Tn5 insertion cannot fall
// beyond the physical end of a chromosome, so if any coordinate observed on a
// chromosome exceeds that chromosome's length in the index, the fragment file
// and index must be built on different genome assemblies (e.g. an mm10 index
// scored against hg38 fragments).
//
// This is checked inline as fragments are read (see `fcount_core`) so the run
// aborts on the first offending coordinate rather than processing the entire
// fragment file first. Every index carries a length for each chromosome in its
// chromosome table; `fcount_core` substitutes `u32::MAX` when the caller passes
// `--no-bounds-check`, which neutralises the comparison without branching in
// the hot loop.
fn assembly_mismatch_error(chrom: &str, observed: u32, index_len: u32) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        format!(
            "Genome mismatch: fragment coordinate {observed} on {chrom} lies beyond the end of \
             that chromosome in the index (length {index_len}), which is impossible unless the \
             fragment file and index were built on different genome assemblies (e.g. an mm10 \
             index scored against hg38 fragments).",
        ),
    )
}

// Verify that the fragment file and the index/cell file actually describe the
// same data.
//
// `matched_cells` counts fragments whose barcode was found in the cell file.
// Those are the only fragments that go on to the chromosome lookup, so it — not
// the raw parsed-line count — is the denominator for the unmatched-chromosome
// rate, and `unmatched_chroms <= matched_cells` always holds.
//
// A completely empty fragment file is left to the caller: it is already obvious
// from the "Finished reading 0 total fragments" line and is not a mismatch.
fn check_input_overlap(
    total_fragments: u64,
    matched_cells: u64,
    unmatched_chroms: u64,
) -> io::Result<()> {
    if total_fragments > 0 && matched_cells == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No fragments matched any cell barcode in the cell file.",
        ));
    }
    if unmatched_chroms == 0 {
        return Ok(());
    }
    if unmatched_chroms == matched_cells {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No fragments matched any chromosome in the index. \
             The fragment file and index likely use different genome assemblies \
             or chromosome naming conventions (e.g. 'chr1' vs '1').",
        ));
    }
    let pct = 100.0 * unmatched_chroms as f64 / matched_cells as f64;
    if pct > 20.0 {
        eprintln!(
            "WARNING: {} of {} fragments with a matching cell barcode ({:.1}%) were on \
             chromosomes not present in the index. This may indicate a genome mismatch \
             between the fragment file and the index.",
            unmatched_chroms, matched_cells, pct
        );
    }
    Ok(())
}

// One fragment file and the cells to score from it.
//
// `cells` maps a barcode as it appears in this fragment file to a row index
// *local to this source* (0..cell_names.len()). The driver turns those into
// global rows by adding the source's offset into the combined count matrix, so
// a source never needs to know where it sits in the output.
struct FragmentSource {
    // Used only to attribute messages and errors to a file; not part of any
    // output cell name unless the caller put it there.
    sample: String,
    path: PathBuf,
    cells: FxHashMap<String, u32>,
    cell_names: Vec<String>,
}

// Counts accumulated while reading one fragment file. Kept separate from the
// matrix so the driver can run the genome/barcode sanity checks per file and
// name the offending sample, rather than reporting one figure for the batch.
#[derive(Default)]
struct FileStats {
    total_fragments: u64,
    matched_cells: u64,
    unmatched_chroms: u64,
    malformed: u64,
    total_lines: u64,
    unterminated_final_line: bool,
}

// Attribute an error to a sample. Single-file runs pass `None` so their
// messages stay byte-identical to previous releases (and to what the tests
// assert); only a manifest run adds the prefix.
fn with_sample_context(sample: Option<&str>, e: io::Error) -> io::Error {
    match sample {
        Some(s) => io::Error::new(e.kind(), format!("sample '{s}': {e}")),
        None => e,
    }
}

// Read a cells file.
//
// Each line is either a bare barcode, or a barcode and the name to give that
// cell in the output, tab-separated:
//
//     AAACGAAAGATTGGCT-1
//     AAACGAAAGATTGGCT-1<TAB>sample1_AAACGAAAGATTGGCT-1
//
// Returns (barcode -> local row index, output names in file order).
fn read_cells_file(cell_file: &Path) -> io::Result<(FxHashMap<String, u32>, Vec<String>)> {
    let cellreader = File::open(cell_file).map(BufReader::new)?;

    // barcode -> output name, to catch a barcode claiming two different output
    // cells; and output name -> barcode, to catch two barcodes claiming one
    // output cell. Either would silently drop a column's counts otherwise.
    let mut by_barcode: FxHashMap<String, String> = FxHashMap::default();
    let mut by_name: FxHashMap<String, String> = FxHashMap::default();

    let mut cells: FxHashMap<String, u32> = FxHashMap::default();
    let mut cell_names: Vec<String> = Vec::new();

    for (lineno, line) in cellreader.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        // split() on a non-empty string always yields at least one field.
        let barcode = fields.next().unwrap();
        let name = fields.next().filter(|s| !s.is_empty()).unwrap_or(barcode);
        if fields.next().is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}: line {} has more than 2 tab-separated fields. A cells file holds \
                     either one barcode per line, or a barcode and an output cell name.",
                    cell_file.display(),
                    lineno + 1
                ),
            ));
        }
        if barcode.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}: line {} has an empty barcode",
                    cell_file.display(),
                    lineno + 1
                ),
            ));
        }

        // A line repeated verbatim is deduplicated silently, as it always has
        // been. A barcode or output name reused with different partners is an
        // error, because one of the two would otherwise be discarded.
        if let Some(prev) = by_barcode.get(barcode) {
            if prev == name {
                continue;
            }
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}: barcode '{}' is mapped to two different output cell names \
                     ('{}' and '{}')",
                    cell_file.display(),
                    barcode,
                    prev,
                    name
                ),
            ));
        }
        if let Some(prev) = by_name.get(name) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}: output cell name '{}' is claimed by two different barcodes \
                     ('{}' and '{}')",
                    cell_file.display(),
                    name,
                    prev,
                    barcode
                ),
            ));
        }

        by_barcode.insert(barcode.to_string(), name.to_string());
        by_name.insert(name.to_string(), barcode.to_string());
        cells.insert(barcode.to_string(), cell_names.len() as u32);
        cell_names.push(name.to_string());
    }

    if cell_names.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Cell file is empty",
        ));
    }

    info!("number of cells: {}", cell_names.len());
    Ok((cells, cell_names))
}

// A parsed manifest row, with paths already resolved.
#[derive(Debug)]
struct ManifestRow {
    sample: String,
    fragments: PathBuf,
    cells: PathBuf,
}

// Resolve a manifest path: absolute paths are taken as given, relative ones
// are interpreted against the manifest's own directory. Resolving against the
// manifest rather than the process working directory is what makes a manifest
// portable - it can sit next to the files it names and be run from anywhere.
fn resolve_manifest_path(base: &Path, raw: &str) -> PathBuf {
    let p = Path::new(raw);
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        base.join(p)
    }
}

// Parse a tab-separated manifest of fragment files:
//
//     sample<TAB>fragments<TAB>cells
//     FN1<TAB>lung/FN1_atac_fragments.tsv.gz<TAB>FN1_cells.txt
//
// A leading header row naming exactly those three columns is optional and
// skipped if present. Blank lines and lines starting with '#' are ignored.
fn read_manifest(manifest_file: &Path) -> io::Result<Vec<ManifestRow>> {
    let base = manifest_file
        .parent()
        .unwrap_or(Path::new("."))
        .to_path_buf();
    let reader = File::open(manifest_file).map(BufReader::new)?;

    let mut rows: Vec<ManifestRow> = Vec::new();
    let mut seen_samples: FxHashMap<String, usize> = FxHashMap::default();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').map(|f| f.trim()).collect();
        if fields.len() != 3 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}: line {} has {} tab-separated field(s); a manifest row must have 3 \
                     (sample, fragments, cells)",
                    manifest_file.display(),
                    lineno + 1,
                    fields.len()
                ),
            ));
        }

        // Optional header, recognised only as the first data row so that a
        // sample legitimately named "sample" further down is not skipped.
        if rows.is_empty()
            && fields[0].eq_ignore_ascii_case("sample")
            && fields[1].eq_ignore_ascii_case("fragments")
            && fields[2].eq_ignore_ascii_case("cells")
        {
            continue;
        }

        if fields.iter().any(|f| f.is_empty()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}: line {} has an empty field; sample, fragments and cells are all \
                     required",
                    manifest_file.display(),
                    lineno + 1
                ),
            ));
        }

        if let Some(prev) = seen_samples.insert(fields[0].to_string(), lineno + 1) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{}: sample '{}' appears on both line {} and line {}; sample names must \
                     be unique",
                    manifest_file.display(),
                    fields[0],
                    prev,
                    lineno + 1
                ),
            ));
        }

        rows.push(ManifestRow {
            sample: fields[0].to_string(),
            fragments: resolve_manifest_path(&base, fields[1]),
            cells: resolve_manifest_path(&base, fields[2]),
        });
    }

    if rows.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{}: manifest contains no rows", manifest_file.display()),
        ));
    }
    Ok(rows)
}

// Load the cells for every manifest row, turning rows into sources.
fn build_sources(rows: Vec<ManifestRow>) -> io::Result<Vec<FragmentSource>> {
    let mut sources = Vec::with_capacity(rows.len());
    for row in rows {
        let (cells, cell_names) = read_cells_file(&row.cells).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!(
                    "sample '{}': cannot read cells file '{}': {e}",
                    row.sample,
                    row.cells.display()
                ),
            )
        })?;
        sources.push(FragmentSource {
            sample: row.sample,
            path: row.fragments,
            cells,
            cell_names,
        });
    }
    Ok(sources)
}

// Lay the sources' cells out end to end in the combined count matrix, and
// reject any cell that appears in more than one file.
//
// A Bacha score is a non-linear function of a cell's insertion counts, so a
// cell whose fragments are split across two files cannot be scored from either
// file alone
//
// Returns (global cell order, per-source start offsets with a final total).
fn global_cell_layout(sources: &[FragmentSource]) -> io::Result<(Vec<String>, Vec<usize>)> {
    let total: usize = sources.iter().map(|s| s.cell_names.len()).sum();
    let mut cell_order: Vec<String> = Vec::with_capacity(total);
    let mut offsets: Vec<usize> = Vec::with_capacity(sources.len() + 1);
    let mut owner: FxHashMap<&str, &str> = FxHashMap::default();

    for source in sources {
        offsets.push(cell_order.len());
        for name in &source.cell_names {
            if let Some(prev) = owner.insert(name.as_str(), source.sample.as_str()) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "cell '{}' is listed for both sample '{}' and sample '{}'. Bacha \
                         scores each fragment file separately, and a score is a non-linear \
                         function of the insertion counts, so per-file counts for one cell \
                         cannot be combined into a whole-cell score. Give the cells distinct \
                         output names if they are different cells, or merge the fragment \
                         files into a single sorted, bgzip-compressed file if they are not.",
                        name, prev, source.sample
                    ),
                ));
            }
            cell_order.push(name.clone());
        }
    }
    offsets.push(cell_order.len());
    Ok((cell_order, offsets))
}

struct IndexTables {
    lapper_map: FxHashMap<String, Lapper<u32, u32>>,
    site_indices: Vec<u32>,
    control_indices: FlatIndVec,
    num_site_intervals: usize,
    tf_list: Vec<String>,
    n_instances: Vec<u32>,
    chrom_names: Vec<String>,
    chrom_to_id: FxHashMap<String, u8>,
    chrom_size_by_id: Vec<u32>,
    n_features: usize,
}

fn load_index_tables(
    index_file: &Path,
    num_threads: usize,
    check_bounds: bool,
) -> io::Result<IndexTables> {
    eprintln!("Reading index file");

    let t_idx = Instant::now();
    let (
        lapper_map,
        site_indices,
        control_indices,
        num_site_intervals,
        tf_list,
        n_instances,
        chrom_sizes,
    ) = match parse_file_to_lapper(index_file, num_threads) {
        Ok(trees) => trees,
        Err(e) => {
            error!("Failed to read index file: {}", e);
            return Err(e);
        }
    };
    if std::env::var("BACHA_PROFILE").is_ok() {
        eprintln!(
            "[PERF] parse_file_to_lapper={:.2}s",
            t_idx.elapsed().as_secs_f64()
        );
    }

    let n_features: usize = tf_list.len() * 2;

    // Create chromosome name to ID mapping.
    //
    // Every chromosome named in the index's chromosome table gets an ID
    // not just the ones that carry intervals
    let mut chrom_names: Vec<String> = lapper_map.keys().cloned().collect();
    {
        let with_intervals: HashSet<&str> = chrom_names.iter().map(|s| s.as_str()).collect();
        let mut interval_free: Vec<String> = chrom_sizes
            .keys()
            .filter(|name| !with_intervals.contains(name.as_str()))
            .cloned()
            .collect();
        interval_free.sort();
        chrom_names.extend(interval_free);
    }
    if chrom_names.len() > u8::MAX as usize + 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Index contains {} chromosomes; the maximum supported is {}.",
                chrom_names.len(),
                u8::MAX as usize + 1
            ),
        ));
    }
    let chrom_to_id: FxHashMap<String, u8> = chrom_names
        .iter()
        .enumerate()
        .map(|(i, name)| (name.clone(), i as u8))
        .collect();

    // Build chromosome size lookup by chrom_id for the genome/assembly bounds
    // check. A length of u32::MAX means "no bound": the index stored no size
    // for this chromosome, so no coordinate can exceed it. When the caller
    // disables the check (--no-bounds-check) we set every length to u32::MAX,
    // which neutralises the inline comparison without branching in the hot loop.
    let chrom_size_by_id: Vec<u32> = chrom_names
        .iter()
        .map(|name| {
            if check_bounds {
                chrom_sizes.get(name).copied().unwrap_or(u32::MAX)
            } else {
                u32::MAX
            }
        })
        .collect();

    Ok(IndexTables {
        lapper_map,
        site_indices,
        control_indices,
        num_site_intervals,
        tf_list,
        n_instances,
        chrom_names,
        chrom_to_id,
        chrom_size_by_id,
        n_features,
    })
}

// Count every source into one combined matrix.
//
// Fragment files are processed one at a time, each using all `num_threads`
// workers. The index is loaded once up front, and each source writes into its
// own contiguous block of rows, so no per-file matrices are allocated and no
// concatenation step is needed: the matrix this returns is already the final
// one.
fn count_sources(
    sources: &[FragmentSource],
    index_file: &Path,
    num_threads: usize,
    check_bounds: bool,
) -> io::Result<(Array2<u32>, Vec<String>, IndexTables)> {
    let multi = sources.len() > 1;
    // Named for messages only when there is more than one file to confuse.
    fn label<'a>(multi: bool, s: &'a FragmentSource) -> Option<&'a str> {
        if multi {
            Some(s.sample.as_str())
        } else {
            None
        }
    }

    // Before the (expensive) index load, so a truncated fragment file is
    // rejected immediately
    //
    // `File::open`'s io::Error carries no path, and a manifest's row paths are
    // resolved against the manifest's own directory, so an error without the
    // path it tried leaves the user unable to tell which file bacha looked for.
    // The -f/-c route gets this from `canonicalize_input`; add it here.
    for source in sources {
        check_fragment_file_complete(&source.path)
            .map_err(|e| {
                io::Error::new(
                    e.kind(),
                    format!(
                        "cannot access fragment file '{}': {e}",
                        source.path.display()
                    ),
                )
            })
            .map_err(|e| with_sample_context(label(multi, source), e))?;
    }

    let (cell_order, offsets) = global_cell_layout(sources)?;
    let tables = load_index_tables(index_file, num_threads, check_bounds)?;

    // One matrix for the whole run, laid out as (n_cells, n_features)
    // row-major. Each cell's counts live contiguously, which makes the
    // multi-feature-per-control write fan-out cache-friendly.
    let mut count_matrix: Array2<u32> = Array2::zeros((cell_order.len(), tables.n_features));

    let t_count = Instant::now();
    for (i, source) in sources.iter().enumerate() {
        if multi {
            eprintln!(
                "[{}/{}] {}: {}",
                i + 1,
                sources.len(),
                source.sample,
                source.path.display()
            );
        }
        let stats = {
            let mut block = count_matrix.slice_mut(ndarray::s![offsets[i]..offsets[i + 1], ..]);
            count_one_file(
                &tables,
                &source.path,
                &source.cells,
                &mut block,
                num_threads,
            )
            .map_err(|e| with_sample_context(label(multi, source), e))?
        };

        report_skipped_lines(
            stats.malformed,
            stats.total_lines,
            stats.unterminated_final_line,
        );
        check_input_overlap(
            stats.total_fragments,
            stats.matched_cells,
            stats.unmatched_chroms,
        )
        .map_err(|e| with_sample_context(label(multi, source), e))?;
    }
    if std::env::var("BACHA_PROFILE").is_ok() {
        eprintln!(
            "[PERF] counting={:.2}s files={} cells={}",
            t_count.elapsed().as_secs_f64(),
            sources.len(),
            cell_order.len(),
        );
    }

    Ok((count_matrix, cell_order, tables))
}

// Count every source, then score and write the combined matrix once.
fn run_count(
    sources: &[FragmentSource],
    index_file: &Path,
    output_path: &Path,
    num_threads: usize,
    lambda: f32,
    pseudocount: f32,
    check_bounds: bool,
) -> io::Result<()> {
    let (count_matrix, cell_order, tables) =
        count_sources(sources, index_file, num_threads, check_bounds)?;

    let t1 = Instant::now();
    compute_and_write_scores(
        &count_matrix,
        &tables.n_instances,
        output_path,
        b',',
        &cell_order,
        &tables.tf_list,
        num_threads,
        lambda,
        pseudocount,
    )?;
    if std::env::var("BACHA_PROFILE").is_ok() {
        eprintln!(
            "[PERF] compute_and_write_scores={:.2}s",
            t1.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

// Count a single fragment file and hand back the raw matrix.
#[cfg(test)]
fn fcount_core(
    frag_file: &Path,
    index_file: &Path,
    cell_file: &Path,
    num_threads: usize,
    check_bounds: bool,
) -> io::Result<(Array2<u32>, Vec<String>, Vec<String>, Vec<u32>)> {
    let (cells, cell_names) = read_cells_file(cell_file)?;
    let sources = vec![FragmentSource {
        sample: "test".to_string(),
        path: frag_file.to_path_buf(),
        cells,
        cell_names,
    }];
    let (count_matrix, cell_order, tables) =
        count_sources(&sources, index_file, num_threads, check_bounds)?;
    let IndexTables {
        tf_list,
        n_instances,
        ..
    } = tables;
    Ok((count_matrix, cell_order, tf_list, n_instances))
}

// Count one fragment file into `count_matrix`, which is the block of rows
// belonging to this file (the whole matrix for a single-file run). `cells` maps
// a barcode to a row index within that block.
//
// Reports the fragment/line tallies rather than acting on them, so the caller
// can run the genome and barcode sanity checks per file and name the sample.
fn count_one_file(
    tables: &IndexTables,
    frag_file: &Path,
    cells: &FxHashMap<String, u32>,
    count_matrix: &mut ndarray::ArrayViewMut2<u32>,
    num_threads: usize,
) -> io::Result<FileStats> {
    let lapper_map = &tables.lapper_map;
    let site_indices: &[u32] = &tables.site_indices;
    let control_indices: &FlatIndVec = &tables.control_indices;
    let num_site_intervals = tables.num_site_intervals;
    let chrom_names = &tables.chrom_names;
    let chrom_to_id = &tables.chrom_to_id;
    let chrom_size_by_id: &[u32] = &tables.chrom_size_by_id;
    let n_features = tables.n_features;
    let n_cells = count_matrix.nrows();

    // PERF instrumentation prints a handful of timing summaries to stderr
    // when set; off by default to keep normal runs quiet. The timers
    // themselves are always collected (negligible overhead).
    let perf = std::env::var("BACHA_PROFILE").is_ok();

    // Split --threads between BGZF inflaters and counting workers; see
    // `thread_budget`. Workers do all line parsing, cell-barcode filtering,
    // chromosome lookup, coordinate parsing, lapper overlap queries and matrix
    // updates; the reader thread is reduced to byte streaming plus line-boundary
    // chunking. Decided per fragment file, so a manifest may mix BGZF and plain
    // gzip inputs and each gets the right split.
    let budget = thread_budget(num_threads, is_bgzf_file(frag_file)?);
    let n_workers = budget.counting;
    let bgzf_workers = budget.bgzf;
    if perf {
        eprintln!(
            "[PERF] thread budget: {} total = {} bgzf + {} counting",
            num_threads.max(1),
            bgzf_workers,
            n_workers,
        );
    }

    // Bounded MPMC channel of raw byte chunks (each ~CHUNK_BYTES of complete
    // lines). The reader is the only sender; workers each hold a clone of the
    // receiver. When the reader's closure returns, its sender is dropped and
    // workers' rx.recv() returns Err, ending their loops.
    const CHUNK_BYTES: usize = 1 << 20; // ~1 MB per chunk
    let (tx, rx) = bounded::<Vec<u8>>(64);

    // Single-threaded fast path. With -t 1 the user is asking for one core;
    // honour that by skipping thread::scope, the MPMC channel, and the
    // noodles BGZF parallel reader entirely. MultiGzDecoder handles BGZF
    // transparently (BGZF is just concatenated gzip blocks) so we don't
    // need to fork on the file format here.
    if n_workers == 1 {
        let file = match File::open(frag_file) {
            Ok(f) => f,
            Err(e) => {
                error!("Failed to open fragment file: {}", e);
                return Err(e);
            }
        };
        eprintln!("Using single-threaded gzip decoder (-t 1)");
        let mut reader = BufReader::with_capacity(8 * 1024 * 1024, MultiGzDecoder::new(file));

        let mut line_buf: Vec<u8> = Vec::with_capacity(256);
        let mut current_chrom_id: Option<u8> = None;
        let mut current_lapper: Option<&Lapper<u32, u32>> = None;
        let mut total_fragments: u64 = 0;
        let mut matched_cells: u64 = 0;
        let mut unmatched_chroms: u64 = 0;
        let mut malformed: u64 = 0;
        let mut total_lines: u64 = 0;
        let mut unterminated_final_line = false;
        let t0 = Instant::now();

        loop {
            line_buf.clear();
            let n = match reader.read_until(b'\n', &mut line_buf) {
                Ok(n) => n,
                Err(e) => {
                    return Err(io::Error::new(
                        e.kind(),
                        format!(
                            "Error reading fragment file after {total_fragments} \
                             fragment(s); the file is truncated or unreadable: {e}"
                        ),
                    ));
                }
            };
            if n == 0 {
                break;
            }

            let mut line: &[u8] = &line_buf;
            if line.last() == Some(&b'\n') {
                line = &line[..line.len() - 1];
            } else {
                // Only the last line of the file can lack its terminator.
                unterminated_final_line = true;
            }
            if line.last() == Some(&b'\r') {
                line = &line[..line.len() - 1];
            }
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            total_lines += 1;

            let fields: SmallVec<[&[u8]; 10]> = line.split(|&b| b == b'\t').collect();
            if fields.len() < 4 {
                malformed += 1;
                continue;
            }

            total_fragments += 1;
            if total_fragments % 1_000_000 == 0 {
                // Progress goes to stderr, keeping stdout free for `bacha
                // view`/`summary` output and for redirection.
                eprint!("\rProcessed {} M fragments", total_fragments / 1_000_000);
                io::stderr().flush().ok();
            }

            let cell_barcode = match std::str::from_utf8(fields[3]) {
                Ok(s) => s,
                Err(_) => {
                    malformed += 1;
                    continue;
                }
            };
            let cell_idx = match cells.get(cell_barcode) {
                Some(&i) => i,
                None => continue,
            };
            matched_cells += 1;

            let seqname = match std::str::from_utf8(fields[0]) {
                Ok(s) => s,
                Err(_) => {
                    malformed += 1;
                    continue;
                }
            };
            let chrom_id = match chrom_to_id.get(seqname) {
                Some(&i) => i,
                None => {
                    unmatched_chroms += 1;
                    continue;
                }
            };

            let start = match parse::<u32>(fields[1]) {
                Ok(n) => n,
                Err(_) => {
                    malformed += 1;
                    continue;
                }
            };
            let end = match parse::<u32>(fields[2]) {
                Ok(n) => n,
                Err(_) => {
                    malformed += 1;
                    continue;
                }
            };

            // Fail fast on a genome/assembly mismatch: a coordinate beyond the
            // chromosome's length in the index is impossible on a matching
            // assembly, so abort immediately rather than reading the rest of
            // the file. (chrom length is u32::MAX for size-less indexes.)
            let chrom_len = chrom_size_by_id[chrom_id as usize];
            if end > chrom_len {
                return Err(assembly_mismatch_error(
                    &chrom_names[chrom_id as usize],
                    end,
                    chrom_len,
                ));
            }

            if current_chrom_id != Some(chrom_id) {
                current_chrom_id = Some(chrom_id);
                current_lapper = lapper_map.get(&chrom_names[chrom_id as usize]);
            }
            if let Some(lapper) = current_lapper {
                let row = count_matrix
                    .row_mut(cell_idx as usize)
                    .into_slice()
                    .expect("count_matrix rows are contiguous (row-major Array2)");
                // BED end is exclusive: the second insertion is at end - 1.
                // saturating_sub guards a malformed end of 0.
                process_fragment(
                    start,
                    end.saturating_sub(1),
                    row,
                    lapper,
                    site_indices,
                    control_indices,
                    num_site_intervals,
                    n_features,
                );
            }
        }

        eprintln!("\nFinished reading {} total fragments", total_fragments);
        if perf {
            eprintln!(
                "[PERF] single-threaded count: wall={:.2}s matched_cells={} ({:.1}%)",
                t0.elapsed().as_secs_f64(),
                matched_cells,
                100.0 * matched_cells as f64 / total_fragments.max(1) as f64,
            );
        }

        return Ok(FileStats {
            total_fragments,
            matched_cells,
            unmatched_chroms,
            malformed,
            total_lines,
            unterminated_final_line,
        });
    }

    // Shared atomic for cross-thread progress reporting. Workers increment in
    // 64K-line batches; the reader thread (mostly idle waiting on BGZF I/O)
    // polls and prints. Atomic contention is negligible at this cadence.
    let progress = AtomicU64::new(0);

    // Genome/assembly mismatch abort. When any worker sees a coordinate beyond
    // its chromosome's length in the index, it records the offending
    // observation in `mismatch` and flips `abort`; the reader and all workers
    // check `abort` and wind down promptly so we fail on the first bad
    // coordinate rather than reading the whole file. The error is built from
    // `mismatch` after the scope joins.
    let abort = AtomicBool::new(false);
    let mismatch: Mutex<Option<(String, u32, u32)>> = Mutex::new(None);

    // Sharded parallel path. Each worker owns a contiguous row range of the
    // shared `count_matrix` (no per-worker private shadow matrix). Workers
    // parse chunks from the global queue and either count in-shard fragments
    // directly into their slice, or forward parsed coordinates to the
    // owning worker via per-worker unbounded channels.
    let shard_size = (n_cells + n_workers - 1) / n_workers.max(1);
    let mut shards: Vec<ndarray::ArrayViewMut2<u32>> = Vec::with_capacity(n_workers);
    {
        let mut remaining = count_matrix.view_mut();
        for w in 0..n_workers {
            let nrows_remaining = remaining.nrows();
            let take = if w == n_workers - 1 {
                nrows_remaining
            } else {
                shard_size.min(nrows_remaining)
            };
            let (head, rest) = remaining.split_at(ndarray::Axis(0), take);
            shards.push(head);
            remaining = rest;
        }
    }

    // Per-worker forwarding channels. `unbounded` avoids deadlock in the
    // cyclic forwarding graph (every worker can send to every peer); the
    // total in-flight footprint is bounded by the upstream chunk queue
    // (~64 MB) and runs to a few tens of MB of records in practice.
    let (forward_txs, forward_rxs): (
        Vec<Sender<Vec<ForwardRecord>>>,
        Vec<Receiver<Vec<ForwardRecord>>>,
    ) = (0..n_workers).map(|_| unbounded()).unzip();

    let (
        total_fragments,
        matched_cells,
        unmatched_chroms,
        malformed,
        total_lines,
        read_err,
        unterminated_final_line,
    ) = thread::scope(|scope| {
        // ---- Reader thread: just bytes, no parsing ----
        let progress_ref = &progress;
        let abort_ref = &abort;
        // Returns any I/O failure so the caller can fail the run instead of
        // scoring the prefix of the file that happened to be readable, plus
        // whether the stream ended mid-line (see `report_skipped_lines`).
        let reader_handle = scope.spawn(move || -> (Option<io::Error>, bool) {
                let mut reader = match open_fragment_reader(frag_file, bgzf_workers) {
                    Ok(r) => r,
                    Err(e) => return (Some(e), false),
                };

                let mut chunk: Vec<u8> = Vec::with_capacity(CHUNK_BYTES + 65536);
                let mut tmp = [0u8; 65536];
                let reader_t0 = Instant::now();
                let mut t_io_total = std::time::Duration::ZERO;
                let mut send_block_total = std::time::Duration::ZERO;
                let mut last_reported_m: u64 = 0;
                let mut eof = false;

                while !eof {
                    // Stop producing if a worker hit a genome mismatch.
                    if abort_ref.load(Ordering::Relaxed) {
                        return (None, false);
                    }
                    // Fill `chunk` to at least CHUNK_BYTES (carrying any
                    // trailing partial line from the previous iteration).
                    while chunk.len() < CHUNK_BYTES {
                        let t = Instant::now();
                        let n = match reader.read(&mut tmp) {
                            Ok(n) => n,
                            Err(e) => return (Some(e), false),
                        };
                        t_io_total += t.elapsed();
                        if n == 0 {
                            eof = true;
                            break;
                        }
                        chunk.extend_from_slice(&tmp[..n]);
                    }

                    if chunk.is_empty() {
                        break;
                    }

                    // Slice the chunk at the last newline so each shipped chunk
                    // contains only complete lines; the tail becomes the carry.
                    let last_nl = chunk.iter().rposition(|&b| b == b'\n');
                    let to_send = if let Some(idx) = last_nl {
                        let mut out = Vec::with_capacity(idx + 1);
                        out.extend_from_slice(&chunk[..=idx]);
                        // Move the tail back to the start of `chunk` for next round.
                        chunk.drain(..=idx);
                        out
                    } else if eof {
                        // No newline at all and EOF: ship whatever's left.
                        std::mem::take(&mut chunk)
                    } else {
                        // No newline yet; keep accumulating.
                        continue;
                    };

                    let t = Instant::now();
                    let send_res = tx.send(to_send);
                    send_block_total += t.elapsed();
                    if send_res.is_err() {
                        return (None, false);
                    }

                    // Progress check (cheap; runs once per ~13K lines).
                    let now = progress_ref.load(Ordering::Relaxed);
                    let m = now / 1_000_000;
                    if m > last_reported_m {
                        last_reported_m = m;
                        eprint!("\rProcessed {} M fragments", m);
                        io::stderr().flush().ok();
                    }
                }

                let unterminated_final_line = !chunk.is_empty();
                if unterminated_final_line && tx.send(std::mem::take(&mut chunk)).is_err() {
                    return (None, true);
                }

                // Drop tx by returning; workers will see disconnect.
                let reader_wall = reader_t0.elapsed();
                if perf {
                    eprintln!(
                        "\n[PERF] reader thread: wall={:.2}s io_read={:.2}s blocked_on_send={:.2}s ({:.1}%)",
                        reader_wall.as_secs_f64(),
                        t_io_total.as_secs_f64(),
                        send_block_total.as_secs_f64(),
                        100.0 * send_block_total.as_secs_f64() / reader_wall.as_secs_f64().max(1e-9),
                    );
                } else {
                    // Make sure the progress line ends with a newline even
                    // when PERF output is suppressed.
                    eprintln!();
                }
                (None, unterminated_final_line)
            });

        // ---- Worker threads: own a shard, parse chunks, route by cell ----
        let lapper_map_ref = lapper_map;
        let chrom_names_ref = chrom_names;
        let site_indices_ref = site_indices;
        let control_indices_ref = control_indices;
        let cells_ref = cells;
        let chrom_to_id_ref = chrom_to_id;
        let chrom_size_ref = chrom_size_by_id;
        let forward_txs_ref = &forward_txs;
        let abort_ref = &abort;
        let mismatch_ref = &mismatch;

        let consumer_t0 = Instant::now();
        let worker_handles: Vec<_> = shards
            .into_iter()
            .enumerate()
            .map(|(wid, mut shard)| {
                let chunk_rx = rx.clone();
                let inbox_rx = forward_rxs[wid].clone();
                let peer_txs: Vec<Sender<Vec<ForwardRecord>>> = forward_txs_ref.clone();
                scope.spawn(move || {
                    let mut local_total: u64 = 0;
                    let mut local_matched: u64 = 0;
                    let mut local_unmatched_chroms: u64 = 0;
                    let mut local_malformed: u64 = 0;
                    let mut local_lines: u64 = 0;
                    let mut since_progress_flush: u64 = 0;
                    let mut recv_block_total = std::time::Duration::ZERO;
                    let worker_t0 = Instant::now();
                    let mut cached_chrom_id: Option<u8> = None;
                    let mut cached_lapper: Option<&Lapper<u32, u32>> = None;

                    // Write one fragment's two endpoint counts into the
                    // local shard. Caches the lapper per chrom across calls
                    // (sorted fragment input makes the hit rate ~100%; the
                    // forwarded-batch drain hits less but the lookup is
                    // a cheap FxHashMap get on miss).
                    let mut apply =
                        |chrom_id: u8,
                         start: u32,
                         end: u32,
                         local_idx: usize,
                         shard: &mut ndarray::ArrayViewMut2<u32>| {
                            if cached_chrom_id != Some(chrom_id) {
                                cached_chrom_id = Some(chrom_id);
                                cached_lapper =
                                    lapper_map_ref.get(&chrom_names_ref[chrom_id as usize]);
                            }
                            if let Some(lapper) = cached_lapper {
                                let row = shard
                                    .row_mut(local_idx)
                                    .into_slice()
                                    .expect("shard rows are contiguous (row-major)");
                                process_fragment(
                                    start,
                                    end.saturating_sub(1),
                                    row,
                                    lapper,
                                    site_indices_ref,
                                    control_indices_ref,
                                    num_site_intervals,
                                    n_features,
                                );
                            }
                        };

                    loop {
                        // Bail out if another worker hit a genome mismatch.
                        if abort_ref.load(Ordering::Relaxed) {
                            break;
                        }
                        let t = Instant::now();
                        let chunk = match chunk_rx.recv() {
                            Ok(c) => {
                                recv_block_total += t.elapsed();
                                c
                            }
                            Err(_) => {
                                recv_block_total += t.elapsed();
                                break;
                            }
                        };

                        // Per-peer outgoing batches. Filled while parsing
                        // this chunk, flushed once at end. This collapses
                        // millions of single-record sends into a handful
                        // of batched ones, and keeps the channel hot.
                        let mut outgoing: Vec<Vec<ForwardRecord>> =
                            (0..n_workers).map(|_| Vec::new()).collect();

                        for line in chunk.split(|&b| b == b'\n') {
                            if line.is_empty() {
                                continue;
                            }
                            let line: &[u8] = if line.last() == Some(&b'\r') {
                                &line[..line.len() - 1]
                            } else {
                                line
                            };
                            if line.is_empty() || line[0] == b'#' {
                                continue;
                            }
                            local_lines += 1;

                            let fields: SmallVec<[&[u8]; 10]> =
                                line.split(|&b| b == b'\t').collect();
                            if fields.len() < 4 {
                                local_malformed += 1;
                                continue;
                            }

                            local_total += 1;
                            since_progress_flush += 1;
                            if since_progress_flush >= 65_536 {
                                progress_ref.fetch_add(since_progress_flush, Ordering::Relaxed);
                                since_progress_flush = 0;
                            }

                            let cell_barcode = match std::str::from_utf8(fields[3]) {
                                Ok(s) => s,
                                Err(_) => {
                                    local_malformed += 1;
                                    continue;
                                }
                            };
                            let cell_idx = match cells_ref.get(cell_barcode) {
                                Some(&i) => i,
                                None => continue,
                            };
                            local_matched += 1;

                            let seqname = match std::str::from_utf8(fields[0]) {
                                Ok(s) => s,
                                Err(_) => {
                                    local_malformed += 1;
                                    continue;
                                }
                            };
                            let chrom_id = match chrom_to_id_ref.get(seqname) {
                                Some(&i) => i,
                                None => {
                                    local_unmatched_chroms += 1;
                                    continue;
                                }
                            };

                            let start = match parse::<u32>(fields[1]) {
                                Ok(n) => n,
                                Err(_) => {
                                    local_malformed += 1;
                                    continue;
                                }
                            };
                            let end = match parse::<u32>(fields[2]) {
                                Ok(n) => n,
                                Err(_) => {
                                    local_malformed += 1;
                                    continue;
                                }
                            };

                            // Fail fast on a genome/assembly mismatch: a
                            // coordinate beyond the chromosome's length in
                            // the index is impossible on a matching
                            // assembly. Record the first offender and set
                            // `abort` so every thread winds down without
                            // reading the rest of the file. (chrom length is
                            // u32::MAX for size-less indexes, never tripping.)
                            let chrom_len = chrom_size_ref[chrom_id as usize];
                            if end > chrom_len {
                                if !abort_ref.swap(true, Ordering::Relaxed) {
                                    *mismatch_ref.lock().unwrap() = Some((
                                        chrom_names_ref[chrom_id as usize].clone(),
                                        end,
                                        chrom_len,
                                    ));
                                }
                                break;
                            }

                            // Route by cell shard. Contiguous row ranges
                            // mean `owner = cell_idx / shard_size`,
                            // clamped to the last shard for any straggler
                            // index past the last-shard boundary.
                            let owner = ((cell_idx as usize) / shard_size).min(n_workers - 1);
                            let local_idx = cell_idx as usize - owner * shard_size;

                            if owner == wid {
                                apply(chrom_id, start, end, local_idx, &mut shard);
                            } else {
                                outgoing[owner].push(ForwardRecord {
                                    local_idx: local_idx as u32,
                                    chrom_id,
                                    start,
                                    end,
                                });
                            }
                        }

                        // Flush outgoing batches.
                        for (peer, batch) in outgoing.into_iter().enumerate() {
                            if peer != wid && !batch.is_empty() {
                                // Unbounded channel; send only fails if all
                                // receivers were dropped (worker scope teardown).
                                let _ = peer_txs[peer].send(batch);
                            }
                        }

                        // Opportunistic non-blocking drain of own inbox.
                        // Keeps memory pressure low and balances latency
                        // with the upstream chunk producer.
                        while let Ok(batch) = inbox_rx.try_recv() {
                            for rec in &batch {
                                apply(
                                    rec.chrom_id,
                                    rec.start,
                                    rec.end,
                                    rec.local_idx as usize,
                                    &mut shard,
                                );
                            }
                        }
                    }

                    // We've stopped consuming chunks: either the reader
                    // finished and the channel closed, or a genome mismatch
                    // set `abort`. Drop our chunk receiver so that once every
                    // worker does the same the reader sees the disconnect and
                    // stops producing (it may be blocked on a full queue).
                    drop(chunk_rx);

                    // On a genome-mismatch abort there is no point draining
                    // forwarded records into a result we're about to discard.
                    // Drop peer_txs and bail so the other workers' blocking
                    // drains (and this one, in peers) can unblock.
                    if abort_ref.load(Ordering::Relaxed) {
                        drop(peer_txs);
                        return (
                            worker_t0.elapsed(),
                            recv_block_total,
                            local_total,
                            local_matched,
                            local_unmatched_chroms,
                            local_malformed,
                            local_lines,
                        );
                    }

                    // Drop our peer_txs so the other workers' inboxes can
                    // eventually close (each inbox closes only when ALL its
                    // senders have dropped).
                    drop(peer_txs);

                    // Blocking drain until every peer has dropped its
                    // tx-to-me. The only senders left at this point are
                    // peers still inside their own chunk loop.
                    loop {
                        let t = Instant::now();
                        let batch = match inbox_rx.recv() {
                            Ok(b) => {
                                recv_block_total += t.elapsed();
                                b
                            }
                            Err(_) => {
                                recv_block_total += t.elapsed();
                                break;
                            }
                        };
                        for rec in &batch {
                            apply(
                                rec.chrom_id,
                                rec.start,
                                rec.end,
                                rec.local_idx as usize,
                                &mut shard,
                            );
                        }
                    }

                    // Flush any remaining progress counter.
                    if since_progress_flush > 0 {
                        progress_ref.fetch_add(since_progress_flush, Ordering::Relaxed);
                    }

                    (
                        worker_t0.elapsed(),
                        recv_block_total,
                        local_total,
                        local_matched,
                        local_unmatched_chroms,
                        local_malformed,
                        local_lines,
                    )
                })
            })
            .collect();

        // Drop our copy of the chunk receiver and the forwarding senders
        // so dependants see disconnect once the reader and workers finish.
        drop(rx);
        drop(forward_txs);
        drop(forward_rxs);

        // Reduce stats. The matrix is already populated in-place via the
        // workers' shard views, so no matrix reduction is needed.
        let mut worker_wall_max = std::time::Duration::ZERO;
        let mut worker_active_total = std::time::Duration::ZERO;
        let mut worker_idle_total = std::time::Duration::ZERO;
        let mut total: u64 = 0;
        let mut matched: u64 = 0;
        let mut unmatched_chroms: u64 = 0;
        let mut malformed: u64 = 0;
        let mut total_lines: u64 = 0;
        for h in worker_handles {
            let (wall, idle, t, m, uc, mf, nl) = h.join().expect("Worker thread panicked");
            if wall > worker_wall_max {
                worker_wall_max = wall;
            }
            worker_active_total += wall.saturating_sub(idle);
            worker_idle_total += idle;
            total += t;
            matched += m;
            unmatched_chroms += uc;
            malformed += mf;
            total_lines += nl;
        }
        let consumer_wall = consumer_t0.elapsed();
        if perf {
            eprintln!(
                    "[PERF] consumer pool: workers={} wall={:.2}s worker_wall_max={:.2}s active_sum={:.2}s idle_sum={:.2}s (effective_parallelism={:.2}x)",
                    n_workers,
                    consumer_wall.as_secs_f64(),
                    worker_wall_max.as_secs_f64(),
                    worker_active_total.as_secs_f64(),
                    worker_idle_total.as_secs_f64(),
                    worker_active_total.as_secs_f64() / consumer_wall.as_secs_f64().max(1e-9),
                );
            eprintln!(
                "Finished reading {} total fragments (matched_cells={} = {:.1}%)",
                total,
                matched,
                100.0 * matched as f64 / total.max(1) as f64,
            );
        } else {
            eprintln!("Finished reading {} total fragments", total);
        }

        let (read_err, unterminated_final_line) =
            reader_handle.join().expect("Reader thread panicked");
        (
            total,
            matched,
            unmatched_chroms,
            malformed,
            total_lines,
            read_err,
            unterminated_final_line,
        )
    });

    // A worker aborted on an out-of-bounds coordinate: surface the genome
    // mismatch. This is checked before the unmatched-chromosome heuristics
    // below because it is a definitive assembly mismatch, not a warning.
    if let Some((chrom, observed, index_len)) = mismatch.into_inner().unwrap() {
        return Err(assembly_mismatch_error(&chrom, observed, index_len));
    }

    // The reader failed
    if let Some(e) = read_err {
        if total_fragments == 0 {
            return Err(e);
        }
        return Err(io::Error::new(
            e.kind(),
            format!(
                "Error reading fragment file after {total_fragments} fragment(s); \
                 the file is truncated or unreadable: {e}"
            ),
        ));
    }

    Ok(FileStats {
        total_fragments,
        matched_cells,
        unmatched_chroms,
        malformed,
        total_lines,
        unterminated_final_line,
    })
}

pub fn compute_and_write_scores(
    count_matrix: &Array2<u32>,
    n_instances: &Vec<u32>,
    path: &Path,
    delimiter: u8,
    cell_order: &[String],
    tf_list: &[String],
    num_threads: usize,
    lambda: f32,
    pseudocount: f32,
) -> io::Result<()> {
    validate_thread_count(num_threads)?;
    validate_scoring_parameters(lambda, pseudocount)?;

    let nsites = n_instances.len() / 2;

    // Pre-calculate site ratios. None marks a TF that cannot be scored
    // (no sites and/or no controls); those rows are emitted as 0.0.
    let site_ratios: Vec<Option<f32>> = (0..nsites)
        .map(|i| {
            let n_sites = n_instances[i] as f32;
            let n_controls = n_instances[nsites + i] as f32;
            if n_sites > 0.0 && n_controls > 0.0 {
                Some(n_controls / n_sites)
            } else {
                None
            }
        })
        .collect();

    let unscoreable: Vec<&String> = site_ratios
        .iter()
        .zip(tf_list.iter())
        .filter_map(|(r, name)| if r.is_none() { Some(name) } else { None })
        .collect();
    if !unscoreable.is_empty() {
        eprintln!(
            "WARNING: {} feature(s) have zero sites or zero controls and will be reported as 0.0: {}",
            unscoreable.len(),
            unscoreable.iter().take(10).map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
        );
    }

    // Set up the output writer. CSV escaping is handled inline by
    // `write_csv_field` below; score columns are always plain "%.3f"
    // floats and never need quoting, so we can write raw bytes for the
    // bulk of each row.
    let file = File::create(path)?;
    let mut writer: Box<dyn io::Write> = if path.extension().unwrap_or_default() == "gz" {
        Box::new(
            ParCompressBuilder::<Gzip>::new()
                .compression_level(Compression::default())
                .num_threads(num_threads)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
                .from_writer(file),
        )
    } else {
        Box::new(file)
    };

    let ncells = count_matrix.nrows();

    // Header: leading empty cell, then one column per cell barcode.
    let mut header_buf: Vec<u8> = Vec::with_capacity(ncells * 20 + 2);
    for cell in cell_order {
        header_buf.push(delimiter);
        write_csv_field(&mut header_buf, cell.as_bytes(), delimiter);
    }
    header_buf.push(b'\n');
    writer.write_all(&header_buf)?;

    // Precision-based Bayesian shrinkage weight, consistent with the Gamma-
    // Poisson prior implied by the pseudocounts in the LFC formula.
    //
    // The LFC is log2((s + p) / (c/r + p)), which corresponds to a Gamma(p, .)
    // prior on each per-bp rate. Because the control window is r-fold larger
    // than the site window, p pseudo-observations on the site footprint equal
    // p*r pseudo-observations on the control footprint. The posterior
    // asymptotic variance of the LFC is then
    //     Var(LFC) ~ 1/(s + p) + 1/(c + p*r),
    // giving data precision
    //     prec = (s + p)(c + p*r) / (s + c + p + p*r).
    // With a Normal prior of variance sigma^2 on log-enrichment, the posterior
    // shrinkage weight is
    //     w = prec / (prec + 1/sigma^2)
    // and we parameterise the prior precision as `lambda` directly.

    // Build feature rows in parallel chunks. Buffering every row at once
    // would cost ~n_features × n_cells × 8 bytes (~500 MB for 800 TFs and
    // 100k cells); a chunk of 32 rows caps the buffered footprint at a
    // few tens of MB while still keeping all cores busy.
    const CHUNK_SIZE: usize = 32;
    for chunk_start in (0..nsites).step_by(CHUNK_SIZE) {
        let chunk_end = (chunk_start + CHUNK_SIZE).min(nsites);
        let rows: Vec<Vec<u8>> = (chunk_start..chunk_end)
            .into_par_iter()
            .map(|i| {
                build_score_row(
                    i,
                    nsites,
                    ncells,
                    count_matrix,
                    site_ratios[i],
                    tf_list[i].as_bytes(),
                    delimiter,
                    lambda,
                    pseudocount,
                )
            })
            .collect();
        for row in rows {
            writer.write_all(&row)?;
        }
    }

    writer.flush()?;
    Ok(())
}

// Append `field` to `out` with minimal CSV escaping: quote and double up
// embedded quotes only if the field contains the delimiter, a quote, or
// a line terminator. Score values bypass this and write raw bytes.
fn write_csv_field(out: &mut Vec<u8>, field: &[u8], delimiter: u8) {
    let needs_quote = field
        .iter()
        .any(|&b| b == delimiter || b == b'"' || b == b'\n' || b == b'\r');
    if needs_quote {
        out.push(b'"');
        for &b in field {
            if b == b'"' {
                out.extend_from_slice(b"\"\"");
            } else {
                out.push(b);
            }
        }
        out.push(b'"');
    } else {
        out.extend_from_slice(field);
    }
}

fn build_score_row(
    i: usize,
    nsites: usize,
    ncells: usize,
    count_matrix: &Array2<u32>,
    site_ratio: Option<f32>,
    tf_name: &[u8],
    delimiter: u8,
    lambda: f32,
    pseudocount: f32,
) -> Vec<u8> {
    use std::fmt::Write as _;
    // Reserve ~8 bytes per value (e.g. "-1.234,"); resizing rarely matters
    // because most rows fall well under this once formatted.
    let mut buf = Vec::with_capacity(tf_name.len() + 1 + ncells * 8);
    write_csv_field(&mut buf, tf_name, delimiter);
    match site_ratio {
        None => {
            for _ in 0..ncells {
                buf.push(delimiter);
                buf.extend_from_slice(b"0.000");
            }
        }
        Some(site_ratio) => {
            // Wrap the Vec<u8> in a small adapter so we can use fmt::Write
            // and the std float formatter (which is the same algorithm used
            // by format!("{:.3}", x) but writes into our buffer with no
            // per-element heap allocation).
            struct ByteWriter<'a>(&'a mut Vec<u8>);
            impl<'a> std::fmt::Write for ByteWriter<'a> {
                fn write_str(&mut self, s: &str) -> std::fmt::Result {
                    self.0.extend_from_slice(s.as_bytes());
                    Ok(())
                }
            }

            let site_col = count_matrix.column(i);
            let control_col = count_matrix.column(nsites + i);
            for (&site_count, &control_count) in site_col.iter().zip(control_col.iter()) {
                let s = site_count as f32;
                let c = control_count as f32;
                let s_eff = s + pseudocount;
                let c_eff = c + pseudocount * site_ratio;
                let raw_lfc = (s_eff / (c / site_ratio + pseudocount)).log2();
                let prec = s_eff * c_eff / (s_eff + c_eff);
                let weight = prec / (prec + lambda);

                buf.push(delimiter);
                let _ = write!(ByteWriter(&mut buf), "{:.3}", weight * raw_lfc);
            }
        }
    }
    buf.push(b'\n');
    buf
}

fn validate_thread_count(num_threads: usize) -> io::Result<()> {
    if num_threads == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "threads must be >= 1",
        ));
    }
    Ok(())
}

fn validate_scoring_parameters(lambda: f32, pseudocount: f32) -> io::Result<()> {
    if !lambda.is_finite() || lambda < 0.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "lambda must be finite and >= 0",
        ));
    }
    if !pseudocount.is_finite() || pseudocount <= 0.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "pseudocount must be finite and > 0",
        ));
    }
    Ok(())
}

type IndexData = (
    FxHashMap<String, Lapper<u32, u32>>, // chromosome_trees
    Vec<u32>,                            // site_indices (1:1 with site intervals)
    FlatIndVec,                          // control_indices (variable per control interval)
    usize,                               // num_site_intervals (to distinguish sites from controls)
    Vec<String>,                         // tf_list
    Vec<u32>,                            // n_instances
    FxHashMap<String, u32>,              // chromosome_sizes
);

// Open an index file and return a reader over its decompressed bytes.
// Auto-detects the wrapper format from the first few bytes:
//   * 0x1F 0x8B + FEXTRA + "BC" subfield → BGZF (parallel-decompressed)
//   * 0x1F 0x8B (plain)                  → gzip (MultiGzDecoder)
//   * 0xFD                               → xz
//   * b'B' (start of "BACH" magic)       → uncompressed binary index file
//
// xz gives the smallest files but decompresses at ~75 MB/s on a single
// stream; gzip/BGZF trade ~30% larger files for ~8x faster decode, and
// BGZF in particular parallelises across cores via noodles. Uncompressed
// is the fastest to load (~0 s) but uses ~1 GB of disk for a typical
// human-genome index.
fn open_index_decompressed(index_file: &Path, bgzf_workers: usize) -> io::Result<Box<dyn Read>> {
    let mut probe = File::open(index_file)?;
    let mut header = [0u8; 16];
    let n = probe.read(&mut header)?;
    drop(probe);

    if n == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Empty index file",
        ));
    }

    let file = File::open(index_file)?;
    match header[0] {
        0x1F if n >= 2 && header[1] == 0x8B => {
            // gzip magic. Check for BGZF: FEXTRA flag (bit 2 of flags byte)
            // and the "BC" extra-subfield ID at offset 12.
            let is_bgzf =
                n >= 16 && (header[3] & 0x04) != 0 && header[12] == b'B' && header[13] == b'C';
            if is_bgzf {
                let workers = std::num::NonZeroUsize::new(bgzf_workers.max(1))
                    .unwrap_or(std::num::NonZeroUsize::MIN);
                let reader =
                    noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, file);
                Ok(Box::new(reader))
            } else {
                Ok(Box::new(MultiGzDecoder::new(file)))
            }
        }
        0xFD => Ok(Box::new(XzDecoder::new(file))),
        b'B' => Ok(Box::new(file)),
        b => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Unknown index format (first byte 0x{:02x})", b),
        )),
    }
}

// In-place inverse of per-chromosome delta encoding for u32 LE start arrays.
// Recovers absolute starts via a cumulative sum that resets at every change
// of chromosome id.
fn cumsum_per_chrom_in_place(buf: &mut [u8], chrom_ids: &[u8]) {
    let n = chrom_ids.len();
    if n == 0 {
        return;
    }
    let mut acc: u32 = 0;
    let mut prev_chrom: i32 = -1;
    for i in 0..n {
        let chrom = chrom_ids[i] as i32;
        let off = i * 4;
        let val = u32::from_le_bytes([buf[off], buf[off + 1], buf[off + 2], buf[off + 3]]);
        if chrom != prev_chrom {
            acc = val;
        } else {
            acc = acc.wrapping_add(val);
        }
        let bytes = acc.to_le_bytes();
        buf[off] = bytes[0];
        buf[off + 1] = bytes[1];
        buf[off + 2] = bytes[2];
        buf[off + 3] = bytes[3];
        prev_chrom = chrom;
    }
}

// In-place inverse of within-row delta encoding for the flattened u16 LE
// control TF index stream. Recovers absolute indices via a cumulative sum
// that resets at every row boundary defined by `n_per_row`.
fn cumsum_within_row_u16_in_place(buf: &mut [u8], n_per_row: &[usize]) {
    let mut tf_cursor: usize = 0;
    for &n_idx in n_per_row {
        let mut acc: u16 = 0;
        for j in 0..n_idx {
            let off = (tf_cursor + j) * 2;
            let val = u16::from_le_bytes([buf[off], buf[off + 1]]);
            if j == 0 {
                acc = val;
            } else {
                acc = acc.wrapping_add(val);
            }
            let bytes = acc.to_le_bytes();
            buf[off] = bytes[0];
            buf[off + 1] = bytes[1];
        }
        tf_cursor += n_idx;
    }
}

// Verify every chromosome ID in a columnar block addresses a row of the
// index's chromosome table. IDs are used as direct vector indices, so an
// out-of-range value in a corrupt file must be an error, not a panic.
fn check_chrom_ids(chrom_ids: &[u8], n_chroms: usize, what: &str) -> io::Result<()> {
    if let Some(&bad) = chrom_ids.iter().find(|&&c| c as usize >= n_chroms) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Malformed index: {what} chromosome ID {bad} is outside the chromosome \
                 table (which has {n_chroms} entries)"
            ),
        ));
    }
    Ok(())
}

fn one_based_start_to_zero_based(start: u32) -> io::Result<u32> {
    start.checked_sub(1).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "Index coordinates are 1-based; start position must be >= 1",
        )
    })
}

// Read and verify the "BACH" magic
fn read_magic(reader: &mut impl Read) -> io::Result<()> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != b"BACH" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Not a Bacha binary index (missing 'BACH' magic bytes)",
        ));
    }
    Ok(())
}

pub fn parse_file_to_lapper(index_file: &Path, num_threads: usize) -> io::Result<IndexData> {
    let mut reader = BufReader::new(open_index_decompressed(
        index_file,
        bgzf_worker_count(num_threads),
    )?);
    read_magic(&mut reader)?;
    parse_binary_index(&mut reader)
}

pub fn view_index(index_file: &Path, header_only: bool) -> io::Result<()> {
    // `view`/`summary` take no thread count; decode with the plain default.
    let mut reader = BufReader::new(open_index_decompressed(
        index_file,
        bgzf_worker_count(DEFAULT_BGZF_WORKERS),
    )?);
    read_magic(&mut reader)?;
    view_binary_index(&mut reader, header_only)
}

fn view_binary_index(reader: &mut impl Read, header_only: bool) -> io::Result<()> {
    fn read_u8(r: &mut impl Read) -> io::Result<u8> {
        let mut buf = [0u8; 1];
        r.read_exact(&mut buf)?;
        Ok(buf[0])
    }
    fn read_u16_le(r: &mut impl Read) -> io::Result<u16> {
        let mut buf = [0u8; 2];
        r.read_exact(&mut buf)?;
        Ok(u16::from_le_bytes(buf))
    }
    fn read_u32_le(r: &mut impl Read) -> io::Result<u32> {
        let mut buf = [0u8; 4];
        r.read_exact(&mut buf)?;
        Ok(u32::from_le_bytes(buf))
    }
    fn read_string(r: &mut impl Read) -> io::Result<String> {
        let len = read_u16_le(r)? as usize;
        let mut buf = vec![0u8; len];
        r.read_exact(&mut buf)?;
        String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    let version = read_u8(reader)?;
    if version != 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Unsupported binary index version: {version} (expected 1)"),
        ));
    }
    let genome_name = read_string(reader)?;
    let motif_set = read_string(reader)?;
    let window = read_u32_le(reader)?;
    let n_tfs = read_u16_le(reader)? as usize;
    let n_chroms = read_u8(reader)? as usize;

    writeln!(out, "## format: binary v{}", version)?;
    writeln!(out, "## genome: {}", genome_name)?;
    writeln!(out, "## motif set: {}", motif_set)?;
    writeln!(out, "## window: {}", window)?;

    // TF names
    let mut tf_names = Vec::with_capacity(n_tfs);
    for _ in 0..n_tfs {
        tf_names.push(read_string(reader)?);
    }

    // Instance counts
    let mut site_counts = Vec::with_capacity(n_tfs);
    for _ in 0..n_tfs {
        site_counts.push(read_u32_le(reader)?);
    }
    let mut control_counts = Vec::with_capacity(n_tfs);
    for _ in 0..n_tfs {
        control_counts.push(read_u32_le(reader)?);
    }

    // Chromosome table
    let mut chrom_names = Vec::with_capacity(n_chroms);
    let mut chrom_sizes = Vec::with_capacity(n_chroms);
    for _ in 0..n_chroms {
        chrom_names.push(read_string(reader)?);
        chrom_sizes.push(read_u32_le(reader)?);
    }

    // Print chromosome sizes
    writeln!(out, "##[chromosome sizes]")?;
    for (name, size) in chrom_names.iter().zip(chrom_sizes.iter()) {
        writeln!(out, "{}\t{}", name, size)?;
    }

    // Print sites header
    writeln!(out, "##[sites]")?;
    writeln!(out, "# id\tindex\tcount")?;
    for (i, name) in tf_names.iter().enumerate() {
        writeln!(out, "{}\t{}\t{}", name, i + 1, site_counts[i])?;
    }

    // Print controls header
    writeln!(out, "##[controls]")?;
    writeln!(out, "# id\tindex\tcount")?;
    for (i, name) in tf_names.iter().enumerate() {
        writeln!(out, "{}\t{}\t{}", name, n_tfs + i + 1, control_counts[i])?;
    }

    if header_only {
        return Ok(());
    }

    // Flush before expensive bulk reads so SIGPIPE fires early when piped to head/less
    out.flush()?;

    // Check if the pipe reader (e.g. head) has already exited before doing expensive bulk reads.
    // Without this, the program blocks in read_exact while no stdout writes happen to trigger SIGPIPE.
    #[cfg(unix)]
    {
        let mut pfd = libc::pollfd {
            fd: 1,
            events: 0,
            revents: 0,
        };
        unsafe { libc::poll(&mut pfd, 1, 0) };
        if pfd.revents & (libc::POLLERR | libc::POLLHUP) != 0 {
            return Ok(());
        }
    }

    // Site positions
    let n_sites = read_u32_le(reader)? as usize;
    let mut site_chrom_ids = vec![0u8; n_sites];
    reader.read_exact(&mut site_chrom_ids)?;
    check_chrom_ids(&site_chrom_ids, chrom_names.len(), "site")?;
    let mut site_starts_buf = vec![0u8; n_sites * 4];
    reader.read_exact(&mut site_starts_buf)?;
    let mut site_tf_buf = vec![0u8; n_sites * 2];
    reader.read_exact(&mut site_tf_buf)?;
    cumsum_per_chrom_in_place(&mut site_starts_buf, &site_chrom_ids);

    writeln!(out, "##[site positions]")?;
    writeln!(out, "# chrom\tstart\tend\tindex")?;
    for i in 0..n_sites {
        let chrom = &chrom_names[site_chrom_ids[i] as usize];
        let start = u32::from_le_bytes([
            site_starts_buf[i * 4],
            site_starts_buf[i * 4 + 1],
            site_starts_buf[i * 4 + 2],
            site_starts_buf[i * 4 + 3],
        ]);
        let tf_idx = u16::from_le_bytes([site_tf_buf[i * 2], site_tf_buf[i * 2 + 1]]);
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            chrom,
            start,
            start + window - 1,
            tf_idx
        )?;
    }

    // Flush before next bulk read
    out.flush()?;

    #[cfg(unix)]
    {
        let mut pfd = libc::pollfd {
            fd: 1,
            events: 0,
            revents: 0,
        };
        unsafe { libc::poll(&mut pfd, 1, 0) };
        if pfd.revents & (libc::POLLERR | libc::POLLHUP) != 0 {
            return Ok(());
        }
    }

    // Control positions
    let n_controls = read_u32_le(reader)? as usize;
    let mut ctrl_chrom_ids = vec![0u8; n_controls];
    reader.read_exact(&mut ctrl_chrom_ids)?;
    check_chrom_ids(&ctrl_chrom_ids, chrom_names.len(), "control")?;
    let mut ctrl_starts_buf = vec![0u8; n_controls * 4];
    reader.read_exact(&mut ctrl_starts_buf)?;
    cumsum_per_chrom_in_place(&mut ctrl_starts_buf, &ctrl_chrom_ids);
    let n_per_row: Vec<usize> = {
        let mut buf = vec![0u8; n_controls * 2];
        reader.read_exact(&mut buf)?;
        (0..n_controls)
            .map(|i| u16::from_le_bytes([buf[i * 2], buf[i * 2 + 1]]) as usize)
            .collect()
    };
    let total_ctrl_indices = read_u32_le(reader)? as usize;
    let n_per_row_total: usize = n_per_row.iter().sum();
    if n_per_row_total != total_ctrl_indices {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Malformed index: control rows declare {n_per_row_total} feature index/indices \
                 in total but the index states {total_ctrl_indices}"
            ),
        ));
    }
    let mut ctrl_tf_buf = vec![0u8; total_ctrl_indices * 2];
    reader.read_exact(&mut ctrl_tf_buf)?;
    cumsum_within_row_u16_in_place(&mut ctrl_tf_buf, &n_per_row);

    writeln!(out, "##[control positions]")?;
    writeln!(out, "# chrom\tstart\tend\tindex")?;
    let mut tf_cursor: usize = 0;
    for i in 0..n_controls {
        let chrom = &chrom_names[ctrl_chrom_ids[i] as usize];
        let start = u32::from_le_bytes([
            ctrl_starts_buf[i * 4],
            ctrl_starts_buf[i * 4 + 1],
            ctrl_starts_buf[i * 4 + 2],
            ctrl_starts_buf[i * 4 + 3],
        ]);
        let n_idx = n_per_row[i];

        // Build comma-separated index list
        let mut indices = String::new();
        for j in 0..n_idx {
            let tf_idx = u16::from_le_bytes([
                ctrl_tf_buf[(tf_cursor + j) * 2],
                ctrl_tf_buf[(tf_cursor + j) * 2 + 1],
            ]);
            if j > 0 {
                indices.push(',');
            }
            indices.push_str(&tf_idx.to_string());
        }
        tf_cursor += n_idx;

        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            chrom,
            start,
            start + window - 1,
            indices
        )?;
    }

    Ok(())
}

#[allow(unused_assignments)]
fn parse_binary_index(reader: &mut impl Read) -> io::Result<IndexData> {
    let perf = std::env::var("BACHA_PROFILE").is_ok();
    let t_total = Instant::now();
    let mut t_phase = Instant::now();
    macro_rules! tick {
        ($name:expr) => {
            if perf {
                eprintln!(
                    "[PERF]   {:<28}={:.3}s",
                    $name,
                    t_phase.elapsed().as_secs_f64()
                );
            }
            t_phase = Instant::now();
        };
    }
    // Helper functions
    fn read_u8(r: &mut impl Read) -> io::Result<u8> {
        let mut buf = [0u8; 1];
        r.read_exact(&mut buf)?;
        Ok(buf[0])
    }
    fn read_u16_le(r: &mut impl Read) -> io::Result<u16> {
        let mut buf = [0u8; 2];
        r.read_exact(&mut buf)?;
        Ok(u16::from_le_bytes(buf))
    }
    fn read_u32_le(r: &mut impl Read) -> io::Result<u32> {
        let mut buf = [0u8; 4];
        r.read_exact(&mut buf)?;
        Ok(u32::from_le_bytes(buf))
    }
    fn read_string(r: &mut impl Read) -> io::Result<String> {
        let len = read_u16_le(r)? as usize;
        let mut buf = vec![0u8; len];
        r.read_exact(&mut buf)?;
        String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    // Version. Only v1 (the released format) is supported; the pre-release
    // formats never shipped, so there is nothing else to read.
    let version = read_u8(reader)?;
    if version != 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Unsupported binary index version: {version} (expected 1)"),
        ));
    }

    // Header
    let genome_name = read_string(reader)?;
    let motif_set = read_string(reader)?;
    let window = read_u32_le(reader)?;
    let n_tfs = read_u16_le(reader)? as usize;
    let n_chroms = read_u8(reader)? as usize;

    eprintln!("Index genome: {}", genome_name);
    eprintln!("Index motif set: {}", motif_set);

    // TF names
    let mut tf_list = Vec::with_capacity(n_tfs);
    for _ in 0..n_tfs {
        tf_list.push(read_string(reader)?);
    }

    // Instance counts: site_counts then control_counts
    let mut n_instances = Vec::with_capacity(n_tfs * 2);
    for _ in 0..(n_tfs * 2) {
        n_instances.push(read_u32_le(reader)?);
    }

    // Chromosome table
    let mut chrom_names = Vec::with_capacity(n_chroms);
    let mut chrom_sizes: FxHashMap<String, u32> = FxHashMap::default();
    for _ in 0..n_chroms {
        let name = read_string(reader)?;
        let size = read_u32_le(reader)?;
        chrom_sizes.insert(name.clone(), size);
        chrom_names.push(name);
    }
    tick!("header");

    // Site positions (columnar)
    let n_sites = read_u32_le(reader)? as usize;
    let mut site_chrom_ids = vec![0u8; n_sites];
    reader.read_exact(&mut site_chrom_ids)?;
    let mut site_starts_buf = vec![0u8; n_sites * 4];
    reader.read_exact(&mut site_starts_buf)?;
    let mut site_tf_buf = vec![0u8; n_sites * 2];
    reader.read_exact(&mut site_tf_buf)?;
    // Chromosome IDs index directly into fixed-size per-chromosome vectors
    // below, so an out-of-range ID in a corrupt file would panic rather than
    // report an error. Validate against the header's chromosome count first.
    check_chrom_ids(&site_chrom_ids, chrom_names.len(), "site")?;
    tick!("read site columnar");
    cumsum_per_chrom_in_place(&mut site_starts_buf, &site_chrom_ids);
    tick!("cumsum sites");

    // Accumulate intervals into Vec<Vec<Interval>> indexed by chrom_id (u8).
    // This avoids per-fragment String::clone + HashMap::entry probes on the
    // hot loop (~41M+ iterations); we materialise the FxHashMap<String, _>
    // shape that build_index_result expects in one bulk pass at the end.
    let n_chroms = chrom_names.len();
    let mut chromosome_intervals: Vec<Vec<Interval<u32, u32>>> =
        (0..n_chroms).map(|_| Vec::new()).collect();
    // Pre-count per-chrom site occurrences so each per-chrom Vec is sized
    // once instead of doubling repeatedly during the push loop below.
    let mut per_chrom_site_counts = vec![0usize; n_chroms];
    for &cid in &site_chrom_ids {
        per_chrom_site_counts[cid as usize] += 1;
    }
    for (cid, &cnt) in per_chrom_site_counts.iter().enumerate() {
        chromosome_intervals[cid].reserve(cnt);
    }

    let mut site_indices: Vec<u32> = Vec::with_capacity(n_sites);
    for i in 0..n_sites {
        let chrom_id = site_chrom_ids[i] as usize;
        let start = u32::from_le_bytes([
            site_starts_buf[i * 4],
            site_starts_buf[i * 4 + 1],
            site_starts_buf[i * 4 + 2],
            site_starts_buf[i * 4 + 3],
        ]);
        let start0 = one_based_start_to_zero_based(start)?;
        let tf_idx = u16::from_le_bytes([site_tf_buf[i * 2], site_tf_buf[i * 2 + 1]]) as u32;

        chromosome_intervals[chrom_id].push(Interval {
            start: start0,
            stop: start0 + window,
            val: i as u32,
        });
        site_indices.push(tf_idx);
    }
    let num_site_intervals = n_sites;
    tick!("build site intervals");

    // Control positions (columnar)
    let n_controls = read_u32_le(reader)? as usize;
    let mut ctrl_chrom_ids = vec![0u8; n_controls];
    reader.read_exact(&mut ctrl_chrom_ids)?;
    let mut ctrl_starts_buf = vec![0u8; n_controls * 4];
    reader.read_exact(&mut ctrl_starts_buf)?;
    check_chrom_ids(&ctrl_chrom_ids, chrom_names.len(), "control")?;
    tick!("read ctrl chrom+starts");
    cumsum_per_chrom_in_place(&mut ctrl_starts_buf, &ctrl_chrom_ids);
    tick!("cumsum ctrl starts");
    let n_per_row: Vec<usize> = {
        let mut buf = vec![0u8; n_controls * 2];
        reader.read_exact(&mut buf)?;
        (0..n_controls)
            .map(|i| u16::from_le_bytes([buf[i * 2], buf[i * 2 + 1]]) as usize)
            .collect()
    };
    tick!("read n_per_row");
    let total_ctrl_indices = read_u32_le(reader)? as usize;
    // `n_per_row` and `total_indices` are independent fields in the file, and
    // the within-row cumsum below walks `n_per_row` across `ctrl_tf_buf`. If
    // they disagree the walk runs off the end of the buffer and panics, so
    // reconcile them before touching it.
    let n_per_row_total: usize = n_per_row.iter().sum();
    if n_per_row_total != total_ctrl_indices {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Malformed index: control rows declare {n_per_row_total} feature index/indices \
                 in total but the index states {total_ctrl_indices}"
            ),
        ));
    }
    let mut ctrl_tf_buf = vec![0u8; total_ctrl_indices * 2];
    reader.read_exact(&mut ctrl_tf_buf)?;
    tick!("read ctrl tf buf");
    cumsum_within_row_u16_in_place(&mut ctrl_tf_buf, &n_per_row);
    tick!("cumsum within row");

    // Pre-size each chrom's vec for the additional control intervals it
    // will receive. Combined with the earlier site reservation this means
    // the chromosome_intervals push loop should hit zero reallocations.
    let mut per_chrom_ctrl_counts = vec![0usize; n_chroms];
    for &cid in &ctrl_chrom_ids {
        per_chrom_ctrl_counts[cid as usize] += 1;
    }
    for (cid, &cnt) in per_chrom_ctrl_counts.iter().enumerate() {
        chromosome_intervals[cid].reserve(cnt);
    }

    // control_data is just the cumsum-corrected ctrl_tf_buf widened to u32.
    // Computing it in one tight loop is much cheaper than the previous
    // nested-per-row push.
    let mut control_data: Vec<u32> = Vec::with_capacity(total_ctrl_indices);
    for i in 0..total_ctrl_indices {
        let v = u16::from_le_bytes([ctrl_tf_buf[i * 2], ctrl_tf_buf[i * 2 + 1]]) as u32;
        control_data.push(v);
    }
    // control_offsets is the prefix sum of n_per_row.
    let mut control_offsets: Vec<u32> = Vec::with_capacity(n_controls + 1);
    control_offsets.push(0);
    let mut acc: u32 = 0;
    for &n in &n_per_row {
        acc += n as u32;
        control_offsets.push(acc);
    }

    tick!("convert ctrl_data u16->u32");

    let site_intervals_offset = n_sites as u32;
    for i in 0..n_controls {
        let chrom_id = ctrl_chrom_ids[i] as usize;
        let start = u32::from_le_bytes([
            ctrl_starts_buf[i * 4],
            ctrl_starts_buf[i * 4 + 1],
            ctrl_starts_buf[i * 4 + 2],
            ctrl_starts_buf[i * 4 + 3],
        ]);
        let start0 = one_based_start_to_zero_based(start)?;

        chromosome_intervals[chrom_id].push(Interval {
            start: start0,
            stop: start0 + window,
            val: site_intervals_offset + i as u32,
        });
    }

    tick!("build ctrl data+ivs");

    // Materialise the FxHashMap<String, Vec<Interval>> shape that
    // build_index_result expects. Empty chroms are dropped.
    let chromosome_trees: FxHashMap<String, Vec<Interval<u32, u32>>> = chromosome_intervals
        .into_iter()
        .enumerate()
        .filter(|(_, ivs)| !ivs.is_empty())
        .map(|(i, ivs)| (chrom_names[i].clone(), ivs))
        .collect();
    tick!("convert to FxHashMap");

    info!("Finished reading binary index file");
    let result = build_index_result(
        chromosome_trees,
        site_indices,
        control_data,
        control_offsets,
        num_site_intervals,
        tf_list,
        n_instances,
        chrom_sizes,
    );
    tick!("build_index_result");
    if perf {
        eprintln!(
            "[PERF]   {:<28}={:.3}s",
            "TOTAL parse_binary_index",
            t_total.elapsed().as_secs_f64()
        );
    }
    result
}

fn build_index_result(
    chromosome_trees: FxHashMap<String, Vec<Interval<u32, u32>>>,
    site_indices: Vec<u32>,
    control_data: Vec<u32>,
    control_offsets: Vec<u32>,
    num_site_intervals: usize,
    tf_order: Vec<String>,
    n_instances: Vec<u32>,
    chrom_sizes: FxHashMap<String, u32>,
) -> io::Result<IndexData> {
    // Feature-index validation. The counting hot loop turns each index
    // straight into a column of the count matrix, so establishing the range
    // here is what lets `process_overlaps` skip any per-hit bounds test.
    // Structural consistency (array lengths, offset arrays) is guaranteed by
    // construction in `parse_binary_index`, which validates the raw header
    // fields against each other as it reads them.
    let n_tfs = tf_order.len();
    if n_tfs == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Malformed index: no features defined",
        ));
    }
    // Sites carry a 1-based feature index; controls carry the same index
    // offset by n_tfs so both share one ID space (see the index format docs).
    if let Some(&bad) = site_indices.iter().find(|&&i| i == 0 || i as usize > n_tfs) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Malformed index: site feature index {bad} is outside the valid range 1..={n_tfs}"
            ),
        ));
    }
    if let Some(&bad) = control_data
        .iter()
        .find(|&&i| i as usize <= n_tfs || i as usize > 2 * n_tfs)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Malformed index: control feature index {bad} is outside the valid range {}..={}",
                n_tfs + 1,
                n_tfs * 2
            ),
        ));
    }

    // Sort each chromosome's intervals and build its Lapper in parallel;
    // chroms are independent and the sort dominates total index-load time.
    let lapper_map: FxHashMap<String, Lapper<u32, u32>> = chromosome_trees
        .into_par_iter()
        .map(|(chr, mut intervals)| {
            info!(
                "Chromosome: {}, Number of intervals: {}",
                chr,
                intervals.len()
            );
            intervals.sort_by_key(|iv| iv.start);
            (chr, Lapper::new(intervals))
        })
        .collect();

    let control_indices = FlatIndVec {
        data: control_data,
        offsets: control_offsets,
    };

    if !chrom_sizes.is_empty() {
        eprintln!(
            "Loaded chromosome sizes for {} chromosomes",
            chrom_sizes.len()
        );
    }

    Ok((
        lapper_map,
        site_indices,
        control_indices,
        num_site_intervals,
        tf_order,
        n_instances,
        chrom_sizes,
    ))
}

#[test]
fn test_scoring_parameter_validation() {
    assert!(validate_scoring_parameters(0.0, 1.0).is_ok());
    assert!(validate_scoring_parameters(10.0, 1.0).is_ok());
    assert!(validate_scoring_parameters(-1.0, 1.0).is_err());
    assert!(validate_scoring_parameters(f32::NAN, 1.0).is_err());
    assert!(validate_scoring_parameters(f32::INFINITY, 1.0).is_err());
    assert!(validate_scoring_parameters(10.0, 0.0).is_err());
    assert!(validate_scoring_parameters(10.0, -1.0).is_err());
    assert!(validate_scoring_parameters(10.0, f32::NAN).is_err());
    assert!(validate_scoring_parameters(10.0, f32::INFINITY).is_err());
}

#[test]
fn test_assembly_mismatch_error_message() {
    let err = assembly_mismatch_error("chr19", 62_000_000, 61_431_566);
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    let msg = err.to_string();
    assert!(msg.contains("chr19"));
    assert!(msg.contains("62000000"));
    assert!(msg.contains("61431566"));
}

// ---------------------------------------------------------------------------
// Test fixtures
// ---------------------------------------------------------------------------

// A binary index builder mirroring what `bachaR::write_index()` produces (see
// the format documentation on `BachaIndex`). It writes exactly what it is
// given without sanitising, so malformed-index handling can be exercised, and
// gzips the result so the reader's decompression path is exercised too.
//
// Coordinates are 1-based starts, matching the on-disk format; every interval
// is `window` bp wide. Feature indices are 1-based, with control indices
// offset by `n_tfs` so sites and controls share one ID space.
#[cfg(test)]
#[derive(Clone)]
struct TestIndex {
    window: u32,
    tf_names: Vec<String>,
    chroms: Vec<(String, u32)>,
    sites: Vec<(u8, u32, u16)>,
    controls: Vec<(u8, u32, Vec<u16>)>,
    // Overrides the `total_indices` header field, to test a file whose
    // per-row counts and total disagree.
    total_indices_override: Option<u32>,
}

#[cfg(test)]
impl TestIndex {
    // One TF on `n_chroms` chromosomes of `chrom_len` bp, each carrying a
    // single 1 bp site and a single 1 bp control at the given 1-based
    // positions. Chromosomes are named chr1..chrN.
    fn simple(n_chroms: u8, chrom_len: u32, site_at: u32, control_at: u32) -> Self {
        TestIndex {
            window: 1,
            tf_names: vec!["TF1".to_string()],
            chroms: (0..n_chroms)
                .map(|i| (format!("chr{}", i + 1), chrom_len))
                .collect(),
            sites: (0..n_chroms).map(|c| (c, site_at, 1)).collect(),
            controls: (0..n_chroms).map(|c| (c, control_at, vec![2])).collect(),
            total_indices_override: None,
        }
    }

    fn write(&self, path: &Path) -> io::Result<()> {
        let n_tfs = self.tf_names.len();
        let mut b: Vec<u8> = Vec::new();

        let put_str = |b: &mut Vec<u8>, s: &str| {
            b.extend_from_slice(&(s.len() as u16).to_le_bytes());
            b.extend_from_slice(s.as_bytes());
        };

        b.extend_from_slice(b"BACH");
        b.push(1u8); // format version
        put_str(&mut b, "TestGenome");
        put_str(&mut b, "TestMotifs");
        b.extend_from_slice(&self.window.to_le_bytes());
        b.extend_from_slice(&(n_tfs as u16).to_le_bytes());
        b.push(self.chroms.len() as u8);
        for name in &self.tf_names {
            put_str(&mut b, name);
        }

        // Instance counts: per-TF site totals, then per-TF control totals.
        let mut site_counts = vec![0u32; n_tfs];
        for &(_, _, tf) in &self.sites {
            if (1..=n_tfs as u16).contains(&tf) {
                site_counts[tf as usize - 1] += 1;
            }
        }
        let mut control_counts = vec![0u32; n_tfs];
        for (_, _, idx) in &self.controls {
            for &tf in idx {
                let one_based = tf as usize;
                if one_based > n_tfs && one_based <= 2 * n_tfs {
                    control_counts[one_based - n_tfs - 1] += 1;
                }
            }
        }
        for c in site_counts.iter().chain(control_counts.iter()) {
            b.extend_from_slice(&c.to_le_bytes());
        }
        for (name, len) in &self.chroms {
            put_str(&mut b, name);
            b.extend_from_slice(&len.to_le_bytes());
        }

        // Starts are delta-encoded within each chromosome run, so entries must
        // be grouped by chromosome; sort to guarantee that.
        let delta_per_chrom = |entries: &[(u8, u32)]| -> Vec<u32> {
            let mut out = Vec::with_capacity(entries.len());
            let mut prev: Option<(u8, u32)> = None;
            for &(c, s) in entries {
                out.push(match prev {
                    Some((pc, ps)) if pc == c => s.wrapping_sub(ps),
                    _ => s,
                });
                prev = Some((c, s));
            }
            out
        };

        let mut sites = self.sites.clone();
        sites.sort_by_key(|&(c, s, _)| (c, s));
        b.extend_from_slice(&(sites.len() as u32).to_le_bytes());
        for &(c, _, _) in &sites {
            b.push(c);
        }
        let site_keys: Vec<(u8, u32)> = sites.iter().map(|&(c, s, _)| (c, s)).collect();
        for d in delta_per_chrom(&site_keys) {
            b.extend_from_slice(&d.to_le_bytes());
        }
        for &(_, _, tf) in &sites {
            b.extend_from_slice(&tf.to_le_bytes());
        }

        let mut controls = self.controls.clone();
        controls.sort_by_key(|(c, s, _)| (*c, *s));
        b.extend_from_slice(&(controls.len() as u32).to_le_bytes());
        for (c, _, _) in &controls {
            b.push(*c);
        }
        let ctrl_keys: Vec<(u8, u32)> = controls.iter().map(|(c, s, _)| (*c, *s)).collect();
        for d in delta_per_chrom(&ctrl_keys) {
            b.extend_from_slice(&d.to_le_bytes());
        }
        for (_, _, idx) in &controls {
            b.extend_from_slice(&(idx.len() as u16).to_le_bytes());
        }
        let total: u32 = controls.iter().map(|(_, _, i)| i.len() as u32).sum();
        b.extend_from_slice(&self.total_indices_override.unwrap_or(total).to_le_bytes());
        // Feature indices are delta-encoded within each control row.
        for (_, _, idx) in &controls {
            let mut prev = 0u16;
            for (j, &tf) in idx.iter().enumerate() {
                b.extend_from_slice(
                    &(if j == 0 { tf } else { tf.wrapping_sub(prev) }).to_le_bytes(),
                );
                prev = tf;
            }
        }

        let file = File::create(path)?;
        let mut encoder = flate2::write::GzEncoder::new(file, Compression::default());
        encoder.write_all(&b)?;
        encoder.finish()?;
        Ok(())
    }
}

#[cfg(test)]
fn write_gzip_text(path: &Path, text: &str) -> io::Result<()> {
    let file = File::create(path)?;
    let mut encoder = flate2::write::GzEncoder::new(file, Compression::default());
    encoder.write_all(text.as_bytes())?;
    encoder.finish()?;
    Ok(())
}

// Write `text` as BGZF: a sequence of complete gzip members terminated by the
// BGZF end-of-file marker. Uses flate2 with the FEXTRA "BC" subfield so the
// result is detected as BGZF by `is_bgzf_file`, which is all these tests need
// (they do not depend on the block-size bookkeeping htslib performs).
#[cfg(test)]
fn write_bgzf_text(path: &Path, text: &str, with_eof_marker: bool) -> io::Result<()> {
    use flate2::write::DeflateEncoder;

    let mut deflated = Vec::new();
    let mut enc = DeflateEncoder::new(&mut deflated, Compression::default());
    enc.write_all(text.as_bytes())?;
    enc.finish()?;

    let mut crc = flate2::Crc::new();
    crc.update(text.as_bytes());

    // BGZF block: gzip header with FEXTRA carrying BSIZE-1, deflate payload,
    // CRC32 and ISIZE trailer.
    let bsize = (12 + 6 + deflated.len() + 8) as u16;
    let mut out: Vec<u8> = vec![0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff];
    out.extend_from_slice(&6u16.to_le_bytes()); // XLEN
    out.extend_from_slice(b"BC");
    out.extend_from_slice(&2u16.to_le_bytes()); // SLEN
    out.extend_from_slice(&(bsize - 1).to_le_bytes());
    out.extend_from_slice(&deflated);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&(text.len() as u32).to_le_bytes());
    if with_eof_marker {
        out.extend_from_slice(&BGZF_EOF);
    }

    std::fs::write(path, &out)
}

#[cfg(test)]
fn temp_dir(name: &str) -> io::Result<std::path::PathBuf> {
    let dir = std::env::temp_dir().join(format!("bacha_{name}_test_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir)?;
    Ok(dir)
}

// Write an index plus fragment and cell files into a fresh directory and
// return their paths.
#[cfg(test)]
fn fixture(
    name: &str,
    index: &TestIndex,
    frag_lines: &str,
    cell_lines: &str,
) -> io::Result<(std::path::PathBuf, std::path::PathBuf, std::path::PathBuf)> {
    let dir = temp_dir(name)?;
    let index_path = dir.join("index.bi");
    let frag_path = dir.join("fragments.tsv.gz");
    let cell_path = dir.join("cells.txt");
    index.write(&index_path)?;
    write_gzip_text(&frag_path, frag_lines)?;
    std::fs::write(&cell_path, cell_lines)?;
    Ok((index_path, frag_path, cell_path))
}

// Multi-file fixture: writes one index, one fragment file and one cells file
// per `samples` entry, plus a manifest naming them all. `samples` is
// (sample name, fragment lines, cells file contents).
//
// Fragment and cells paths are written into the manifest relative to the
// manifest's own directory, which also exercises manifest-relative path
// resolution. Returns (index path, manifest path); the index comes first so
// the `remove_dir_all(index.parent())` cleanup idiom used by other tests works
// unchanged.
#[cfg(test)]
fn multi_fixture(
    name: &str,
    index: &TestIndex,
    samples: &[(&str, &str, &str)],
) -> io::Result<(std::path::PathBuf, std::path::PathBuf)> {
    let dir = temp_dir(name)?;
    let index_path = dir.join("index.bi");
    index.write(&index_path)?;

    let manifest_path = dir.join("manifest.tsv");
    let mut manifest = String::from("sample\tfragments\tcells\n");
    for (sample, frag_lines, cell_lines) in samples {
        let frag_rel = format!("{sample}_fragments.tsv.gz");
        let cell_rel = format!("{sample}_cells.txt");
        write_gzip_text(&dir.join(&frag_rel), frag_lines)?;
        std::fs::write(dir.join(&cell_rel), cell_lines)?;
        manifest.push_str(&format!("{sample}\t{frag_rel}\t{cell_rel}\n"));
    }
    std::fs::write(&manifest_path, manifest)?;
    Ok((index_path, manifest_path))
}

// Count a manifest and hand back the raw matrix, mirroring `fcount_core` for
// the multi-file path.
#[cfg(test)]
fn fcount_manifest(
    manifest_file: &Path,
    index_file: &Path,
    num_threads: usize,
    check_bounds: bool,
) -> io::Result<(Array2<u32>, Vec<String>, Vec<String>, Vec<u32>)> {
    let sources = build_sources(read_manifest(manifest_file)?)?;
    let (count_matrix, cell_order, tables) =
        count_sources(&sources, index_file, num_threads, check_bounds)?;
    let IndexTables {
        tf_list,
        n_instances,
        ..
    } = tables;
    Ok((count_matrix, cell_order, tf_list, n_instances))
}

// ---------------------------------------------------------------------------
// Multiple fragment files
// ---------------------------------------------------------------------------

// The point of the multi-file path: counting N files in one run must give
// exactly what counting them separately and placing the results side by side
// would, with the cells in manifest order.
#[test]
fn test_manifest_matches_separate_single_file_runs() -> io::Result<()> {
    let idx = TestIndex::simple(2, 1000, 100, 200);

    // Three "samples" with disjoint cells and different fragment counts, so a
    // mix-up of row blocks or cell offsets would show up as wrong counts
    // rather than merely reordered ones.
    let samples: [(&str, &str, &str); 3] = [
        (
            "S1",
            "chr1\t99\t200\tc1\nchr1\t99\t200\tc1\nchr2\t99\t200\tc2\n",
            "c1\nc2\n",
        ),
        ("S2", "chr1\t99\t200\td1\n", "d1\n"),
        (
            "S3",
            "chr2\t99\t200\te1\nchr2\t199\t300\te2\nchr1\t99\t200\te2\n",
            "e1\ne2\ne3\n",
        ),
    ];

    for threads in [1usize, 4] {
        let (index, manifest) =
            multi_fixture(&format!("manifest_vs_separate_{threads}"), &idx, &samples)?;
        let (combined, cell_order, tf_list, _) = fcount_manifest(&manifest, &index, threads, true)?;

        assert_eq!(tf_list.len(), 1);
        assert_eq!(
            cell_order,
            vec!["c1", "c2", "d1", "e1", "e2", "e3"],
            "-t {threads}: cells must be in manifest order, then file order"
        );
        assert_eq!(combined.nrows(), 6);

        // Now the same three files one at a time, stacked in the same order.
        let mut row = 0usize;
        for (sample, frag_lines, cell_lines) in &samples {
            let (sindex, sfrag, scells) = fixture(
                &format!("manifest_vs_separate_{threads}_{sample}"),
                &idx,
                frag_lines,
                cell_lines,
            )?;
            let (single, _, _, _) = fcount_core(&sfrag, &sindex, &scells, threads, true)?;
            for i in 0..single.nrows() {
                for j in 0..single.ncols() {
                    assert_eq!(
                        combined[[row + i, j]],
                        single[[i, j]],
                        "-t {threads}: sample {sample} row {i} col {j}"
                    );
                }
            }
            row += single.nrows();
            std::fs::remove_dir_all(sindex.parent().unwrap())?;
        }
        assert_eq!(row, combined.nrows());
        std::fs::remove_dir_all(index.parent().unwrap())?;
    }
    Ok(())
}

// A cell may belong to only one fragment file
#[test]
fn test_manifest_rejects_a_cell_in_two_files() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let (index, manifest) = multi_fixture(
        "manifest_shared_cell",
        &idx,
        &[
            ("S1", "chr1\t99\t200\tshared\n", "shared\n"),
            ("S2", "chr1\t99\t200\tshared\n", "shared\n"),
        ],
    )?;

    let err = fcount_manifest(&manifest, &index, 4, true)
        .expect_err("a cell listed for two samples must be rejected");
    assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    let msg = err.to_string();
    assert!(
        msg.contains("shared"),
        "message should name the cell: {msg}"
    );
    assert!(
        msg.contains("S1") && msg.contains("S2"),
        "should name both samples: {msg}"
    );

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

// The same raw barcode in two samples is routine (shared barcode whitelists).
// It is fine as long as the cells files give the two cells distinct output
// names, and each must then be scored only from its own file.
#[test]
fn test_manifest_allows_shared_barcode_with_distinct_output_names() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    // Both files use barcode "bc1"; S1 has two fragments over the site, S2 one.
    let (index, manifest) = multi_fixture(
        "manifest_shared_barcode",
        &idx,
        &[
            (
                "S1",
                "chr1\t99\t200\tbc1\nchr1\t99\t200\tbc1\n",
                "bc1\tS1_bc1\n",
            ),
            ("S2", "chr1\t99\t200\tbc1\n", "bc1\tS2_bc1\n"),
        ],
    )?;

    for threads in [1usize, 4] {
        let (m, cell_order, _, _) = fcount_manifest(&manifest, &index, threads, true)?;
        assert_eq!(cell_order, vec!["S1_bc1", "S2_bc1"]);
        // Site column: 2 counts for S1, 1 for S2 - no bleed between files.
        assert_eq!(m[[0, 0]], 2, "-t {threads}: S1 site");
        assert_eq!(m[[1, 0]], 1, "-t {threads}: S2 site");
        assert_eq!(m[[0, 1]], 2, "-t {threads}: S1 control");
        assert_eq!(m[[1, 1]], 1, "-t {threads}: S2 control");
    }

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

// A failure in any file must name the sample it came from
#[test]
fn test_manifest_error_names_the_offending_sample() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    // S2's fragments are all on a chromosome the index does not know, which is
    // the "no fragments matched any chromosome" hard error.
    let (index, manifest) = multi_fixture(
        "manifest_sample_named",
        &idx,
        &[
            ("good", "chr1\t99\t200\tc1\n", "c1\n"),
            ("bad", "contigX\t99\t200\td1\n", "d1\n"),
        ],
    )?;

    let err = fcount_manifest(&manifest, &index, 4, true)
        .expect_err("a sample whose fragments match no chromosome must fail");
    let msg = err.to_string();
    assert!(
        msg.contains("sample 'bad'"),
        "error should name the sample: {msg}"
    );
    assert!(
        msg.contains("different genome assemblies"),
        "original message should be preserved: {msg}"
    );

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

// Every fragment file is checked for truncation before the (expensive) index
// load, so a bad last file fails immediately instead of after counting all the
// earlier ones.
#[test]
fn test_manifest_rejects_truncated_file_before_counting() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let dir = temp_dir("manifest_truncated")?;
    let index_path = dir.join("index.bi");
    idx.write(&index_path)?;

    // First file is fine; the second is BGZF without its EOF marker.
    write_gzip_text(&dir.join("A_fragments.tsv.gz"), "chr1\t99\t200\tc1\n")?;
    std::fs::write(dir.join("A_cells.txt"), "c1\n")?;
    write_bgzf_text(
        &dir.join("B_fragments.tsv.gz"),
        "chr1\t99\t200\td1\n",
        false,
    )?;
    std::fs::write(dir.join("B_cells.txt"), "d1\n")?;
    std::fs::write(
        dir.join("manifest.tsv"),
        "A\tA_fragments.tsv.gz\tA_cells.txt\nB\tB_fragments.tsv.gz\tB_cells.txt\n",
    )?;

    let err = fcount_manifest(&dir.join("manifest.tsv"), &index_path, 4, true)
        .expect_err("a truncated fragment file must be rejected");
    let msg = err.to_string();
    assert!(msg.contains("sample 'B'"), "should name the sample: {msg}");
    assert!(msg.contains("truncated"), "should say truncated: {msg}");

    std::fs::remove_dir_all(&dir)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Manifest and cells file parsing
// ---------------------------------------------------------------------------

#[test]
fn test_manifest_parsing() -> io::Result<()> {
    let dir = temp_dir("manifest_parsing")?;
    std::fs::write(dir.join("f.tsv.gz"), "")?;
    std::fs::write(dir.join("c.txt"), "")?;

    let write = |name: &str, body: &str| -> io::Result<std::path::PathBuf> {
        let p = dir.join(name);
        std::fs::write(&p, body)?;
        Ok(p)
    };

    // Header, comments and blank lines are skipped; relative paths resolve
    // against the manifest's directory.
    let p = write(
        "ok.tsv",
        "sample\tfragments\tcells\n# a comment\n\nS1\tf.tsv.gz\tc.txt\nS2\tf.tsv.gz\tc.txt\n",
    )?;
    let rows = read_manifest(&p)?;
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].sample, "S1");
    assert_eq!(rows[0].fragments, dir.join("f.tsv.gz"));
    assert_eq!(rows[1].cells, dir.join("c.txt"));

    // A header is optional.
    let p = write("noheader.tsv", "S1\tf.tsv.gz\tc.txt\n")?;
    assert_eq!(read_manifest(&p)?.len(), 1);

    // An absolute path is taken as given, not joined to the manifest dir.
    let abs = dir.join("f.tsv.gz");
    let p = write(
        "abs.tsv",
        &format!("S1\t{}\t{}\n", abs.display(), abs.display()),
    )?;
    assert_eq!(read_manifest(&p)?[0].fragments, abs);

    // A row named "sample" that is not the first row is data, not a header.
    let p = write(
        "late.tsv",
        "S1\tf.tsv.gz\tc.txt\nsample\tfragments\tcells\n",
    )?;
    assert_eq!(read_manifest(&p)?.len(), 2);

    for (name, body, want) in [
        ("wrongcols.tsv", "S1\tf.tsv.gz\n", "must have 3"),
        ("empty.tsv", "# only a comment\n", "no rows"),
        ("emptyfield.tsv", "S1\t\tc.txt\n", "empty field"),
        (
            "dupsample.tsv",
            "S1\tf.tsv.gz\tc.txt\nS1\tf.tsv.gz\tc.txt\n",
            "must be unique",
        ),
    ] {
        let p = write(name, body)?;
        let err = read_manifest(&p).expect_err(name);
        assert!(
            err.to_string().contains(want),
            "{name}: expected {want:?}, got {}",
            err
        );
    }

    std::fs::remove_dir_all(&dir)?;
    Ok(())
}

#[test]
fn test_cells_file_parsing() -> io::Result<()> {
    let dir = temp_dir("cells_parsing")?;
    let write = |name: &str, body: &str| -> io::Result<std::path::PathBuf> {
        let p = dir.join(name);
        std::fs::write(&p, body)?;
        Ok(p)
    };

    // One column: the output name is the barcode. Order is preserved.
    let (cells, names) = read_cells_file(&write("one.txt", "b\na\nc\n")?)?;
    assert_eq!(names, vec!["b", "a", "c"]);
    assert_eq!(cells["b"], 0);
    assert_eq!(cells["a"], 1);
    assert_eq!(cells["c"], 2);

    // Two columns: barcode keyed by column 1, output name from column 2.
    let (cells, names) = read_cells_file(&write("two.txt", "bc1\tout1\nbc2\tout2\n")?)?;
    assert_eq!(names, vec!["out1", "out2"]);
    assert_eq!(cells["bc1"], 0);
    assert_eq!(cells["bc2"], 1);
    assert!(!cells.contains_key("out1"), "lookup must be by barcode");

    // An empty second column falls back to the barcode.
    let (_, names) = read_cells_file(&write("emptycol2.txt", "bc1\t\n")?)?;
    assert_eq!(names, vec!["bc1"]);

    // A line repeated verbatim is deduplicated, as it always has been.
    let (_, names) = read_cells_file(&write("dup.txt", "a\na\nb\n")?)?;
    assert_eq!(names, vec!["a", "b"]);
    let (_, names) = read_cells_file(&write("dup2.txt", "a\tx\na\tx\n")?)?;
    assert_eq!(names, vec!["x"]);

    // Blank lines are skipped.
    let (_, names) = read_cells_file(&write("blank.txt", "a\n\nb\n")?)?;
    assert_eq!(names, vec!["a", "b"]);

    for (name, body, want) in [
        ("empty.txt", "", "empty"),
        ("threecol.txt", "a\tb\tc\n", "more than 2"),
        ("nobarcode.txt", "\tname\n", "empty barcode"),
        // One barcode cannot name two output cells: one would be dropped.
        (
            "twonames.txt",
            "a\tx\na\ty\n",
            "two different output cell names",
        ),
        // Nor can two barcodes share an output cell.
        ("twobarcodes.txt", "a\tx\nb\tx\n", "two different barcodes"),
    ] {
        let p = write(name, body)?;
        let err = read_cells_file(&p).expect_err(name);
        assert!(
            err.to_string().contains(want),
            "{name}: expected {want:?}, got {}",
            err
        );
    }

    std::fs::remove_dir_all(&dir)?;
    Ok(())
}

// The 2-column cells file is also available to single-file runs, where it
// renames output columns without any manifest.
#[test]
fn test_single_file_two_column_cells_file_renames_output() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let (index, fragments, cells) = fixture(
        "two_col_cells_single",
        &idx,
        "chr1\t99\t200\tbc1\n",
        "bc1\trenamed\n",
    )?;

    for threads in [1usize, 4] {
        let (m, cell_order, _, _) = fcount_core(&fragments, &index, &cells, threads, true)?;
        assert_eq!(cell_order, vec!["renamed"]);
        assert_eq!(m[[0, 0]], 1, "-t {threads}: site");
        assert_eq!(m[[0, 1]], 1, "-t {threads}: control");
    }

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Counting
// ---------------------------------------------------------------------------

// The decode and counting pools run at the same time, so they have to come out
// of one --threads budget instead of being summed.
#[test]
fn test_thread_budget_is_not_additive() {
    if std::env::var_os("BACHA_BGZF_WORKERS").is_some() {
        // An explicit override is a documented tuning escape hatch.
        return;
    }
    for total in [2usize, 3, 4, 6, 12, 64] {
        let b = thread_budget(total, true);
        assert_eq!(
            b.bgzf + b.counting,
            total,
            "-t {total}: the split must add up to the request, not exceed it"
        );
        assert!(
            b.counting >= 1,
            "-t {total}: a counting worker must survive"
        );
        assert!(
            b.bgzf <= DEFAULT_BGZF_WORKERS,
            "-t {total}: decode share is capped"
        );
    }

    // The shipped default.
    assert_eq!(
        (thread_budget(6, true).bgzf, thread_budget(6, true).counting),
        (2, 4)
    );
    // Too small to give decode its full share without starving counting.
    assert_eq!(
        (thread_budget(2, true).bgzf, thread_budget(2, true).counting),
        (1, 1)
    );
    // Plain gzip is one deflate stream inflated inline by the reader thread, so
    // there is no pool to fund and counting gets everything.
    assert_eq!(
        (
            thread_budget(6, false).bgzf,
            thread_budget(6, false).counting
        ),
        (0, 6)
    );
    // -t 1 takes the single-threaded path, which has no decode pool at all.
    assert_eq!(
        (thread_budget(1, true).bgzf, thread_budget(1, true).counting),
        (0, 1)
    );
    // 0 is rejected by `precheck`, but must never produce a zero-thread run.
    assert_eq!(thread_budget(0, true).counting, 1);
}

// A chromosome the index knows about but has no intervals on -- one excluded
// when the index was built, like the chrM/chrY/chrX that `BachaIndex()` drops
// by default -- is a known chromosome, not a genome mismatch. Its fragments are
// skipped silently instead of counting toward the unmatched-chromosome rate,
// which used to warn (and, when every fragment was on such a chromosome, fail).
#[test]
fn test_interval_free_chromosome_in_table_is_not_a_mismatch() -> io::Result<()> {
    let idx = TestIndex {
        window: 1,
        tf_names: vec!["TF1".to_string()],
        // chrM is in the chromosome table but carries no site or control.
        chroms: vec![("chr1".to_string(), 1000), ("chrM".to_string(), 16569)],
        sites: vec![(0, 100, 1)],
        controls: vec![(0, 200, vec![2])],
        total_indices_override: None,
    };

    // Every fragment on the interval-free chromosome: previously indistinguishable
    // from an assembly mismatch, so the run aborted.
    let all_chrm = "chrM\t99\t200\tcell1\n\
                    chrM\t299\t400\tcell1\n";
    let (index, fragments, cells) = fixture("chrm_only", &idx, all_chrm, "cell1\n")?;
    for threads in [1usize, 4usize] {
        let (matrix, _, _, _) = fcount_core(&fragments, &index, &cells, threads, true)
            .unwrap_or_else(|e| panic!("-t {threads}: chrM is a known chromosome: {e}"));
        assert_eq!(matrix[[0, 0]], 0, "-t {threads}: site counts");
        assert_eq!(matrix[[0, 1]], 0, "-t {threads}: control counts");
    }
    std::fs::remove_dir_all(index.parent().unwrap())?;

    // Mixed input: the chr1 fragments still count, the chrM ones are ignored.
    let mixed = "chr1\t99\t200\tcell1\n\
                 chrM\t99\t200\tcell1\n\
                 chrM\t299\t400\tcell1\n";
    let (index, fragments, cells) = fixture("chrm_mixed", &idx, mixed, "cell1\n")?;
    for threads in [1usize, 4usize] {
        let (matrix, _, _, _) = fcount_core(&fragments, &index, &cells, threads, true)
            .unwrap_or_else(|e| panic!("-t {threads}: should not fail: {e}"));
        assert_eq!(matrix[[0, 0]], 1, "-t {threads}: site counts");
        assert_eq!(matrix[[0, 1]], 1, "-t {threads}: control counts");
    }
    std::fs::remove_dir_all(index.parent().unwrap())?;

    // A chromosome absent from the table entirely is still a mismatch.
    let unknown = "contigX\t99\t200\tcell1\n";
    let (index, fragments, cells) = fixture("chrm_unknown", &idx, unknown, "cell1\n")?;
    let err = fcount_core(&fragments, &index, &cells, 1, true)
        .expect_err("a chromosome outside the index table is a mismatch");
    assert!(
        err.to_string()
            .contains("No fragments matched any chromosome"),
        "unexpected error: {err}"
    );
    std::fs::remove_dir_all(index.parent().unwrap())?;

    Ok(())
}

// `-t 0` reaches the gzip writer only after the whole fragment file has been
// counted, so it has to be rejected up front: otherwise a long run is thrown
// away at the last moment and leaves an empty output file that blocks the retry.
#[test]
fn test_zero_threads_is_rejected_before_counting() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let (index, fragments, cells) =
        fixture("zero_threads", &idx, "chr1\t99\t200\tcell1\n", "cell1\n")?;
    let out = index.parent().unwrap().join("scores.csv.gz");

    let err = bacha(&fragments, &index, &cells, &out, 0, 10.0, 1.0, true)
        .expect_err("-t 0 must be rejected");
    assert!(
        err.to_string().contains("threads must be >= 1"),
        "unexpected error: {err}"
    );
    // Nothing was written, so re-running with a valid -t is not blocked by the
    // "output file already exists" check.
    assert!(
        !out.exists(),
        "a rejected run must not leave an output file"
    );

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

#[test]
fn test_counting_treats_index_coordinates_as_one_based() -> io::Result<()> {
    // Site at 1-based chr1:10 => 0-based [9, 10); control at 1-based 20 =>
    // 0-based [19, 20). The fragment's insertions are at `start` and `end - 1`,
    // so 9..20 puts one on each.
    let idx = TestIndex::simple(1, 100, 10, 20);
    let (index, fragments, cells) =
        fixture("one_based_index", &idx, "chr1\t9\t20\tcell1\n", "cell1\n")?;

    let (count_matrix, cell_order, tf_list, n_instances) =
        fcount_core(&fragments, &index, &cells, 1, true)?;

    assert_eq!(cell_order, vec!["cell1"]);
    assert_eq!(tf_list, vec!["TF1"]);
    assert_eq!(n_instances, vec![1, 1]);
    assert_eq!(count_matrix.shape(), &[1, 2]);
    assert_eq!(count_matrix[[0, 0]], 1, "site feature");
    assert_eq!(count_matrix[[0, 1]], 1, "control feature");

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

#[test]
fn test_sharded_counting_matches_single_threaded() -> io::Result<()> {
    // 10 cells across 2 chromosomes. With 4 workers most cells are owned by a
    // different worker than the one parsing their chunk, exercising the
    // cross-worker forwarding path.
    let idx = TestIndex::simple(2, 1000, 100, 200);
    let n_cells = 10usize;
    let mut frag_lines = String::new();
    let mut cell_lines = String::new();
    for i in 0..n_cells {
        let chrom = if i % 2 == 0 { "chr1" } else { "chr2" };
        frag_lines.push_str(&format!("{chrom}\t99\t200\tcell{i}\n"));
        cell_lines.push_str(&format!("cell{i}\n"));
    }
    let (index, fragments, cells) = fixture("sharded_parity", &idx, &frag_lines, &cell_lines)?;

    let (single, cells_s, tfs_s, _) = fcount_core(&fragments, &index, &cells, 1, true)?;
    let (sharded, cells_p, tfs_p, _) = fcount_core(&fragments, &index, &cells, 4, true)?;

    assert_eq!(cells_s, cells_p);
    assert_eq!(tfs_s, tfs_p);
    assert_eq!(
        single, sharded,
        "sharded (n=4) output must match the single-threaded reference"
    );
    assert_eq!(single.shape(), &[n_cells, 2]);
    for row in 0..n_cells {
        assert_eq!(single[[row, 0]], 1, "cell {row} site count");
        assert_eq!(single[[row, 1]], 1, "cell {row} control count");
    }

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

// Paired insertion counting: a fragment adds at most one count to any interval,
// however many of its insertions land inside. Every other counting test uses
// 1 bp intervals, where two insertions can never share an interval, so none of
// them can distinguish this from counting the insertions independently.
#[test]
fn test_fragment_counted_once_per_interval() -> io::Result<()> {
    // Site at 1-based [1001, 1100] => 0-based [1000, 1100); control at
    // 1-based [2001, 2100] => 0-based [2000, 2100).
    let idx = TestIndex {
        window: 100,
        tf_names: vec!["TF1".to_string()],
        chroms: vec![("chr1".to_string(), 10_000)],
        sites: vec![(0, 1001, 1)],
        controls: vec![(0, 2001, vec![2])],
        total_indices_override: None,
    };

    // (label, fragment, expected site count, expected control count).
    // Insertions are at `start` and `end - 1`.
    let cases: [(&str, &str, u32, u32); 5] = [
        // Both insertions inside one site window: one molecule, one count.
        (
            "both insertions in the site",
            "chr1\t1010\t1081\tcell1\n",
            1,
            0,
        ),
        // Degenerate 1 bp fragment: both insertions on the same base.
        (
            "both insertions on one base",
            "chr1\t1010\t1011\tcell1\n",
            1,
            0,
        ),
        ("only ins1 in the site", "chr1\t1050\t1501\tcell1\n", 1, 0),
        ("only ins2 in the site", "chr1\t500\t1051\tcell1\n", 1, 0),
        // Distinct intervals each get their own count -- PIC caps per interval,
        // not per fragment.
        ("one insertion in each", "chr1\t1050\t2051\tcell1\n", 1, 1),
    ];

    for (i, (label, frag, want_site, want_ctrl)) in cases.iter().enumerate() {
        let (index, fragments, cells) = fixture(&format!("pic_{i}"), &idx, frag, "cell1\n")?;
        for threads in [1usize, 4] {
            let (m, _, _, _) = fcount_core(&fragments, &index, &cells, threads, true)?;
            assert_eq!(m[[0, 0]], *want_site, "{label} (-t {threads}): site");
            assert_eq!(m[[0, 1]], *want_ctrl, "{label} (-t {threads}): control");
        }
        std::fs::remove_dir_all(index.parent().unwrap())?;
    }

    Ok(())
}

// A control region shared by several features credits each of them, and still
// only once per fragment.
#[test]
fn test_shared_control_counted_once_per_feature_per_fragment() -> io::Result<()> {
    let idx = TestIndex {
        window: 100,
        tf_names: vec!["TF1".to_string(), "TF2".to_string()],
        chroms: vec![("chr1".to_string(), 10_000)],
        sites: vec![(0, 1001, 1), (0, 3001, 2)],
        // One control region selected by both features (1-based indices
        // n_tfs + 1 and n_tfs + 2).
        controls: vec![(0, 2001, vec![3, 4])],
        total_indices_override: None,
    };
    // Both insertions inside the shared control window.
    let (index, fragments, cells) = fixture(
        "shared_control_pic",
        &idx,
        "chr1\t2010\t2081\tcell1\n",
        "cell1\n",
    )?;

    for threads in [1usize, 4] {
        let (m, _, tf_list, _) = fcount_core(&fragments, &index, &cells, threads, true)?;
        assert_eq!(tf_list.len(), 2);
        assert_eq!(m[[0, 0]], 0, "-t {threads}: TF1 site");
        assert_eq!(m[[0, 1]], 0, "-t {threads}: TF2 site");
        assert_eq!(m[[0, 2]], 1, "-t {threads}: TF1 control");
        assert_eq!(m[[0, 3]], 1, "-t {threads}: TF2 control");
    }

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

#[test]
fn test_second_insertion_is_at_end_minus_one() -> io::Result<()> {
    // Site at 1-based [1001, 1100] => 0-based [1000, 1100). Its last base is
    // 0-based 1099, so a fragment ending at 1100 has its second insertion
    // exactly on that base, and one ending at 1101 falls just past the window.
    let idx = TestIndex {
        window: 100,
        tf_names: vec!["TF1".to_string()],
        chroms: vec![("chr1".to_string(), 10_000)],
        sites: vec![(0, 1001, 1)],
        controls: vec![(0, 5001, vec![2])],
        total_indices_override: None,
    };

    // Each fragment starts far from the site, so only its second insertion can
    // land in the window and the count isolates that endpoint.
    let cases: [(&str, &str, u32); 3] = [
        ("last base of the window", "chr1\t500\t1100\tcell1\n", 1),
        ("one past the window", "chr1\t500\t1101\tcell1\n", 0),
        ("first base of the window", "chr1\t500\t1001\tcell1\n", 1),
    ];

    for (i, (label, frag, want)) in cases.iter().enumerate() {
        let (index, fragments, cells) =
            fixture(&format!("end_insertion_{i}"), &idx, frag, "cell1\n")?;
        for threads in [1usize, 4] {
            let (m, _, _, _) = fcount_core(&fragments, &index, &cells, threads, true)?;
            assert_eq!(m[[0, 0]], *want, "{label} (-t {threads})");
        }
        std::fs::remove_dir_all(index.parent().unwrap())?;
    }

    Ok(())
}

// A truncated BGZF fragment file
#[test]
fn test_truncated_bgzf_fragment_file_is_rejected() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let dir = temp_dir("truncated_bgzf")?;
    let index_path = dir.join("index.bi");
    let cell_path = dir.join("cells.txt");
    idx.write(&index_path)?;
    std::fs::write(&cell_path, "cell1\n")?;

    let frags = "chr1\t99\t199\tcell1\nchr1\t99\t199\tcell1\n";
    let complete = dir.join("complete.tsv.gz");
    let truncated = dir.join("truncated.tsv.gz");
    write_bgzf_text(&complete, frags, true)?;
    write_bgzf_text(&truncated, frags, false)?;

    for threads in [1usize, 4] {
        // The complete file still counts normally.
        let (matrix, _, _, _) = fcount_core(&complete, &index_path, &cell_path, threads, true)?;
        assert_eq!(matrix[[0, 0]], 2, "-t {threads}: complete file site counts");

        // The truncated one is refused
        let err = fcount_core(&truncated, &index_path, &cell_path, threads, true)
            .expect_err("truncated BGZF must be rejected");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof, "-t {threads}");
        let msg = err.to_string();
        assert!(msg.contains("truncated"), "-t {threads}: {msg}");
    }

    std::fs::remove_dir_all(&dir)?;
    Ok(())
}

// A fragment file whose last line has no terminating newline is counted the
// same way regardless of thread count
#[test]
fn test_final_line_without_newline_is_counted_on_every_path() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let (index, fragments, cells) = fixture(
        "unterminated_final_line",
        &idx,
        "chr1\t99\t199\tcell1\nchr1\t99\t199\tcell1",
        "cell1\n",
    )?;

    let (single, _, _, _) = fcount_core(&fragments, &index, &cells, 1, true)?;
    let (sharded, _, _, _) = fcount_core(&fragments, &index, &cells, 4, true)?;

    assert_eq!(
        single[[0, 0]],
        2,
        "single-threaded must count the final line"
    );
    assert_eq!(single, sharded, "thread count must not change the counts");

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

#[test]
fn test_bounds_check_aborts_on_out_of_bounds_coordinate() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let (index, fragments, cells) = fixture(
        "bounds_check",
        &idx,
        "chr1\t50\t150\tcell1\nchr1\t900\t5000\tcell1\n",
        "cell1\n",
    )?;

    for threads in [1usize, 4usize] {
        let err = fcount_core(&fragments, &index, &cells, threads, true).expect_err(&format!(
            "expected a genome mismatch error with -t {threads}"
        ));
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput, "-t {threads}");
        let msg = err.to_string();
        assert!(msg.contains("chr1"), "-t {threads}: {msg}");
        assert!(msg.contains("5000"), "-t {threads}: {msg}");
    }

    for threads in [1usize, 4usize] {
        let (matrix, _, _, _) = fcount_core(&fragments, &index, &cells, threads, false)
            .unwrap_or_else(|e| panic!("-t {threads}: --no-bounds-check should not error: {e}"));
        assert_eq!(matrix.shape(), &[1, 2], "-t {threads}");
    }

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Input-overlap checks
// ---------------------------------------------------------------------------

// A fragment file whose chromosome names never appear in the index (the
// classic 'chr1' vs '1' mismatch) must be a hard error. When the unmatched
// count was folded into its own denominator this reported "50.0%" and
// returned an all-zero score matrix with exit status 0.
#[test]
fn test_errors_when_no_fragment_matches_a_chromosome() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let frags: String = (1..=5)
        .map(|i| format!("1\t{}\t{}\tcell1\n", i * 100, i * 100 + 50))
        .collect();
    let (index, fragments, cells) = fixture("no_chrom_match", &idx, &frags, "cell1\n")?;

    for threads in [1usize, 4usize] {
        let err = fcount_core(&fragments, &index, &cells, threads, true).expect_err(&format!(
            "expected a genome mismatch error with -t {threads}"
        ));
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput, "-t {threads}");
        assert!(
            err.to_string()
                .contains("No fragments matched any chromosome"),
            "-t {threads}: {err}"
        );
    }

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

// A partial chromosome mismatch stays a warning: fragments that do match must
// still be counted rather than aborting the run.
#[test]
fn test_partial_chromosome_mismatch_still_counts_matching_fragments() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    // 2 fragments on chr1 (each hitting site and control, via insertions at
    // `start` and `end - 1`) and 3 on an unknown contig => 60% unmatched, above
    // the 20% warning threshold.
    let frags = "chr1\t99\t200\tcell1\n\
                 chr1\t99\t200\tcell1\n\
                 contigX\t99\t200\tcell1\n\
                 contigX\t99\t200\tcell1\n\
                 contigX\t99\t200\tcell1\n";
    let (index, fragments, cells) = fixture("partial_chrom_match", &idx, frags, "cell1\n")?;

    for threads in [1usize, 4usize] {
        let (matrix, _, _, _) = fcount_core(&fragments, &index, &cells, threads, true)
            .unwrap_or_else(|e| panic!("-t {threads}: should warn, not fail: {e}"));
        assert_eq!(matrix[[0, 0]], 2, "-t {threads}: site counts");
        assert_eq!(matrix[[0, 1]], 2, "-t {threads}: control counts");
    }

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

// The sibling failure mode: every barcode in the fragment file is absent from
// the cell file (typically a '-1' suffix mismatch). Without a check this also
// yields a clean all-zero matrix and exit status 0.
#[test]
fn test_errors_when_no_fragment_matches_a_cell_barcode() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let frags = "chr1\t99\t199\tcellA-1\n\
                 chr1\t99\t199\tcellB-1\n\
                 chr1\t99\t199\tcellC-1\n";
    let (index, fragments, cells) =
        fixture("no_barcode_match", &idx, frags, "cellA\ncellB\ncellC\n")?;

    for threads in [1usize, 4usize] {
        let err = fcount_core(&fragments, &index, &cells, threads, true).expect_err(&format!(
            "expected a barcode mismatch error with -t {threads}"
        ));
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput, "-t {threads}");
        assert!(
            err.to_string()
                .contains("No fragments matched any cell barcode"),
            "-t {threads}: {err}"
        );
    }

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Index validation
// ---------------------------------------------------------------------------

#[test]
fn test_text_index_is_rejected() -> io::Result<()> {
    let dir = temp_dir("text_index_rejected")?;
    let index = dir.join("legacy.bacha.gz");
    let fragments = dir.join("fragments.tsv.gz");
    let cells = dir.join("cells.txt");
    write_gzip_text(
        &index,
        "##[sites]\nTF1\t1\t1\n##[controls]\nTF1\t2\t1\n\
         ##[site positions]\nchr1\t100\t100\t1\n\
         ##[control positions]\nchr1\t200\t200\t2\n",
    )?;
    write_gzip_text(&fragments, "chr1\t99\t199\tcell1\n")?;
    std::fs::write(&cells, "cell1\n")?;

    let err = fcount_core(&fragments, &index, &cells, 1, true)
        .expect_err("a plain-text index must be rejected");
    assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    assert!(
        err.to_string().contains("BACH"),
        "the error should name the missing magic bytes, got: {err}"
    );

    std::fs::remove_dir_all(&dir)?;
    Ok(())
}

// Malformed binary indexes must produce a clear error at load time. Each case
// below would otherwise either panic (out-of-range array indices) or silently
// mis-attribute counts to the wrong feature column.
#[test]
fn test_malformed_index_is_rejected_at_load() -> io::Result<()> {
    let dir = temp_dir("malformed_index")?;
    let fragments = dir.join("fragments.tsv.gz");
    let cells = dir.join("cells.txt");
    write_gzip_text(&fragments, "chr1\t99\t199\tcell1\n")?;
    std::fs::write(&cells, "cell1\n")?;

    let base = TestIndex::simple(1, 1000, 100, 200);

    let mut no_features = base.clone();
    no_features.tf_names.clear();

    let mut site_index_zero = base.clone();
    site_index_zero.sites[0].2 = 0;

    // A site index reaching into the control half of the shared ID space
    // would have quietly incremented a control column.
    let mut site_in_control_range = base.clone();
    site_in_control_range.sites[0].2 = 2;

    let mut control_in_site_range = base.clone();
    control_in_site_range.controls[0].2 = vec![1];

    // Chromosome IDs are used as direct vector indices.
    let mut chrom_id_out_of_range = base.clone();
    chrom_id_out_of_range.sites[0].0 = 7;

    // Per-row counts disagreeing with the declared total sends the within-row
    // cumsum off the end of its buffer.
    let mut inconsistent_totals = base.clone();
    inconsistent_totals.total_indices_override = Some(5);

    let cases: [(&str, &TestIndex, &str); 6] = [
        ("no_features", &no_features, "no features defined"),
        ("site_index_zero", &site_index_zero, "site feature index"),
        (
            "site_in_control_range",
            &site_in_control_range,
            "site feature index",
        ),
        (
            "control_in_site_range",
            &control_in_site_range,
            "control feature index",
        ),
        (
            "chrom_id_out_of_range",
            &chrom_id_out_of_range,
            "chromosome ID",
        ),
        (
            "inconsistent_totals",
            &inconsistent_totals,
            "feature index/indices in total",
        ),
    ];

    for (name, idx, expected) in cases {
        let index = dir.join(format!("{name}.bi"));
        idx.write(&index)?;
        let err = fcount_core(&fragments, &index, &cells, 1, true)
            .expect_err(&format!("{name}: expected the index to be rejected"));
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "{name}: {err}");
        assert!(
            err.to_string().contains(expected),
            "{name}: message {:?} does not mention {expected:?}",
            err.to_string()
        );
    }

    std::fs::remove_dir_all(&dir)?;
    Ok(())
}

// An existing output path must be reported as an error the caller can handle,
// not terminate the process from inside the library. It must also leave the
// existing file untouched.
#[test]
fn test_existing_output_file_is_an_error() -> io::Result<()> {
    let idx = TestIndex::simple(1, 1000, 100, 200);
    let (index, fragments, cells) =
        fixture("existing_output", &idx, "chr1\t99\t199\tcell1\n", "cell1\n")?;
    let output = index.parent().unwrap().join("scores.csv");
    std::fs::write(&output, "pre-existing contents\n")?;

    let err = bacha(&fragments, &index, &cells, &output, 1, 10.0, 1.0, true)
        .expect_err("an existing output file must be rejected");
    let msg = err.to_string();
    assert!(msg.contains("already exists"), "{msg}");
    assert_eq!(
        err.downcast_ref::<io::Error>().map(|e| e.kind()),
        Some(io::ErrorKind::AlreadyExists),
        "should surface as an io::Error the caller can match on: {msg}"
    );
    assert_eq!(
        std::fs::read_to_string(&output)?,
        "pre-existing contents\n",
        "the existing file must not be overwritten"
    );

    // The same run succeeds once the path is free.
    std::fs::remove_file(&output)?;
    bacha(&fragments, &index, &cells, &output, 1, 10.0, 1.0, true)
        .unwrap_or_else(|e| panic!("scoring should succeed with a free output path: {e}"));
    assert!(std::fs::read_to_string(&output)?.starts_with(",cell1"));

    std::fs::remove_dir_all(index.parent().unwrap())?;
    Ok(())
}

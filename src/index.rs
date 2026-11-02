// Build bacha index file
// Rust port of BachaIndex() - GC content + mean Tn5 bias as matching features
// outputs gzip-compressed text file or binary

use bio::io::fasta::IndexedReader;
use hdf5::File as hdf5_file;
use log::info;
use nalgebra::DMatrix;
use ndarray::Array1;
use rand_distr::{Distribution, Normal};
use rust_lapper::{Interval, Lapper};
use rustc_hash::FxHashMap;
use std::collections::HashSet;
use std::{
    error::Error, fs::File, io, io::BufRead, io::BufReader, io::BufWriter, io::Write, path::Path,
};
use xz2::write::XzEncoder;
// for kd-tree
use kiddo::KdTree;
use kiddo::SquaredEuclidean;
// use kiddo::NearestNeighbour;
use flate2::write::GzEncoder;
use flate2::Compression;

use clap::ValueEnum;
#[derive(Clone, Copy, PartialEq, Debug, ValueEnum)]
pub enum IndexFormat {
    Text,
    Binary,
    Both,
}
// Define structs for writing bacha index file
//struct Position {chrom: String, start: u32, end: u32, index: usize,}
//struct IndexEntry {id: String, index: usize, count: usize,}
// Define structs for collecting motif and bg candidate sites
#[derive(Clone)]
pub struct MotifSite {
    pub chrom: String,
    pub start: usize,
    pub stop: usize,
    //    pub tf_id: String,
    pub tf_index: usize,
    pub features: [f32; 2], // GC, mean_bias
}
#[derive(Clone)]
pub struct BgInterval {
    pub chrom: String,
    pub start: usize,
    pub stop: usize,
    pub features: [f32; 2],
}

pub fn build_index(
    motifmatch: &str,
    ocr: &str,
    bias: &str,
    fasta_path: &str,
    genome_name: &str,
    outfile: &str,
    format: IndexFormat,
    motif_set_name: &str,
    k: usize,
    min_dist: usize,
) -> Result<(), Box<dyn Error>> {
    // bedfile: chr, start, end, motif_id, site_integer (full_motif_match.bed)
    let motif_match_file = Path::new(motifmatch)
        .canonicalize()
        .expect("Can't find path to input BED file");
    info!("Received motif match BED file: {:?}", motif_match_file);

    // OCR where the background regions will be picked from ENCODE: GRCh38-cCREs.bed
    // /home/users/astar/gis/limtw/scratch/proj_TF/motif/toolkit/
    let ocr_file = Path::new(ocr)
        .canonicalize()
        .expect("Can't find path to input OCR BED file");
    info!("Received OCR BED file: {:?}", ocr_file);

    //let outfile = outfile;

    let bias_h5 = bias;

    let faidx = load_fasta(fasta_path)?;
    // Some parameters
    let window_size: usize = 50;
    //let k: usize = 10; // num of controls per motif instance
    //let min_dist: usize = 1000; // min gap between a control and any motif site

    // prepare containers
    //let mut site_positions: Vec<Position> = Vec::new();
    //let mut control_positions: Vec<Posiiton> = Vec::new();

    //let mut site_index_table: Vec<IndexEntry> = Vec::new();
    //l/et mut control_index_table: Vec<IndexEntry> = Vec::new();

    //let mut site_index = 1usize;
    //let mut control_index = 1usize;

    // 1. Load motif sites
    let (mut motif_sites, tf_names) = load_motif_sites(&motif_match_file, fasta_path, bias_h5)?;
    let n_tfs = tf_names.len();
    info!(
        "Loaded {} motif sites across {} TFs",
        motif_sites.len(),
        n_tfs
    );

    // 2. Generate candidate background regions
    let mut candidate_bg =
        generate_background_intervals(&ocr_file, fasta_path, bias_h5, window_size)?;
    let n_bg = candidate_bg.len();
    info!("Generated {} background intervals", n_bg);

    // edit: change global exclusion to exclusion per TF motif
    // BUILD AN INTERVAL TREE OVER BACKGROUND,
    // for each TF walk only the TF's sites and mark the controls they cover
    let bg_lapper: FxHashMap<String, Lapper<usize, usize>> = {
        let mut m: FxHashMap<String, Vec<Interval<usize, usize>>> = FxHashMap::default();
        for (i, b) in candidate_bg.iter().enumerate() {
            m.entry(b.chrom.clone()).or_default().push(Interval {
                start: b.start,
                stop: b.stop,
                val: i,
            });
        }
        m.into_iter()
            .map(|(ch, mut ivs)| {
                ivs.sort_by_key(|iv| iv.start);
                (ch, Lapper::new(ivs))
            })
            .collect()
    };
    // group motif row indices by TF (replace the previous per-TF linear scan)
    let mut sites_by_tf: Vec<Vec<usize>> = vec![Vec::new(); n_tfs];
    for (i, s) in motif_sites.iter().enumerate() {
        sites_by_tf[s.tf_index - 1].push(i);
    }

    // add gitter if not it will panick "Too many items with the same position on one axis. Bucket size must be increased to at least 1 more than the number of items with the same position on one axis."
    let gc_jitter_sd = 0.5f32 / 100.0; // cus i am using fraction not %
                                       // thread_rng is often the most convenient source of randomness:
                                       //let mut rng = thread_rng();
    let mut rng = rand::thread_rng();
    let normal = Normal::new(0.0, gc_jitter_sd as f64).unwrap();

    for b in candidate_bg.iter_mut() {
        b.features[0] += normal.sample(&mut rng) as f32;
    }
    for s in motif_sites.iter_mut() {
        s.features[0] += normal.sample(&mut rng) as f32;
    }

    // 4. Whitening Transformation use nalgebra::{DMatrix, Cholesky};
    // build feature matrix
    let d = candidate_bg[0].features.len(); //num_features used -> here its 2 : GC and Bias
    let bg_matrix = DMatrix::from_fn(n_bg, d, |i, j| candidate_bg[i].features[j] as f64); //let bg_mat = DMatrix::from_fn(n, d, |i, j| bg_features[i][j]);
                                                                                          //compute covariance
    let col_mean = bg_matrix.column_mean(); // why is it row_mean and not
    let centered_bg = DMatrix::from_fn(n_bg, d, |i, j| bg_matrix[(i, j)] - col_mean[j]);
    let cov = (&centered_bg.transpose() * &centered_bg) / (n_bg as f64 - 1.0);
    let l_inv = cov
        .cholesky()
        .expect("background covariance matrix is not positive definite")
        .l()
        .try_inverse()
        .expect("Cholesky L is singular");
    // whitened = centered * L^{-T}  (R: t(forwardsolve(t(chol(cov)), t(mat))))
    let bg_white = &centered_bg * l_inv.transpose();
    let motif_matrix = DMatrix::from_fn(motif_sites.len(), 2, |i, j| {
        motif_sites[i].features[j] as f64
    });
    let centered_m = DMatrix::from_fn(motif_sites.len(), 2, |i, j| {
        motif_matrix[(i, j)] - col_mean[j]
    });
    let motif_white = &centered_m * l_inv.transpose();

    // 5. kNN search (equivalent to the RANN in the R code)
    // https://docs.rs/kiddo/latest/kiddo/ claims to be fatest k-d tree library in the world lol
    // good for low-ish numbers of dimensions
    let bg_entries: Vec<[f64; 2]> = (0..n_bg)
        .map(|i| [bg_white[(i, 0)], bg_white[(i, 1)]])
        .collect();
    let kdtree: KdTree<f64, 2> = (&bg_entries).into(); // put in the background
                                                       // Loop and query per motif instance
                                                       // motif_id[bg_idx] accumulates the 1-based TF indices this control is chosen for
    let mut motif_id: Vec<Vec<usize>> = vec![Vec::new(); n_bg]; // collect motif_id for each background site
    let mut n_motif_sites_per_tf: Vec<usize> = vec![0; n_tfs];
    for s in &motif_sites {
        n_motif_sites_per_tf[s.tf_index - 1] += 1;
    }
    for tf in 1..=n_tfs {
        // TODO: parallelize this
        // collect all motif sites belonging to this tf
        let instances = &sites_by_tf[tf - 1];
        if instances.is_empty() {
            continue;
        }
        // per-TF exclusion: controls within min_dist of the TF's sites only
        let mut excluded = vec![false; n_bg];
        for &mi in instances {
            let s = &motif_sites[mi];
            if let Some(l) = bg_lapper.get(&s.chrom) {
                let lo = s.start.saturating_sub(min_dist);
                let hi = s.stop + min_dist;
                for iv in l.find(lo, hi) {
                    excluded[iv.val] = true;
                }
            }
        }
        let n_eligible_total = n_bg - excluded.iter().filter(|&&e| e).count();
        if n_eligible_total == 0 {
            eprintln!(
                "warning: no eligible control regions for TF {}",
                tf_names[tf - 1]
            );
            continue;
        }
        let n_cap = (n_eligible_total as f64 * 0.2) as usize;
        let k_use = k
            .min(n_eligible_total)
            .min((n_cap / instances.len().max(1)).max(1));
        let k_query = ((2 * k_use).max(20)).min(n_eligible_total);
        let mut used = excluded; // moved instead of cloned

        for &mi in instances {
            let query = [motif_white[(mi, 0)], motif_white[(mi, 1)]];
            let mut assigned = 0usize;
            let mut fetch = k_query;
            loop {
                let neighbors = kdtree.nearest_n::<SquaredEuclidean>(&query, fetch);
                let got = neighbors.len();
                for nn in neighbors {
                    let bi = nn.item as usize;
                    if used[bi] {
                        continue;
                    }
                    motif_id[bi].push(tf); // store 1-based TF index
                    used[bi] = true;
                    assigned += 1;
                    if assigned >= k_use {
                        break;
                    }
                }
                if assigned >= k_use || fetch >= n_bg || got < fetch {
                    break;
                }
                fetch = (fetch * 4).min(n_bg);
            }
        }
    }

    //drop bg intervals not selected
    let (controls, ctrl_motif_id): (Vec<BgInterval>, Vec<Vec<usize>>) = candidate_bg
        .into_iter()
        .zip(motif_id.into_iter())
        .filter(|(_, ids)| !ids.is_empty())
        .unzip();

    let mut n_controls_per_tf: Vec<usize> = vec![0; n_tfs];
    for ids in &ctrl_motif_id {
        for &id in ids {
            n_controls_per_tf[id - 1] += 1;
        }
    }

    // Compute the chromosome lengths in faidx  (actually only needed for binary format)
    let fai_lengths: Vec<(String, u32)> = faidx
        .index
        .sequences()
        .iter()
        .map(|s| {
            let name = format!("chr{}", s.name.strip_prefix("chr").unwrap_or(&s.name));
            (name, s.len as u32)
        })
        .collect();

    // let user choose from indexformat
    match format {
        IndexFormat::Text => write_index_text(
            &format!("{outfile}.txt.gz"),
            &motif_sites,
            &controls,
            &ctrl_motif_id,
            &tf_names,
            &n_motif_sites_per_tf,
            &n_controls_per_tf,
            n_tfs,
        )?,
        IndexFormat::Binary => write_index_binary(
            &format!("{outfile}.bacha.xz"),
            &motif_sites,
            &controls,
            &ctrl_motif_id,
            &tf_names,
            &n_motif_sites_per_tf,
            &n_controls_per_tf,
            &fai_lengths,
            window_size,
            genome_name,
            motif_set_name,
        )?,
        IndexFormat::Both => {
            write_index_text(
                &format!("{outfile}.txt.gz"),
                &motif_sites,
                &controls,
                &ctrl_motif_id,
                &tf_names,
                &n_motif_sites_per_tf,
                &n_controls_per_tf,
                n_tfs,
            )?;
            write_index_binary(
                &format!("{outfile}.bacha.xz"),
                &motif_sites,
                &controls,
                &ctrl_motif_id,
                &tf_names,
                &n_motif_sites_per_tf,
                &n_controls_per_tf,
                &fai_lengths,
                window_size,
                genome_name,
                motif_set_name,
            )?;
        }
    }

    eprintln!("Done. Index written to {}", outfile);

    Ok(())
}

fn load_motif_sites(
    motif_match_file: &Path,
    ref_path: &str,
    bias_h5: &str,
) -> io::Result<(Vec<MotifSite>, Vec<String>)> {
    //let mut flank_properties_vector : FxHashMap<usize, Vec<Vec<f32>>> = FxHashMap::default(); -> changed to be part of MotifSite
    let mut sites: Vec<MotifSite> = Vec::new();
    let mut tf_index_map: FxHashMap<String, usize> = FxHashMap::default();
    let mut tf_names: Vec<String> = Vec::new();

    // read predicted bias
    let bias_factors = read_bias_factors(bias_h5);
    // read genome fasta // need to put this as argument
    let mut faidx = load_fasta(ref_path).expect("Could not load ref fasta");
    let mut seq = Vec::new();

    // Read motif matches & compute properties // put as a separate function?
    let file = File::open(motif_match_file)?;
    let reader = BufReader::new(file);
    for (index, line) in reader.lines().enumerate() {
        let line = line.expect("Error reading line");
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }
        let chr = fields[0].to_string();
        let chr_trimmed = &chr[3..];
        let motif_start: usize = fields[1]
            .parse()
            .expect(&format!("Invalid start position {}", index + 1));
        let motif_end: usize = fields[2]
            .parse()
            .expect(&format!("Invalid end position {}", index + 1)); // non inclusive
        let tf_id = fields[3].to_string();
        // need to make a tf_index
        let tf_index = if let Some(&idx) = tf_index_map.get(&tf_id) {
            idx
        } else {
            let idx = tf_names.len() + 1; // 1-based to match the R version
            tf_index_map.insert(tf_id.clone(), idx);
            tf_names.push(tf_id.clone());
            idx
        };

        //// match_features = c("GC", "meanBias") // can change later to be user-customized
        let gc_content: f32 = {
            seq.clear();
            faidx
                .fetch(&chr_trimmed, motif_start as u64, motif_end as u64)
                .expect("Couldn't fetch interval");
            faidx.read(&mut seq).expect("Couldnt read the interval");
            gc_percent(&seq)
        };
        let mean_bias: f32 = cal_mean_bias(&chr, motif_start, motif_end, &bias_factors);

        sites.push(MotifSite {
            chrom: chr,
            start: motif_start,
            stop: motif_end,
            tf_index,
            features: [gc_content, mean_bias],
        });
    }
    Ok((sites, tf_names))
}

fn write_index_text(
    path: &str,
    sites: &[MotifSite],
    controls: &[BgInterval],
    motif_id: &[Vec<usize>],
    tf_names: &[String],
    n_motif: &[usize],
    n_controls: &[usize],
    n_tfs: usize,
) -> std::io::Result<()> {
    let file = File::create(path)?;
    let gz = GzEncoder::new(file, Compression::default());
    let mut w = BufWriter::new(gz);

    // write sites
    writeln!(w, "##[sites]")?;
    writeln!(w, "# id\tindex\tcount")?;
    for (i, name) in tf_names.iter().enumerate() {
        writeln!(w, "{}\t{}\t{}", name, i + 1, n_motif[i])?;
    }
    // write controls
    writeln!(w, "##[controls]")?;
    writeln!(w, "# id\tindex\tcount")?;
    for (i, name) in tf_names.iter().enumerate() {
        writeln!(w, "{}\t{}\t{}", name, i + 1 + n_tfs, n_controls[i])?;
    }
    // write site positions
    writeln!(w, "##[site positions]")?;
    writeln!(w, "# chrom\tstart\tend\tindex")?;
    let mut site_order: Vec<usize> = (0..sites.len()).collect();
    site_order.sort_by(|&a, &b| {
        sites[a]
            .chrom
            .cmp(&sites[b].chrom)
            .then(sites[a].start.cmp(&sites[b].start))
    });
    for &si in &site_order {
        let s = &sites[si];
        writeln!(w, "{}\t{}\t{}\t{}", s.chrom, s.start, s.stop, s.tf_index)?;
    }
    // write control positions
    writeln!(w, "##[control positions]")?;
    writeln!(w, "# chrom\tstart\tend\tindex")?;
    let mut ctrl_order: Vec<usize> = (0..controls.len()).collect();
    ctrl_order.sort_by(|&a, &b| {
        controls[a]
            .chrom
            .cmp(&controls[b].chrom)
            .then(controls[a].start.cmp(&controls[b].start))
    });
    for &ci in &ctrl_order {
        let c = &controls[ci];
        let ids: Vec<String> = motif_id[ci]
            .iter()
            .map(|x| (x + n_tfs).to_string())
            .collect();
        writeln!(w, "{}\t{}\t{}\t{}", c.chrom, c.start, c.stop, ids.join(","))?;
    }

    // explicitly flush and finish:
    w.into_inner().map_err(|e| e.into_error())?.finish()?;
    Ok(())
}

// Binary writer
fn write_str<W: Write>(w: &mut W, s: &str) -> io::Result<()> {
    // there is no \n \t in binary
    // so each string is written as its length, then its bytes
    // reader does: read 2 bytes -> get n -> read exactly n bytes -> get the string
    let b = s.as_bytes();
    w.write_all(&(b.len() as u16).to_le_bytes())?;
    w.write_all(b)
}
fn delta_per_chrom(chroms: &[u8], starts: &[u32]) -> Vec<u32> {
    // positions sorted by chrom, then sorted by start
    // instead of storing the absolute position, we just store the difference
    // cus diff is much smaller than the absolute
    let mut out = starts.to_vec();
    for i in (1..starts.len()).rev() {
        if chroms[i] == chroms[i - 1] {
            out[i] = starts[i] - starts[i - 1]
        }
    }
    out
}
fn write_index_binary(
    path: &str,
    sites: &[MotifSite],
    controls: &[BgInterval],
    motif_id: &[Vec<usize>],
    tf_names: &[String],
    n_motif: &[usize],
    n_controls: &[usize],
    fai_chrom_lengths: &[(String, u32)],
    window: usize,
    genome_name: &str,
    motif_set_name: &str,
) -> io::Result<()> {
    let n_tfs = tf_names.len();

    // binary format doesnt store chr names on each row, Site and control positions carry a u8 chromsomeome id
    // this chrom id table will be stored in the header
    /*** in R:
    chrom_lengths <- seqlengths(x = genome)
    chroms_used <- unique(x = c(as.character(x=seqnames(x= motifs)), as.character(x=seqnames(x=controls))))
    chrom_lengths <- chrom_lengths[intersect(names(x = chrom_lengths), chrom_used)]
    chrom_lengths <- chrom_lengths[!is.na(x=chrom_lengths)]
    */
    let mut chroms_used: HashSet<&str> = sites
        .iter()
        .map(|s| s.chrom.as_str())
        .chain(controls.iter().map(|c| c.chrom.as_str()))
        .collect();
    // subset fai_chrom_lengths to chrom_used
    let chrom_lengths: Vec<&(String, u32)> = fai_chrom_lengths
        .iter()
        .filter(|(n, _)| chroms_used.remove(n.as_str()))
        .collect();
    // raise concern if there is chromosomes from bed file not in the genome reference.
    if !chroms_used.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("Genome reference does not contain entry for Chromomsome {chroms_used:?}"),
        ));
    }
    let chrom_to_id: FxHashMap<&str, u8> = chrom_lengths
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (n.as_str(), i as u8))
        .collect();
    let f = File::create(path)?;
    let mut w = BufWriter::new(XzEncoder::new(f, 6));
    // magic + version
    w.write_all(b"BACH")?;
    w.write_all(&[1u8])?; // change the version here lol

    // header
    write_str(&mut w, genome_name)?;
    write_str(&mut w, motif_set_name)?;
    w.write_all(&(window as u32).to_le_bytes())?;
    w.write_all(&(n_tfs as u16).to_le_bytes())?;
    w.write_all(&[chrom_lengths.len() as u8])?;
    for name in tf_names {
        write_str(&mut w, name)?;
    }
    for &c in n_motif {
        w.write_all(&(c as u32).to_le_bytes())?;
    }
    for &c in n_controls {
        w.write_all(&(c as u32).to_le_bytes())?;
    }
    for (name, len) in chrom_lengths {
        write_str(&mut w, &name)?;
        w.write_all(&len.to_le_bytes())?;
    }

    // site positions, sorted by chrom id then start
    let mut ord: Vec<usize> = (0..sites.len()).collect();
    ord.sort_by_key(|&i| (chrom_to_id[sites[i].chrom.as_str()], sites[i].start));
    let s_chrom: Vec<u8> = ord
        .iter()
        .map(|&i| chrom_to_id[sites[i].chrom.as_str()])
        .collect();
    let s_start: Vec<u32> = ord.iter().map(|&i| sites[i].start as u32 + 1).collect();
    // + 1 in s_start cause R is 1-based so the reader in bacha.rs subtracts 1 from every start ....
    // + 1 also in the c_start
    let s_tf: Vec<u16> = ord.iter().map(|&i| sites[i].tf_index as u16).collect();
    let s_delta = delta_per_chrom(&s_chrom, &s_start);
    w.write_all(&(sites.len() as u32).to_le_bytes())?;
    w.write_all(&s_chrom)?;
    for v in &s_delta {
        w.write_all(&v.to_le_bytes())?;
    }
    for v in &s_tf {
        w.write_all(&v.to_le_bytes())?;
    }

    // control positions
    let mut cord: Vec<usize> = (0..controls.len()).collect();
    cord.sort_by_key(|&i| (chrom_to_id[controls[i].chrom.as_str()], controls[i].start));
    let c_chrom: Vec<u8> = cord
        .iter()
        .map(|&i| chrom_to_id[controls[i].chrom.as_str()])
        .collect();
    let c_start: Vec<u32> = cord.iter().map(|&i| controls[i].start as u32 + 1).collect();
    let c_delta = delta_per_chrom(&c_chrom, &c_start);
    let n_per_row: Vec<u16> = cord.iter().map(|&i| motif_id[i].len() as u16).collect();

    // flatten
    // rows must be ascending for the deltas to be decodable
    let mut flat: Vec<u16> = Vec::new();
    for &i in &cord {
        // revisit when parallizing the TF-for-loop
        // let mut ids: Vec<usize> = motif_id[i].clone();
        //ids.sort_unstable();
        let mut prev: Option<u16> = None;
        for &tf in &motif_id[i] {
            // for tf in ids {
            let abs = (tf + n_tfs) as u16;
            flat.push(match prev {
                Some(p) => abs - p,
                None => abs,
            });
            prev = Some(abs);
        }
    }
    w.write_all(&(controls.len() as u32).to_le_bytes())?;
    w.write_all(&c_chrom)?;
    for v in &c_delta {
        w.write_all(&v.to_le_bytes())?;
    }
    for v in &n_per_row {
        w.write_all(&v.to_le_bytes())?;
    }
    w.write_all(&(flat.len() as u32).to_le_bytes())?;
    for v in &flat {
        w.write_all(&v.to_le_bytes())?;
    }
    w.into_inner().map_err(|e| e.into_error())?.finish()?;
    Ok(())
}

// helper functions
fn generate_background_intervals(
    peak_file: &Path, // OCR --> background regions
    ref_path: &str,
    bias_h5: &str,
    window: usize, //50
) -> io::Result<Vec<BgInterval>> {
    /////////////
    let bias_factors = read_bias_factors(bias_h5);
    let mut faidx = load_fasta(ref_path).expect("Could not load ref fasta");
    let mut seq = Vec::new();
    //////////
    let mut bg_intervals: Vec<BgInterval> = Vec::new();
    let step = window / 2;
    let file = File::open(peak_file)?;
    let reader = BufReader::new(file);

    for (index, line) in reader.lines().enumerate() {
        let line = line.expect("Error reading line in Peak File");
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }
        let chr = fields[0].to_string();
        let start: usize = fields[1]
            .parse()
            .expect(&format!("Invalid start position {}", index + 1));
        let end: usize = fields[2]
            .parse()
            .expect(&format!("Invalid end position {}", index + 1));
        let peak_len = end - start;
        if peak_len == 0 {
            continue;
        }

        // extend OCR (resize(.., fix = "center"))
        let center = start + peak_len / 2;
        let mut ext_start = center.saturating_sub(peak_len / 2 + window);
        let mut ext_end = center + peak_len / 2 + window;

        // chromosome trimming (trim)
        let chr_trimmed = &chr[3..];
        //let chr_len = faidx.seqlen(chr_trimmed).expect("unknown chromosome") as usize;
        let chr_len = faidx
            .index
            .sequences()
            .iter()
            .find(|s| s.name == chr_trimmed)
            .map(|s| s.len as usize)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::Other,
                    format!("chromosome {} not found in FASTA index", chr_trimmed),
                )
            })?;
        ext_start = ext_start.min(chr_len);
        ext_end = ext_end.min(chr_len);

        // sliding window instead of random seed (slidingWindow)
        // new edit: fetch per peak instead of per window
        seq.clear();
        faidx
            .fetch(chr_trimmed, ext_start as u64, ext_end as u64)
            .expect("Couldn't fetch interval");
        faidx.read(&mut seq).expect("Couldn't read the interval");
        let mut pos = ext_start;
        while pos + window <= ext_end {
            let off = pos - ext_start;
            let gc_content = gc_percent(&seq[off..off + window]);
            let mean_bias = cal_mean_bias(&chr, pos, pos + window, &bias_factors);
            bg_intervals.push(BgInterval {
                chrom: chr.clone(),
                start: pos,
                stop: pos + window,
                features: [gc_content, mean_bias],
            });
            pos += step;
        }
    }
    Ok(bg_intervals)
}
// read_bias_factor & cal_mean_bias & gc_percent
pub fn read_bias_factors(h5_path: &str) -> FxHashMap<String, Vec<f32>> {
    let file = hdf5_file::open(h5_path).expect("Unable to open HDF5 file");
    let mut bias_factors: FxHashMap<String, Vec<f32>> = FxHashMap::default();
    // define range of chrom to use // throw the rest to save memory
    let chromosomes: Vec<String> = (1..=22).map(|i| format!("chr{}", i)).collect();
    for key in chromosomes {
        let dataset = file.dataset(&key).expect("Unable to open dataset");
        // Read the dataset into an ndarray Array
        let array: Array1<f32> = dataset.read().expect("Unable to read dataset");
        // Convert the ndarray Array into a Vec<f32>
        let bias: Vec<f32> = array.iter().copied().collect();
        // let bias: Vec<f32> = dataset.read().expect("Unable to read dataset");
        bias_factors.insert(key, bias);
    }
    bias_factors
}
pub fn cal_mean_bias(
    chrom: &String,
    start: usize,
    end: usize,
    bias_factors: &FxHashMap<String, Vec<f32>>,
) -> f32 {
    if let Some(bias) = bias_factors.get(chrom) {
        let bias_sum: f32 = bias[start..end].iter().copied().sum();
        let bias_mean = bias_sum / ((end - start) as f32);
        bias_mean
    } else {
        0.167
    }
}

// Compute GC% for a DNA sequence
pub fn gc_percent(seq: &[u8]) -> f32 {
    let gc = seq.iter().filter(|&&c| c == b'G' || c == b'C').count();
    let len = seq.len();
    if len == 0 {
        return 0.0;
    }
    gc as f32 / len as f32
}

pub fn load_fasta(fasta_path: &str) -> Result<IndexedReader<File>, std::io::Error> {
    let fasta_path = Path::new(fasta_path);
    let fai_path = fasta_path.with_file_name(format!(
        "{}.fai",
        fasta_path.file_name().unwrap().to_string_lossy()
    ));
    println!("fai_path:{}", fai_path.display());

    let fasta_file = File::open(fasta_path)?;
    println!("opened fasta file");
    let fai_file = File::open(fai_path)?;
    println!("opened fai file");
    let faidx = IndexedReader::new(fasta_file, fai_file)?;
    // if dont want result<>, just replace "?" with .unwrap()

    Ok(faidx)
    //https://docs.rs/bio/latest/bio/io/fasta/
}

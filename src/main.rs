#[cfg(not(target_os = "windows"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use clap::{Parser, Subcommand};
use std::error::Error;
use std::path::{Path, PathBuf};
mod bacha;
// The index builder is WIP
#[allow(dead_code)]
mod index;

/// Canonicalize a user-supplied input path, turning a missing or unreadable
/// path into a clean error (surfaced by `main` via `Display`) instead of a
/// panic with a backtrace.
fn canonicalize_input(path: &str, label: &str) -> Result<PathBuf, Box<dyn Error>> {
    Path::new(path)
        .canonicalize()
        .map_err(|e| format!("cannot access {label} '{path}': {e}").into())
}

#[derive(Parser, Debug)]
#[command(
    name = env!("CARGO_PKG_NAME"),
    author = env!("CARGO_PKG_AUTHORS"),
    about = env!("CARGO_PKG_DESCRIPTION"),
    version = env!("CARGO_PKG_VERSION"),
    long_about = None
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Count Tn5 insertions and compute activity scores
    Count {
        /// Path to the fragment file
        #[arg(
            short = 'f',
            long = "fragments",
            required_unless_present = "manifest",
            conflicts_with = "manifest",
            requires = "cells"
        )]
        fragments: Option<String>,

        /// Bacha index file
        #[arg(short = 'i', long = "index")]
        index: String,

        /// File containing cell barcodes to include
        #[arg(
            short = 'c',
            long = "cells",
            required_unless_present = "manifest",
            conflicts_with = "manifest",
            requires = "fragments",
            long_help = "File containing cell barcodes to include. One barcode per line, or \
                         two tab-separated columns where the first is the barcode as it \
                         appears in the fragment file and the second is the name to give \
                         that cell in the output."
        )]
        cells: Option<String>,

        /// Tab-separated manifest of several fragment files to score together
        #[arg(
            short = 'm',
            long = "manifest",
            conflicts_with_all = ["fragments", "cells"],
            long_help = "Tab-separated manifest listing several fragment files to score into \
                         one output matrix, as an alternative to -f/-c. Columns are sample, \
                         fragments and cells; a header row naming exactly those three is \
                         optional. Blank lines and lines beginning with '#' are ignored, and \
                         relative paths are resolved against the manifest's own directory. \
                         Files are counted one at a time, each using all --threads workers, \
                         and the Bacha index is loaded once for the whole run. Every cell \
                         must belong to exactly one fragment file: a Bacha score is a \
                         non-linear function of a cell's insertion counts, so counts for one \
                         cell split across two files cannot be combined."
        )]
        manifest: Option<String>,

        /// Output file name (CSV or .gz)
        #[arg(
            short = 'o',
            long = "output",
            long_help = "Output file name. Output will be a CSV file. If the file extension is .gz, the file will be compressed with gzip."
        )]
        output: String,

        /// Number of threads to use
        #[arg(
            short = 't',
            long = "threads",
            default_value_t = 6,
            long_help = "Total number of threads to use. Against a BGZF \
                         (bgzip-compressed) fragment file the budget is split: 2 threads \
                         inflate the file and the rest count it, so the two pools share \
                         --threads rather than adding to it. A plain gzip file is inflated \
                         inline by the reader thread, so every thread goes to counting. Also \
                         sets the number of gzip compression threads when the output path \
                         ends in .gz."
        )]
        threads: usize,

        /// Shrinkage strength parameter (lambda)
        #[arg(short = 'l', long = "lambda", default_value_t = 10.0)]
        lambda: f32,

        /// Pseudocount added to counts before log fold change calculation
        /// to avoid log(0).
        #[arg(short = 'p', long = "pseudocount", default_value_t = 1.0)]
        pseudocount: f32,

        /// Skip the genome/assembly bounds check that aborts when a fragment
        /// coordinate falls beyond the end of its chromosome in the index.
        #[arg(
            long = "no-bounds-check",
            long_help = "Disable the genome/assembly mismatch check. By default the run aborts \
                         on the first fragment coordinate that lies beyond the end of its \
                         chromosome in the index (an impossible position unless the fragment \
                         file and index were built on different genome assemblies). Pass this \
                         flag to skip the check and process such coordinates anyway."
        )]
        no_bounds_check: bool,
    },
    /// View index file contents in plain text
    View {
        /// Bacha index file to view
        index: String,

        /// Print header information only
        #[arg(short = 'H', long = "header-only")]
        header_only: bool,
    },
    /// Print summary information about an index file (header + per-TF site
    /// and control counts) without reading the site position blocks.
    Summary {
        /// Bacha index file to summarise
        index: String,
    },
}

fn main() {
    if let Err(e) = run() {
        // Print the Display form (not Debug) so users see a clean message
        // rather than `Custom { kind: InvalidInput, error: "..." }`.
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    // Reset SIGPIPE to default behavior so piping to `less`/`head` exits cleanly
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }

    pretty_env_logger::init_timed();

    let cli = Cli::parse();

    match cli.command {
        Commands::Count {
            fragments,
            index,
            cells,
            manifest,
            output,
            threads,
            lambda,
            pseudocount,
            no_bounds_check,
        } => {
            let index_file = canonicalize_input(&index, "index file")?;
            let outfile_path = Path::new(&output);

            match manifest {
                Some(manifest) => {
                    let manifest_file = canonicalize_input(&manifest, "manifest file")?;
                    bacha::bacha_manifest(
                        &manifest_file,
                        &index_file,
                        outfile_path,
                        threads,
                        lambda,
                        pseudocount,
                        !no_bounds_check,
                    )?;
                }
                None => {
                    // clap guarantees both are present when --manifest is absent.
                    let fragments = fragments.expect("clap requires -f without -m");
                    let cells = cells.expect("clap requires -c without -m");
                    let frag_file = canonicalize_input(&fragments, "fragment file")?;
                    let cell_file = canonicalize_input(&cells, "cell file")?;
                    bacha::bacha(
                        &frag_file,
                        &index_file,
                        &cell_file,
                        outfile_path,
                        threads,
                        lambda,
                        pseudocount,
                        !no_bounds_check,
                    )?;
                }
            }
        }
        Commands::View { index, header_only } => {
            let index_file = canonicalize_input(&index, "index file")?;
            bacha::view_index(&index_file, header_only)?;
        }
        Commands::Summary { index } => {
            let index_file = canonicalize_input(&index, "index file")?;
            bacha::view_index(&index_file, true)?;
        }
    }

    Ok(())
}

/// Keeps the CLI reference in README.md in sync with the clap definition, which
/// is the single source of truth for arguments. The reference lives between the
/// marker comments below; `regen_readme` renders it from `Cli` and either checks
/// it (default, so CI fails on drift) or rewrites it when `UPDATE_README` is set.
#[cfg(test)]
mod readme_sync {
    use super::Cli;

    const BEGIN: &str =
        "<!-- BEGIN CLI REFERENCE (generated by `UPDATE_README=1 cargo test regen_readme`) -->";
    const END: &str = "<!-- END CLI REFERENCE -->";

    fn readme_path() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("README.md")
    }

    /// Render the CLI reference block (markers included) from the clap command.
    fn expected_block() -> String {
        let options = clap_markdown::MarkdownOptions::new()
            .show_footer(false)
            .show_table_of_contents(false);
        let md = clap_markdown::help_markdown_custom::<Cli>(&options);
        format!("{BEGIN}\n\n{}\n\n{END}", fit_headings(&md))
    }

    /// Adapt clap-markdown's output to nest under the README's `### Command
    /// line` heading: drop its top-level document title and boilerplate intro
    /// line, then demote every remaining heading two levels (clamped at h6, the
    /// deepest markdown heading) so per-command sections become h4.
    fn fit_headings(md: &str) -> String {
        let mut out: Vec<String> = Vec::new();
        for line in md.lines() {
            let hashes = line.chars().take_while(|&c| c == '#').count();
            let is_heading = hashes > 0 && line[hashes..].starts_with(' ');
            if is_heading {
                // Drop clap-markdown's single `#` document title outright.
                if hashes == 1 {
                    continue;
                }
                let level = (hashes + 2).min(6);
                out.push(format!("{} {}", "#".repeat(level), &line[hashes + 1..]));
            } else if line.starts_with("This document contains the help content") {
                continue;
            } else {
                out.push(line.to_string());
            }
        }
        out.join("\n").trim().to_string()
    }

    /// Replace the marked region of `readme` with `block`.
    fn splice(readme: &str, block: &str) -> String {
        let start = readme
            .find(BEGIN)
            .expect("README.md is missing the BEGIN CLI REFERENCE marker");
        let end = readme[start..]
            .find(END)
            .map(|i| start + i + END.len())
            .expect("README.md is missing the END CLI REFERENCE marker");
        format!("{}{}{}", &readme[..start], block, &readme[end..])
    }

    #[test]
    fn regen_readme() {
        let path = readme_path();
        let current = std::fs::read_to_string(&path).expect("read README.md");
        let updated = splice(&current, &expected_block());

        if std::env::var_os("UPDATE_README").is_some() {
            if updated != current {
                std::fs::write(&path, &updated).expect("write README.md");
            }
        } else {
            assert!(
                updated == current,
                "README.md CLI reference is out of date with the clap definition. \
                 Run `UPDATE_README=1 cargo test regen_readme` to regenerate it."
            );
        }
    }
}

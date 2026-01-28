//! CLI entry point for rgmatch.
//!
//! This provides a command-line interface matching the Python implementation.

use anyhow::{bail, Context, Result};
use clap::Parser;
use std::path::PathBuf;

use rgmatch::config::Config;
use rgmatch::matcher::match_regions_to_genes;
use rgmatch::output::write_results;
use rgmatch::parser::{parse_bed, parse_gtf};
use rgmatch::types::ReportLevel;

/// Genomic region-to-gene matching tool.
///
/// Maps genomic regions from a BED file to gene annotations from a GTF file.
#[derive(Parser, Debug)]
#[command(name = "rgmatch")]
#[command(author, version, about, long_about = None)]
struct Args {
    /// GTF annotation file (required)
    #[arg(short = 'g', long = "gtf")]
    gtf: PathBuf,

    /// Region BED file (required)
    #[arg(short = 'b', long = "bed")]
    bed: PathBuf,

    /// Output file (required)
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    /// Report level: exon, transcript, or gene
    #[arg(short = 'r', long = "report", default_value = "exon")]
    report: String,

    /// Maximum distance in kb to report associations
    #[arg(short = 'q', long = "distance", default_value = "10")]
    distance: i64,

    /// TSS region distance in bp
    #[arg(short = 't', long = "tss", default_value = "200")]
    tss: i64,

    /// TTS region distance in bp
    #[arg(short = 's', long = "tts", default_value = "0")]
    tts: i64,

    /// Promoter region distance in bp
    #[arg(short = 'p', long = "promoter", default_value = "1300")]
    promoter: i64,

    /// Percentage of the area overlap threshold (0-100)
    #[arg(short = 'v', long = "perc_area", default_value = "90")]
    perc_area: f64,

    /// Percentage of the region overlap threshold (0-100)
    #[arg(short = 'w', long = "perc_region", default_value = "50")]
    perc_region: f64,

    /// Priority rules (comma-separated)
    #[arg(short = 'R', long = "rules", default_value = "TSS,1st_EXON,PROMOTER,TTS,INTRON,GENE_BODY,UPSTREAM,DOWNSTREAM")]
    rules: String,

    /// GTF tag for gene ID
    #[arg(short = 'G', long = "gene", default_value = "gene_id")]
    gene_tag: String,

    /// GTF tag for transcript ID
    #[arg(short = 'T', long = "transcript", default_value = "transcript_id")]
    transcript_tag: String,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Validate inputs
    if !args.gtf.exists() {
        bail!("GTF file not found: {}", args.gtf.display());
    }
    if !args.bed.exists() {
        bail!("BED file not found: {}", args.bed.display());
    }

    // Parse report level
    let level = ReportLevel::from_str(&args.report).context(
        "Report can only be one of the following: exon, transcript or gene",
    )?;

    // Build configuration
    let mut config = Config::new();
    config.level = level;

    // Set distance (convert from kb to bp)
    if args.distance >= 0 {
        config.set_distance_kb(args.distance);
    }

    // Set TSS distance
    if args.tss >= 0 {
        config.tss = args.tss as f64;
    } else {
        bail!("The TSS distance cannot be lower than 0 bps.");
    }

    // Set TTS distance
    if args.tts >= 0 {
        config.tts = args.tts as f64;
    } else {
        bail!("The TTS distance cannot be lower than 0 bps.");
    }

    // Set promoter distance
    if args.promoter >= 0 {
        config.promoter = args.promoter as f64;
    } else {
        bail!("The promoter distance cannot be lower than 0 bps.");
    }

    // Set percentage thresholds
    if args.perc_area >= 0.0 && args.perc_area <= 100.0 {
        config.perc_area = args.perc_area;
    } else {
        bail!("The percentage of area defined was wrong. It should range between 0 and 100.");
    }

    if args.perc_region >= 0.0 && args.perc_region <= 100.0 {
        config.perc_region = args.perc_region;
    } else {
        bail!("The percentage of region defined was wrong. It should range between 0 and 100.");
    }

    // Parse rules
    if !config.parse_rules(&args.rules) {
        bail!("Rules not properly passed.");
    }

    // Set GTF tags
    config.gene_id_tag = args.gene_tag;
    config.transcript_id_tag = args.transcript_tag;

    // Parse GTF file
    eprintln!("Parsing GTF file: {}", args.gtf.display());
    let gtf_data = parse_gtf(&args.gtf, &config.gene_id_tag, &config.transcript_id_tag)?;

    // Parse BED file
    eprintln!("Parsing BED file: {}", args.bed.display());
    let bed_data = parse_bed(&args.bed)?;

    // Process each chromosome
    eprintln!("Matching regions to genes...");
    let mut all_results = Vec::new();

    // Process each chromosome in sorted order for deterministic output
    let mut chroms: Vec<_> = bed_data.regions_by_chrom.keys().collect();
    chroms.sort();

    for chrom in chroms {
        let regions = &bed_data.regions_by_chrom[chrom];
        if let Some(mut genes) = gtf_data.genes_by_chrom.get(chrom).cloned() {
            // Sort regions by start position
            let mut sorted_regions = regions.clone();
            sorted_regions.sort_by_key(|r| r.start);

            let results = match_regions_to_genes(&sorted_regions, &mut genes, &config);
            all_results.extend(results);
        } else {
            eprintln!("Warning: {} not found in genes", chrom);
        }
    }

    // Write output
    eprintln!("Writing output to: {}", args.output.display());
    write_results(&args.output, &all_results, bed_data.num_meta_columns)?;

    eprintln!("Done!");
    Ok(())
}

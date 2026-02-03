# rgmatch-rs

[![Build Status](https://img.shields.io/github/actions/workflow/status/TianYuan-Liu/rgmatch-rs/ci.yml?branch=master)](https://github.com/TianYuan-Liu/rgmatch-rs/actions)
![Crate](https://img.shields.io/badge/crates.io-coming%20soon-yellow)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**A high-performance Rust implementation of the RGmatch tool for genomic interval matching.**

`rgmatch-rs` is a specialized bioinformatics tool designed to associate genomic regions (provided in BED format) with proximal gene features (from GTF annotation files). It provides flexible, rule-based annotation at the exon, transcript, or gene level, making it essential for integrating omics data such as ChIP-seq, ATAC-seq, or SMP data.

## Features

- **High Performance**: Optimized Rust implementation offers significant speedups over the original Python version.
- **Flexible Reporting**: Output associations at the **exon**, **transcript**, or **gene** level.
- **Detailed Annotations**: Identifies overlaps with exons, introns, promoters, TSS, TTS, and intergenic regions.
- **Customizable Rules**: Users can define priority rules for overlapping features (e.g., prioritize TSS over Exons).
- **Parallel Processing**: multi-threaded execution for handling large datasets efficiently.
- **Streaming Support**: Capable of processing large genomic files with constant memory usage.

## Credits

- **Original Author**: Pedro Furió-Tarí
- **Current Maintainer (Rust Version)**: Tianyuan Liu

## Citation

If you use `rgmatch-rs` in your research, please cite the original publication:

> Furió-Tarí P, Conesa A, Tarazona S. **RGmatch: matching genomic regions to proximal genes in omics data integration.** *BMC Bioinformatics*. 2016;17(Suppl 15):427.
>
> **DOI**: [10.1186/s12859-016-1293-1](https://doi.org/10.1186/s12859-016-1293-1) | **PMID**: [28185573](https://pubmed.ncbi.nlm.nih.gov/28185573/) | **PMCID**: PMC5133492

## Installation

### From Source

Ensure you have [Rust](https://www.rust-lang.org/tools/install) installed (version 1.70 or later).

```bash
# Clone the repository
git clone https://github.com/TianYuan-Liu/rgmatch-rs.git
cd rgmatch-rs

# Build in release mode
cargo build --release

# The binary will be located at:
./target/release/rgmatch
```

## Usage

### Basic Command

```bash
rgmatch -g annotations.gtf.gz -b regions.bed -o output.txt
```

### Options

| Support | Option | Description | Default |
|:-------:|:-------|:------------|:--------|
| **Input** | `-g`, `--gtf` | Path to GTF annotation file (supports .gz) | Required |
| **Input** | `-b`, `--bed` | Path to BED file with regions | Required |
| **Output** | `-o`, `--output` | Output file path | Required |
| **Mode** | `-r`, `--report` | Report level: `exon`, `transcript`, or `gene` | `exon` |
| **Parallel**| `-j`, `--threads` | Number of worker threads | `8` |
| **Config** | `-q`, `--distance`| Max distance (kb) for upstream/downstream | `10` |
| **Config** | `-t`, `--tss` | TSS region size (bp) | `200` |
| **Config** | `-s`, `--tts` | TTS region size (bp) | `0` |
| **Config** | `-p`, `--promoter`| Promoter region size (bp) | `1300` |
| **Filter** | `-v`, `--perc_area`| Min % of feature covered | `90` |
| **Filter** | `-w`, `--perc_region`| Min % of region covered | `50` |
| **Rules** | `-R`, `--rules` | Priority rules (comma-separated) | *See below* |

### Priority Rules

The `--rules` flag controls the priority when a region overlaps multiple features.
**Default Configuration:**
```text
TSS > 1st_EXON > GENE_BODY > PROMOTER > INTRON > TTS > UPSTREAM > DOWNSTREAM
```
You can customize this order, e.g., to prioritize Promoters over TSS:
`-R PROMOTER,TSS,1st_EXON,...`

### Output Format

The output is a tab-separated file containing the original BED fields followed by `rgmatch` annotations:

| Column | Description |
|:-------|:------------|
| `AREA` | Feature type (e.g., TSS, EXON, INTRON) |
| `GENE` | Gene ID |
| `TRANSCRIPT` | Transcript ID |
| `EXON_NR` | Exon number(s) |
| `STRAND` | Strand (`+` or `-`) |
| `DISTANCE` | Distance to feature (0 if overlapping) |
| `TSS_DISTANCE` | Distance to Transcription Start Site |
| `PCTG_DHS` | Percentage of the input region covered |
| `PCTG_AREA` | Percentage of the genomic feature covered |

## Testing

Run the comprehensive test suite to ensure correctness:

```bash
# Run all tests (library and integration)
cargo test
```

## Comparisons

`rgmatch-rs` is designed to be a drop-in high-performance replacement for the original Python implementation.

- **Speed**: significantly faster due to native compilation and parallelization.
- **Memory**: Optimized to handle large datasets with low memory footprint using streaming.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

//! BED file parser with gzip support.
//!
//! Parses BED (Browser Extensible Data) files containing genomic regions.

use ahash::AHashMap;
use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::types::Region;

/// Result of parsing a BED file.
pub struct BedData {
    /// Regions organized by chromosome.
    pub regions_by_chrom: AHashMap<String, Vec<Region>>,
    /// Number of metadata columns found.
    pub num_meta_columns: usize,
}

/// Parse a BED file and return organized region data.
///
/// Supports both plain text and gzip-compressed BED files.
pub fn parse_bed(path: &Path) -> Result<BedData> {
    let file = File::open(path).context("Failed to open BED file")?;

    let reader: Box<dyn BufRead> = if path
        .to_string_lossy()
        .ends_with(".gz")
    {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    parse_bed_reader(reader)
}

/// Parse BED data from a reader.
fn parse_bed_reader<R: BufRead>(reader: R) -> Result<BedData> {
    let mut regions_by_chrom: AHashMap<String, Vec<Region>> = AHashMap::new();
    let mut num_meta_columns = 0;

    for line_result in reader.lines() {
        let line = line_result.context("Failed to read BED line")?;

        // Skip empty lines
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        // Need at least 3 columns: chrom, start, end
        if fields.len() < 3 {
            continue;
        }

        let chrom = fields[0].to_string();

        // Try to parse start and end as integers
        // If they fail (e.g., header line), skip this line
        let start: i64 = match fields[1].parse() {
            Ok(v) => v,
            Err(_) => continue, // Skip header lines
        };
        let end: i64 = match fields[2].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        // Extract up to 9 additional BED columns as metadata
        let metadata: Vec<String> = fields
            .iter()
            .skip(3)
            .take(9)
            .map(|s| s.to_string())
            .collect();

        // Track the maximum number of metadata columns
        if metadata.len() > num_meta_columns {
            num_meta_columns = metadata.len();
        }

        let region = Region::new(chrom.clone(), start, end, metadata);
        regions_by_chrom.entry(chrom).or_default().push(region);
    }

    Ok(BedData {
        regions_by_chrom,
        num_meta_columns,
    })
}

/// Get standard BED column headers for metadata columns.
pub fn get_bed_headers(num_columns: usize) -> Vec<&'static str> {
    let all_headers = [
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
        "blockCount",
        "blockSizes",
        "blockStarts",
    ];

    all_headers.iter().take(num_columns).copied().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bed_basic() {
        let bed_content = "chr1\t100\t200\nchrom2\t300\t400\n";

        let reader = BufReader::new(bed_content.as_bytes());
        let result = parse_bed_reader(reader).unwrap();

        assert!(result.regions_by_chrom.contains_key("chr1"));
        assert!(result.regions_by_chrom.contains_key("chrom2"));

        let chr1_regions = &result.regions_by_chrom["chr1"];
        assert_eq!(chr1_regions.len(), 1);
        assert_eq!(chr1_regions[0].start, 100);
        assert_eq!(chr1_regions[0].end, 200);
        assert!(chr1_regions[0].metadata.is_empty());
    }

    #[test]
    fn test_parse_bed_with_metadata() {
        let bed_content = "chr1\t100\t200\tregion1\t500\t+\n";

        let reader = BufReader::new(bed_content.as_bytes());
        let result = parse_bed_reader(reader).unwrap();

        let regions = &result.regions_by_chrom["chr1"];
        assert_eq!(regions[0].metadata.len(), 3);
        assert_eq!(regions[0].metadata[0], "region1");
        assert_eq!(regions[0].metadata[1], "500");
        assert_eq!(regions[0].metadata[2], "+");
        assert_eq!(result.num_meta_columns, 3);
    }

    #[test]
    fn test_parse_bed_skip_header() {
        let bed_content = "chrom\tstart\tend\tname\nchr1\t100\t200\tregion1\n";

        let reader = BufReader::new(bed_content.as_bytes());
        let result = parse_bed_reader(reader).unwrap();

        // Should skip header line (can't parse 'start' as int)
        assert!(result.regions_by_chrom.contains_key("chr1"));
        assert!(!result.regions_by_chrom.contains_key("chrom"));
    }

    #[test]
    fn test_parse_bed_empty_lines() {
        let bed_content = "\nchr1\t100\t200\n\nchr1\t300\t400\n\n";

        let reader = BufReader::new(bed_content.as_bytes());
        let result = parse_bed_reader(reader).unwrap();

        let regions = &result.regions_by_chrom["chr1"];
        assert_eq!(regions.len(), 2);
    }

    #[test]
    fn test_get_bed_headers() {
        assert_eq!(get_bed_headers(0), Vec::<&str>::new());
        assert_eq!(get_bed_headers(3), vec!["name", "score", "strand"]);
        assert_eq!(
            get_bed_headers(9),
            vec![
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts"
            ]
        );
    }

    #[test]
    fn test_region_methods() {
        let region = Region::new("chr1".to_string(), 100, 200, vec!["test".to_string()]);

        assert_eq!(region.length(), 101);
        assert_eq!(region.midpoint(), 150);
        assert_eq!(region.id(), "chr1_100_200");
    }

    #[test]
    fn test_region_midpoint_integer_division() {
        // Test that midpoint uses integer division
        let region = Region::new("chr1".to_string(), 100, 201, vec![]);
        // (100 + 201) / 2 = 150 (integer division)
        assert_eq!(region.midpoint(), 150);
    }
}

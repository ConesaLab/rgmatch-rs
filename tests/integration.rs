//! Integration tests comparing Rust output against Python golden output.
//!
//! These tests verify that the Rust implementation produces identical results
//! to the Python reference implementation.

use std::collections::HashSet;
use std::fs;
use std::process::Command;

/// Helper to get path to test data directory
fn test_data_dir() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("test_data")
}

/// Helper to check if cargo is available and the binary is built
fn ensure_binary_built() -> Option<std::path::PathBuf> {
    let manifest_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let binary_path = manifest_dir.join("target").join("release").join("rgmatch");

    // Try debug build first
    let debug_path = manifest_dir.join("target").join("debug").join("rgmatch");
    if debug_path.exists() {
        return Some(debug_path);
    }
    if binary_path.exists() {
        return Some(binary_path);
    }
    None
}

/// Parse output file into a set of normalized lines for comparison
fn parse_output_lines(content: &str) -> HashSet<String> {
    content
        .lines()
        .skip(1) // Skip header
        .filter(|line| !line.is_empty())
        .map(|line| normalize_line(line))
        .collect()
}

/// Normalize a line for comparison (handle floating point formatting)
fn normalize_line(line: &str) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 10 {
        return line.to_string();
    }

    // Fields that need float normalization (indices 8 and 9 are PercRegion and PercArea)
    let mut normalized = Vec::new();
    for (i, field) in fields.iter().enumerate() {
        if i == 8 || i == 9 {
            // Normalize float to 2 decimal places
            if let Ok(f) = field.parse::<f64>() {
                normalized.push(format!("{:.2}", f));
            } else {
                normalized.push(field.to_string());
            }
        } else {
            normalized.push(field.to_string());
        }
    }
    normalized.join("\t")
}

#[test]
#[ignore] // Run with --ignored flag since it requires building the binary
fn test_golden_output_exon_level() {
    let binary = match ensure_binary_built() {
        Some(b) => b,
        None => {
            eprintln!("Binary not built, skipping integration test");
            return;
        }
    };

    let test_data = test_data_dir();
    let gtf = test_data.join("benchmark.gtf");
    let bed = test_data.join("benchmark.bed");
    let golden = test_data.join("golden_output.txt");

    if !gtf.exists() || !bed.exists() || !golden.exists() {
        eprintln!("Test data files not found, skipping integration test");
        return;
    }

    // Create temp output file
    let temp_output = std::env::temp_dir().join("rgmatch_test_output.txt");

    // Run the Rust binary
    let output = Command::new(&binary)
        .args([
            "-g",
            gtf.to_str().unwrap(),
            "-b",
            bed.to_str().unwrap(),
            "-o",
            temp_output.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to execute rgmatch");

    if !output.status.success() {
        eprintln!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("rgmatch failed to run");
    }

    // Compare outputs
    let rust_output = fs::read_to_string(&temp_output).expect("Failed to read output");
    let golden_output = fs::read_to_string(&golden).expect("Failed to read golden output");

    let rust_lines = parse_output_lines(&rust_output);
    let golden_lines = parse_output_lines(&golden_output);

    // Check for differences
    let only_in_rust: Vec<_> = rust_lines.difference(&golden_lines).collect();
    let only_in_golden: Vec<_> = golden_lines.difference(&rust_lines).collect();

    if !only_in_rust.is_empty() || !only_in_golden.is_empty() {
        eprintln!("Lines only in Rust output ({}):", only_in_rust.len());
        for line in only_in_rust.iter().take(5) {
            eprintln!("  {}", line);
        }
        eprintln!("Lines only in golden output ({}):", only_in_golden.len());
        for line in only_in_golden.iter().take(5) {
            eprintln!("  {}", line);
        }
        panic!(
            "Output mismatch: {} extra lines in Rust, {} missing from golden",
            only_in_rust.len(),
            only_in_golden.len()
        );
    }

    // Cleanup
    let _ = fs::remove_file(temp_output);
}

#[test]
#[ignore]
fn test_golden_output_transcript_level() {
    let binary = match ensure_binary_built() {
        Some(b) => b,
        None => {
            eprintln!("Binary not built, skipping integration test");
            return;
        }
    };

    let test_data = test_data_dir();
    let gtf = test_data.join("benchmark.gtf");
    let bed = test_data.join("benchmark.bed");
    let golden = test_data.join("golden_output_transcript.txt");

    if !gtf.exists() || !bed.exists() || !golden.exists() {
        eprintln!("Test data files not found, skipping integration test");
        return;
    }

    let temp_output = std::env::temp_dir().join("rgmatch_test_transcript_output.txt");

    let output = Command::new(&binary)
        .args([
            "-g",
            gtf.to_str().unwrap(),
            "-b",
            bed.to_str().unwrap(),
            "-o",
            temp_output.to_str().unwrap(),
            "-r",
            "transcript",
        ])
        .output()
        .expect("Failed to execute rgmatch");

    if !output.status.success() {
        eprintln!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("rgmatch failed to run");
    }

    let rust_output = fs::read_to_string(&temp_output).expect("Failed to read output");
    let golden_output = fs::read_to_string(&golden).expect("Failed to read golden output");

    let rust_lines = parse_output_lines(&rust_output);
    let golden_lines = parse_output_lines(&golden_output);

    let only_in_rust: Vec<_> = rust_lines.difference(&golden_lines).collect();
    let only_in_golden: Vec<_> = golden_lines.difference(&rust_lines).collect();

    if !only_in_rust.is_empty() || !only_in_golden.is_empty() {
        eprintln!(
            "Transcript level mismatch: {} extra, {} missing",
            only_in_rust.len(),
            only_in_golden.len()
        );
        panic!("Output mismatch at transcript level");
    }

    let _ = fs::remove_file(temp_output);
}

#[test]
#[ignore]
fn test_golden_output_gene_level() {
    let binary = match ensure_binary_built() {
        Some(b) => b,
        None => {
            eprintln!("Binary not built, skipping integration test");
            return;
        }
    };

    let test_data = test_data_dir();
    let gtf = test_data.join("benchmark.gtf");
    let bed = test_data.join("benchmark.bed");
    let golden = test_data.join("golden_output_gene.txt");

    if !gtf.exists() || !bed.exists() || !golden.exists() {
        eprintln!("Test data files not found, skipping integration test");
        return;
    }

    let temp_output = std::env::temp_dir().join("rgmatch_test_gene_output.txt");

    let output = Command::new(&binary)
        .args([
            "-g",
            gtf.to_str().unwrap(),
            "-b",
            bed.to_str().unwrap(),
            "-o",
            temp_output.to_str().unwrap(),
            "-r",
            "gene",
        ])
        .output()
        .expect("Failed to execute rgmatch");

    if !output.status.success() {
        eprintln!("stderr: {}", String::from_utf8_lossy(&output.stderr));
        panic!("rgmatch failed to run");
    }

    let rust_output = fs::read_to_string(&temp_output).expect("Failed to read output");
    let golden_output = fs::read_to_string(&golden).expect("Failed to read golden output");

    let rust_lines = parse_output_lines(&rust_output);
    let golden_lines = parse_output_lines(&golden_output);

    let only_in_rust: Vec<_> = rust_lines.difference(&golden_lines).collect();
    let only_in_golden: Vec<_> = golden_lines.difference(&rust_lines).collect();

    if !only_in_rust.is_empty() || !only_in_golden.is_empty() {
        eprintln!(
            "Gene level mismatch: {} extra, {} missing",
            only_in_rust.len(),
            only_in_golden.len()
        );
        eprintln!("Lines only in Rust output (Top 5):");
        for line in only_in_rust.iter().take(5) {
            eprintln!("  {}", line);
        }
        eprintln!("Lines only in golden output (Top 5):");
        for line in only_in_golden.iter().take(5) {
            eprintln!("  {}", line);
        }
        panic!("Output mismatch at gene level");
    }

    let _ = fs::remove_file(temp_output);
}

// Unit tests that don't require the binary
#[test]
fn test_normalize_line() {
    let line = "chr1_100_200\t150\tG1\tT1\t1\tTSS\t0\t100\t99.999\t50.001\tname";
    let normalized = normalize_line(line);
    assert!(normalized.contains("100.00"));
    assert!(normalized.contains("50.00"));
}

#[test]
fn test_parse_output_lines() {
    let content = "Header\nline1\nline2\n";
    let lines = parse_output_lines(content);
    assert_eq!(lines.len(), 2);
    assert!(!lines.contains("Header"));
}

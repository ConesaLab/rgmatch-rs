use assert_cmd::Command;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Helper function to run rgmatch and compare output against a golden file.
fn run_golden_test(
    resolution: &str,
    golden_filename: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let cargo_manifest_dir = env!("CARGO_MANIFEST_DIR");
    let base_dir = Path::new(cargo_manifest_dir);
    let data_dir = base_dir.join("tests").join("data");

    let gtf_path = data_dir.join("subset_genome.gtf");
    let bed_path = data_dir.join("subset_peaks.bed");
    let golden_path = data_dir.join(golden_filename);

    // Output file in data dir for inspection if it fails
    let output_path = data_dir.join(format!("integration_test_output_{}.txt", resolution));

    // Ensure test data exists
    if !gtf_path.exists() || !bed_path.exists() {
        panic!(
            "Test data not found. GTF: {:?}, BED: {:?}",
            gtf_path, bed_path
        );
    }

    if !golden_path.exists() {
        panic!("Golden output not found: {:?}", golden_path);
    }

    println!("Running integration test for resolution: {}", resolution);
    println!("GTF: {:?}", gtf_path);
    println!("BED: {:?}", bed_path);
    println!("Golden: {:?}", golden_path);

    // Run the binary
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rgmatch"));
    cmd.arg("-g")
        .arg(&gtf_path)
        .arg("-b")
        .arg(&bed_path)
        .arg("-r")
        .arg(resolution)
        .arg("-o")
        .arg(&output_path)
        .assert()
        .success();

    println!("Binary finished. Comparing outputs...");

    // Open files for streaming comparison
    let output_file = File::open(&output_path)?;
    let golden_file = File::open(&golden_path)?;

    let output_reader = BufReader::new(output_file);
    let golden_reader = BufReader::new(golden_file);

    // Compare line by line
    let mut output_lines = output_reader.lines();
    let mut golden_lines = golden_reader.lines();
    let mut line_num = 0;

    loop {
        line_num += 1;
        match (output_lines.next(), golden_lines.next()) {
            (Some(Ok(out_line)), Some(Ok(gold_line))) => {
                if out_line != gold_line {
                    panic!(
                        "Mismatch at line {}: \nExpected: {}\nActual:   {}",
                        line_num, gold_line, out_line
                    );
                }
            }
            (None, None) => break, // Both files ended
            (Some(_), None) => {
                panic!(
                    "Output file has more lines than golden file (golden ended at line {})",
                    line_num
                );
            }
            (None, Some(_)) => {
                panic!(
                    "Golden file has more lines than output file (output ended at line {})",
                    line_num
                );
            }
            (Some(Err(e)), _) => {
                panic!("Error reading output file at line {}: {}", line_num, e);
            }
            (_, Some(Err(e))) => {
                panic!("Error reading golden file at line {}: {}", line_num, e);
            }
        }
    }

    println!("All {} lines matched!", line_num - 1);

    // Clean up output file on success
    let _ = std::fs::remove_file(output_path);

    Ok(())
}

#[test]
fn test_golden_output_exon() -> Result<(), Box<dyn std::error::Error>> {
    run_golden_test("exon", "subset_golden_output_exon.txt")
}

#[test]
fn test_golden_output_transcript() -> Result<(), Box<dyn std::error::Error>> {
    run_golden_test("transcript", "subset_golden_output_transcript.txt")
}

#[test]
fn test_golden_output_gene() -> Result<(), Box<dyn std::error::Error>> {
    run_golden_test("gene", "subset_golden_output_gene.txt")
}

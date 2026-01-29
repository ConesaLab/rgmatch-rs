# Unit Test Coverage Task Notes

## Current Status
- Added 58 unit tests to `tests/unit_tests.rs` covering:
  - types.rs (Strand, Area, Exon, Transcript, Gene, Candidate, Region, ReportLevel)
  - config.rs (Config, parse_rules, set_distance_kb, max_lookback_distance)
  - matcher/tss.rs (check_tss with various strand/boundary conditions)
  - matcher/tts.rs (check_tts with various strand/boundary conditions)
  - matcher/rules.rs (apply_rules, select_transcript)
  - matcher/overlap.rs (find_search_start_index, process_candidates_for_output, match_regions_to_genes)

## What's Been Tested
All public functions in the crate now have unit test coverage either:
1. In-module tests (inside each source file's `#[cfg(test)] mod tests`)
2. Integration tests in `tests/unit_tests.rs`

## Next Steps for 100% Coverage
To achieve true 100% line/branch coverage, consider:

1. **Edge cases in overlap.rs**: The `match_region_to_genes` function has many branches (Cases 1-6). Add tests for:
   - Case 4: Exon overlapping region shifted right
   - Case 5: Region completely within exon
   - Negative strand edge cases for all overlap scenarios

2. **Parser edge cases**:
   - GTF files with missing attributes
   - BED files with malformed lines
   - Gzip file handling (requires actual .gz test files)

3. **Error path testing**:
   - Test error messages from ParseStrandError, ParseAreaError, ParseReportLevelError

4. **Run coverage tool**: Use `cargo tarpaulin` or `cargo llvm-cov` to get actual coverage percentages and identify uncovered lines.

## Commands
```bash
# Run all tests
cargo test

# Run only unit_tests.rs
cargo test --test unit_tests

# Get coverage (install tarpaulin first: cargo install cargo-tarpaulin)
cargo tarpaulin --out Html
```

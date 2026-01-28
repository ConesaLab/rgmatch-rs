//! TTS (Transcription Termination Site) overlap checking.
//!
//! This module implements the checkTTS logic with coordinate mirroring
//! for positive strand genes (opposite of TSS!).

use crate::types::Strand;

/// Result of a TTS check: (area_tag, pctg_dhs, pctg_area).
pub type TtsResult = (String, f64, f64);

/// Helper struct to pass exon-like data to checkTTS.
pub struct TtsExonInfo {
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub distance: i64,
}

/// Check overlap with TTS (Transcription Termination Site) region.
///
/// Calculates the overlap between a DHS region and the TTS/downstream
/// regions after the last exon. Handles strand orientation
/// by coordinate transformation.
///
/// CRITICAL: For POSITIVE strand (opposite of TSS!), coordinates are mirrored!
///
/// # Arguments
/// * `dhs_start` - Start coordinate of the DHS region
/// * `dhs_end` - End coordinate of the DHS region
/// * `exon_info` - Exon information including position, strand, and distance
/// * `tts_distance` - TTS region distance (default 0bp)
///
/// # Returns
/// A vector of (area_tag, pctg_dhs, pctg_area) tuples for each overlapping region type.
pub fn check_tts(
    dhs_start: i64,
    dhs_end: i64,
    exon_info: &TtsExonInfo,
    tts_distance: f64,
) -> Vec<TtsResult> {
    let mut exon_start = exon_info.start;
    let distance_val = exon_info.distance;
    let mut actual_dhs_start = dhs_start;
    let mut actual_dhs_end = dhs_end;

    // CRITICAL: Coordinate mirroring for POSITIVE strand (opposite of TSS!)
    // For positive strand, we flip the coordinates to make the code strand-invariant
    if exon_info.strand == Strand::Positive {
        let aux = actual_dhs_end;
        actual_dhs_end = 2 * exon_info.end - actual_dhs_start;
        actual_dhs_start = 2 * exon_info.end - aux;
        exon_start = exon_info.end; // TTS is at exon END for positive strand
    }

    let dhs_length = actual_dhs_end - actual_dhs_start + 1;

    // Zero-length region check - must be <= 0, not < 0
    if dhs_length <= 0 {
        return vec![];
    }

    let mut results = Vec::new();
    let dhs_length_f = dhs_length as f64;

    if distance_val as f64 <= tts_distance {
        // Region is within TTS distance

        // DOWNSTREAM       TTS        last exon
        // ..........|...............|----------->

        if (exon_start - actual_dhs_start) as f64 <= tts_distance {
            // Region is entirely within TTS zone
            // DOWNSTREAM        TTS          last exon
            // ..........|................|----------->
            //                      DHS
            //                  |-------------

            let overlap_end = std::cmp::min(exon_start - 1, actual_dhs_end);
            let overlap = overlap_end - actual_dhs_start + 1;
            let pctg_dhs = (overlap as f64 / dhs_length_f) * 100.0;
            let pctg_tts = (overlap as f64 / tts_distance) * 100.0;
            results.push(("TTS".to_string(), pctg_dhs, pctg_tts));
        } else {
            // Region spans TTS and extends into DOWNSTREAM
            // DOWNSTREAM         TTS          last exon
            // ............|..............|----------->
            //               DHS
            //       --------------

            // TTS portion
            let tts_start = exon_start - tts_distance as i64;
            let overlap_end = std::cmp::min(exon_start - 1, actual_dhs_end);
            let tts_overlap = overlap_end - tts_start + 1;
            let pctg_dhs_tts = (tts_overlap as f64 / dhs_length_f) * 100.0;
            let pctg_tts = (tts_overlap as f64 / tts_distance) * 100.0;
            results.push(("TTS".to_string(), pctg_dhs_tts, pctg_tts));

            // DOWNSTREAM portion
            let downstream_overlap = tts_start - actual_dhs_start;
            let pctg_dhs_downstream = (downstream_overlap as f64 / dhs_length_f) * 100.0;
            results.push(("DOWNSTREAM".to_string(), pctg_dhs_downstream, -1.0));
        }
    } else {
        // Region is entirely in DOWNSTREAM zone
        results.push(("DOWNSTREAM".to_string(), 100.0, -1.0));
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pos_strand_tts() {
        // Exon [1000, 2000]. TTS @ 2000.
        // Downstream > 2000.

        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };

        // Case 1: Downstream 100bp [2100, 2150]
        // Logic flip:
        // dhs_start' = 2*2000 - 2150 = 1850
        // dhs_end' = 2*2000 - 2100 = 1900
        // exon_start' = 2000 (from exon.getEnd())
        // Check: exon_start' - dhs_start' = 2000 - 1850 = 150 <= 200
        // Returns TTS
        let res = check_tts(2100, 2150, &exon, 200.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TTS"));
    }

    #[test]
    fn test_neg_strand_tts() {
        // Exon [1000, 2000]. Strand "-".
        // TTS @ 1000 (Start). Downstream < 1000.

        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Negative,
            distance: 0,
        };

        // Case 1: Downstream 100bp [850, 900]
        // No flip in code for "-".
        // exon_start = 1000
        // dhs_start = 850
        // exon_start - dhs_start = 1000 - 850 = 150 <= 200
        // Returns TTS
        let res = check_tts(850, 900, &exon, 200.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TTS"));
    }

    #[test]
    fn test_zero_tts_value() {
        // Test checkTTS when tts is set to 0
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 100,
        };
        let res = check_tts(2100, 2150, &exon, 0.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"DOWNSTREAM"));
    }

    #[test]
    fn test_zero_length_region_check_tts() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };
        let res = check_tts(2100, 2099, &exon, 200.0);
        assert!(res.is_empty());
    }

    #[test]
    fn test_large_tts_value() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 3000,
        };
        let res = check_tts(5000, 5100, &exon, 5000.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TTS"));
    }

    #[test]
    fn test_region_spanning_tts_downstream() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };
        // Region from 2100 to 2300 (spans TTS at 200bp)
        let res = check_tts(2100, 2300, &exon, 200.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"TTS"));
        assert!(tags.contains(&"DOWNSTREAM"));
    }

    #[test]
    fn test_neg_strand_percentage_calculations() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Negative,
            distance: 0,
        };
        let res = check_tts(850, 950, &exon, 200.0);
        for (tag, pctg_dhs, pctg_tts) in &res {
            if tag == "TTS" {
                assert!(*pctg_dhs >= 0.0 && *pctg_dhs <= 100.0);
                assert!(*pctg_tts >= 0.0 && *pctg_tts <= 100.0);
            }
        }
    }
}

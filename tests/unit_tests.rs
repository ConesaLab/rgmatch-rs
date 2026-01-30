//! Unit tests ported from Python test_unit.py
//!
//! These tests verify the core logic of rgmatch, especially coordinate mirroring
//! and priority rule application.

use rgmatch::config::Config;
use rgmatch::matcher::overlap::{
    find_search_start_index, match_region_to_genes, match_regions_to_genes,
    process_candidates_for_output,
};
use rgmatch::matcher::rules::{apply_rules, select_transcript};
use rgmatch::matcher::tss::{check_tss, TssExonInfo};
use rgmatch::matcher::tts::{check_tts, TtsExonInfo};
use rgmatch::types::{Area, Candidate, ReportLevel, Strand, Transcript};

// -------------------------------------------------------------------------
// Helper functions
// -------------------------------------------------------------------------

fn make_candidate(
    area: Area,
    pctg_region: f64,
    pctg_area: f64,
    transcript: &str,
    gene: &str,
    exon_number: &str,
) -> Candidate {
    Candidate::new(
        100,
        200,
        Strand::Positive,
        exon_number.to_string(),
        area,
        transcript.to_string(),
        gene.to_string(),
        0,
        pctg_region,
        pctg_area,
        100,
    )
}

fn default_rules() -> Vec<Area> {
    vec![
        Area::Tss,
        Area::FirstExon,
        Area::Promoter,
        Area::Tts,
        Area::Intron,
        Area::GeneBody,
        Area::Upstream,
        Area::Downstream,
    ]
}

// -------------------------------------------------------------------------
// 1. Coordinate Arithmetic (Fragile Math) Tests
// -------------------------------------------------------------------------

mod test_check_tss {
    use super::*;

    #[test]
    fn test_pos_strand_boundaries() {
        // Exon: [2000, 3000]. TSS @ 2000.
        // TSS zone: [1800, 2000].
        // Promoter zone: [500, 1800) (1300bp long).

        // Case 1: Exactly at TSS boundary (200bp upstream) -> [1800, 1810]
        // 2000 - 1800 = 200. <= 200 TSS.
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        let res = check_tss(1800, 1810, &exon, 200.0, 1300.0);
        assert!(
            res.iter().any(|(tag, _, _)| tag == "TSS"),
            "1800 should be TSS: {:?}",
            res
        );

        // Case 2: Just outside TSS boundary -> [1799, 1810]
        // 2000 - 1799 = 201. > 200. Should be PROMOTER.
        let res = check_tss(1799, 1810, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));
        assert!(tags.contains(&"TSS"));

        // Case 3: Far upstream (UPSTREAM tag)
        let exon_far = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 1800,
        };
        let res = check_tss(100, 200, &exon_far, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"UPSTREAM"));
        assert!(!tags.contains(&"TSS"));
        assert!(!tags.contains(&"PROMOTER"));
    }

    #[test]
    fn test_neg_strand_mirror() {
        // Exon: [2000, 3000]. Strand "-".
        // TSS @ 3000. Upstream > 3000.

        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Negative,
            distance: 0,
        };

        // Case 1: Region [3200, 3210] should be PROMOTER
        let res = check_tss(3200, 3210, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));

        // Case 2: TSS Zone Inside [3100, 3150]
        let res = check_tss(3100, 3150, &exon, 200.0, 1300.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TSS"));
    }

    #[test]
    fn test_integer_division_check_tss() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        let res = check_tss(1801, 1810, &exon, 200.0, 1300.0);
        assert!(!res.is_empty());
    }

    #[test]
    fn test_zero_length_region_check_tss() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        let res = check_tss(1900, 1899, &exon, 200.0, 1300.0);
        assert!(res.is_empty());
    }

    #[test]
    fn test_region_entirely_in_promoter() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 500, // Distance > tss but < tss + promoter
        };
        // Region in promoter zone (between TSS and UPSTREAM)
        let res = check_tss(1400, 1500, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));
        assert!(!tags.contains(&"TSS"));
    }

    #[test]
    fn test_region_spanning_promoter_upstream() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 500,
        };
        // Region spanning promoter and upstream zones
        let res = check_tss(100, 600, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));
        assert!(tags.contains(&"UPSTREAM"));
    }

    #[test]
    fn test_region_entirely_upstream() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 2000, // Far upstream
        };
        let res = check_tss(100, 200, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"UPSTREAM"));
        assert!(!tags.contains(&"TSS"));
        assert!(!tags.contains(&"PROMOTER"));
    }

    #[test]
    fn test_percentage_calculations_tss() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        // Region [1900, 1999] is 100bp long, entirely in 200bp TSS zone
        let res = check_tss(1900, 1999, &exon, 200.0, 1300.0);
        assert_eq!(res.len(), 1);
        let (tag, pctg_dhs, pctg_area) = &res[0];
        assert_eq!(tag, "TSS");
        assert_eq!(*pctg_dhs, 100.0); // 100bp / 100bp = 100%
        assert_eq!(*pctg_area, 50.0); // 100bp / 200bp = 50%
    }

    #[test]
    fn test_negative_strand_entirely_in_tss() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Negative,
            distance: 0,
        };
        // For negative strand, TSS is at exon end (3000)
        // TSS zone is [3001, 3200]
        let res = check_tss(3050, 3100, &exon, 200.0, 1300.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TSS"));
    }

    #[test]
    fn test_negative_strand_entirely_upstream() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Negative,
            distance: 2000,
        };
        let res = check_tss(5000, 5100, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"UPSTREAM"));
    }

    #[test]
    fn test_region_spanning_tss_promoter_upstream() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 0,
        };
        // Large region spanning all zones: TSS (200bp), Promoter (1300bp), and into Upstream
        // TSS zone: [1800, 2000]
        // Promoter zone: [500, 1800]
        // Upstream: < 500
        let res = check_tss(100, 1950, &exon, 200.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"TSS"));
        assert!(tags.contains(&"PROMOTER"));
        assert!(tags.contains(&"UPSTREAM"));
    }

    #[test]
    fn test_large_tss_distance() {
        let exon = TssExonInfo {
            start: 20000,
            end: 30000,
            strand: Strand::Positive,
            distance: 5000,
        };
        let res = check_tss(15000, 15100, &exon, 10000.0, 1300.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TSS"));
    }

    #[test]
    fn test_zero_tss_distance() {
        let exon = TssExonInfo {
            start: 2000,
            end: 3000,
            strand: Strand::Positive,
            distance: 500,
        };
        // With tss=0, should go directly to promoter logic
        let res = check_tss(1500, 1600, &exon, 0.0, 1300.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"PROMOTER"));
    }
}

mod test_check_tts {
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
        let res = check_tts(850, 900, &exon, 200.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TTS"));
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
    fn test_region_entirely_in_tts() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };
        // Region [2050, 2100] is entirely within TTS zone (200bp)
        let res = check_tts(2050, 2100, &exon, 200.0);
        assert_eq!(res.len(), 1);
        let (tag, _, _) = &res[0];
        assert_eq!(tag, "TTS");
    }

    #[test]
    fn test_region_entirely_downstream() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 500, // Beyond TTS distance
        };
        // Region far downstream
        let res = check_tts(2500, 2600, &exon, 200.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"DOWNSTREAM"));
        assert!(!tags.contains(&"TTS"));
    }

    #[test]
    fn test_region_spanning_tts_downstream() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };
        // Region from 2100 to 2300 spans TTS and DOWNSTREAM
        let res = check_tts(2100, 2300, &exon, 200.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"TTS"));
        assert!(tags.contains(&"DOWNSTREAM"));
    }

    #[test]
    fn test_percentage_calculations_tts() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 0,
        };
        // Region [2001, 2100] is 100bp long, entirely in 200bp TTS zone
        let res = check_tts(2001, 2100, &exon, 200.0);
        assert_eq!(res.len(), 1);
        let (tag, pctg_dhs, pctg_tts) = &res[0];
        assert_eq!(tag, "TTS");
        assert_eq!(*pctg_dhs, 100.0); // 100bp / 100bp = 100%
        assert_eq!(*pctg_tts, 50.0); // 100bp / 200bp = 50%
    }

    #[test]
    fn test_negative_strand_entirely_in_tts() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Negative,
            distance: 0,
        };
        // For negative strand, TTS is at exon start (1000)
        // TTS zone is [800, 999]
        let res = check_tts(850, 950, &exon, 200.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TTS"));
    }

    #[test]
    fn test_negative_strand_entirely_downstream() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Negative,
            distance: 500,
        };
        // For negative strand, downstream is < 1000
        let res = check_tts(400, 500, &exon, 200.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"DOWNSTREAM"));
    }

    #[test]
    fn test_zero_tts_distance() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 100,
        };
        // With tts=0, should go directly to downstream logic
        let res = check_tts(2100, 2200, &exon, 0.0);
        let tags: Vec<&str> = res.iter().map(|(tag, _, _)| tag.as_str()).collect();
        assert!(tags.contains(&"DOWNSTREAM"));
    }

    #[test]
    fn test_large_tts_distance() {
        let exon = TtsExonInfo {
            start: 1000,
            end: 2000,
            strand: Strand::Positive,
            distance: 3000,
        };
        // Region far from exon but within large TTS zone
        let res = check_tts(5000, 5100, &exon, 5000.0);
        assert!(res.iter().any(|(tag, _, _)| tag == "TTS"));
    }

    #[test]
    fn test_negative_strand_percentage_calculations() {
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

// -------------------------------------------------------------------------
// 2. Filtering Logic Tests
// -------------------------------------------------------------------------

mod test_apply_rules {
    use super::*;
    use ahash::AHashMap;

    #[test]
    fn test_priority_logic() {
        let rules = default_rules();

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1"); // Should win
        let c3 = make_candidate(Area::GeneBody, 100.0, 100.0, "T1", "G1", "1");

        let candidates = vec![c1, c2, c3];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("trans1".to_string(), vec![0, 1, 2]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_priority_logic_custom_rules() {
        // Change rules order - Intron now higher priority
        let rules = vec![Area::Intron, Area::Tss];

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("trans1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Intron);
    }

    #[test]
    fn test_pctg_region_threshold() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Intron, 60.0, 100.0, "T1", "G1", "1"); // Passes
        let c2 = make_candidate(Area::Tss, 40.0, 100.0, "T1", "G1", "1"); // Fails threshold

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Intron);
    }

    #[test]
    fn test_all_below_region_threshold() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Intron, 30.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 40.0, 100.0, "T1", "G1", "1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 90.0, 90.0, &rules);

        // Should still pick one based on rules priority
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_max_pctg_region_tiebreaker() {
        let rules = vec![Area::Tss];

        let c1 = make_candidate(Area::Tss, 80.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 90.0, 100.0, "T2", "G1", "1"); // Higher

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].pctg_region, 90.0);
    }

    #[test]
    fn test_same_area_same_pctg_region_tie() {
        let rules = vec![Area::Tss];

        let c1 = make_candidate(Area::Tss, 80.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 80.0, 100.0, "T2", "G1", "1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // Both should be reported (tie)
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_pctg_area_threshold() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Intron, 60.0, 95.0, "T1", "G1", "1"); // High area %
        let c2 = make_candidate(Area::Tss, 60.0, 80.0, "T1", "G1", "1"); // Low area %

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // c1 passes area threshold, c2 doesn't, so c1 should win
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Intron);
    }

    #[test]
    fn test_all_below_area_threshold() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Intron, 60.0, 80.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 60.0, 85.0, "T1", "G1", "1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // Both below area threshold, so fallback to rules priority
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_empty_candidates() {
        let rules = default_rules();
        let candidates: Vec<Candidate> = vec![];
        let grouped_by = AHashMap::new();

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);
        assert!(result.is_empty());
    }

    #[test]
    fn test_multiple_groups() {
        let rules = vec![Area::Tss, Area::Intron];

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c3 = make_candidate(Area::Promoter, 100.0, 100.0, "T2", "G2", "2");

        let candidates = vec![c1, c2, c3];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);
        grouped_by.insert("T2".to_string(), vec![2]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // Should have 2 results - one from each group
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_single_candidate_group() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Downstream, 30.0, 20.0, "T1", "G1", "1");

        let candidates = vec![c1];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // Single candidate is always returned regardless of thresholds
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Downstream);
    }

    #[test]
    fn test_area_not_in_rules() {
        // Rules that don't include all area types
        let rules = vec![Area::Tss, Area::Promoter];

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::GeneBody, 100.0, 100.0, "T1", "G1", "1");

        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("T1".to_string(), vec![0, 1]);

        let result = apply_rules(&candidates, &grouped_by, 50.0, 90.0, &rules);

        // No rules match, so nothing is reported
        assert!(result.is_empty());
    }
}

// -------------------------------------------------------------------------
// 3. Class Logic Tests - Transcript exon numbering
// -------------------------------------------------------------------------

mod test_transcript {
    use super::*;
    use rgmatch::types::Exon;

    #[test]
    fn test_check_exon_numbers_pos() {
        let mut t = Transcript::new("t1".to_string());
        // Add exons in random order
        t.add_exon(Exon::new(500, 600));
        t.add_exon(Exon::new(100, 200));
        t.add_exon(Exon::new(300, 400));

        t.renumber_exons(Strand::Positive);

        // Should be sorted by start: e1, e2, e3
        // Numbering: 1, 2, 3
        assert_eq!(t.exons[0].start, 100);
        assert_eq!(t.exons[0].exon_number, Some("1".to_string()));
        assert_eq!(t.exons[1].start, 300);
        assert_eq!(t.exons[1].exon_number, Some("2".to_string()));
        assert_eq!(t.exons[2].start, 500);
        assert_eq!(t.exons[2].exon_number, Some("3".to_string()));
    }

    #[test]
    fn test_check_exon_numbers_neg() {
        let mut t = Transcript::new("t1".to_string());
        t.add_exon(Exon::new(300, 400));
        t.add_exon(Exon::new(100, 200));

        t.renumber_exons(Strand::Negative);

        // Sorted by start: [100-200, 300-400]
        // For negative strand: first (lowest) gets N, last (highest) gets 1
        assert_eq!(t.exons[0].start, 100);
        assert_eq!(t.exons[0].exon_number, Some("2".to_string()));
        assert_eq!(t.exons[1].start, 300);
        assert_eq!(t.exons[1].exon_number, Some("1".to_string()));
    }
}

// -------------------------------------------------------------------------
// 4. selectTranscript Tests
// -------------------------------------------------------------------------

mod test_select_transcript {
    use super::*;
    use ahash::AHashMap;

    #[test]
    fn test_single_candidate_per_gene() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let candidates = vec![c1];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_multiple_candidates_different_areas() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T2", "G1", "1");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_multiple_candidates_same_area_tie() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 80.0, 70.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 90.0, 60.0, "T2", "G1", "2");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert_eq!(result.len(), 1);
        // Should contain merged transcript info
        assert!(result[0].transcript.contains("T1"));
        assert!(result[0].transcript.contains("T2"));
    }

    #[test]
    fn test_merged_output_pctg_values() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 80.0, 70.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 90.0, 60.0, "T2", "G1", "3");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert_eq!(result[0].pctg_region, 90.0); // max of 80, 90
        assert_eq!(result[0].pctg_area, 70.0); // max of 70, 60
    }

    #[test]
    fn test_merged_exon_numbers() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 100.0, 100.0, "T2", "G1", "3");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert!(result[0].exon_number.contains("1"));
        assert!(result[0].exon_number.contains("3"));
    }

    #[test]
    fn test_empty_candidates() {
        let rules = default_rules();
        let candidates: Vec<Candidate> = vec![];
        let grouped_by = AHashMap::new();

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert!(result.is_empty());
    }

    #[test]
    fn test_multiple_genes() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Intron, 100.0, 100.0, "T2", "G2", "2");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0]);
        grouped_by.insert("G2".to_string(), vec![1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_area_fallback_when_no_rules_match() {
        // Rules that don't include the candidate's area
        let rules = vec![Area::Tss, Area::Promoter];

        let c1 = make_candidate(Area::Intron, 100.0, 100.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::GeneBody, 100.0, 100.0, "T2", "G1", "2");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);
        // Should fallback to first candidate's area
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Intron);
    }

    #[test]
    fn test_merged_transcript_preserves_fields() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 80.0, 70.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 90.0, 60.0, "T2", "G1", "3");
        let candidates = vec![c1, c2];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1]);

        let result = select_transcript(&candidates, &grouped_by, &rules);

        // Verify merged result preserves reference candidate fields
        assert_eq!(result[0].gene, "G1");
        assert_eq!(result[0].area, Area::Tss);
        assert_eq!(result[0].strand, Strand::Positive);
    }

    #[test]
    fn test_three_way_merge() {
        let rules = default_rules();
        let c1 = make_candidate(Area::Tss, 70.0, 60.0, "T1", "G1", "1");
        let c2 = make_candidate(Area::Tss, 80.0, 50.0, "T2", "G1", "2");
        let c3 = make_candidate(Area::Tss, 90.0, 55.0, "T3", "G1", "3");
        let candidates = vec![c1, c2, c3];
        let mut grouped_by = AHashMap::new();
        grouped_by.insert("G1".to_string(), vec![0, 1, 2]);

        let result = select_transcript(&candidates, &grouped_by, &rules);

        assert_eq!(result.len(), 1);
        assert!(result[0].transcript.contains("T1"));
        assert!(result[0].transcript.contains("T2"));
        assert!(result[0].transcript.contains("T3"));
        assert_eq!(result[0].pctg_region, 90.0); // max of 70, 80, 90
        assert_eq!(result[0].pctg_area, 60.0); // max of 60, 50, 55
    }
}

// -------------------------------------------------------------------------
// 5. Config/Rules Tests
// -------------------------------------------------------------------------

mod test_config {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = Config::default();
        assert_eq!(config.rules.len(), 8);
        assert_eq!(config.perc_area, 90.0);
        assert_eq!(config.perc_region, 50.0);
        assert_eq!(config.tss, 200.0);
        assert_eq!(config.tts, 0.0);
        assert_eq!(config.promoter, 1300.0);
        assert_eq!(config.distance, 10000);
    }

    #[test]
    fn test_parse_rules_valid() {
        let mut config = Config::new();
        let result =
            config.parse_rules("DOWNSTREAM,UPSTREAM,GENE_BODY,INTRON,TTS,PROMOTER,1st_EXON,TSS");
        assert!(result);
        assert_eq!(config.rules.len(), 8);
        assert_eq!(config.rules[0], Area::Downstream);
        assert_eq!(config.rules[7], Area::Tss);
    }

    #[test]
    fn test_parse_rules_missing_tags() {
        let mut config = Config::new();
        let result = config.parse_rules("TSS,1st_EXON,PROMOTER");
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_duplicate_tags() {
        let mut config = Config::new();
        let result = config.parse_rules("TSS,TSS,TSS,TSS,TSS,TSS,TSS,TSS");
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_case_sensitive() {
        let mut config = Config::new();
        let result =
            config.parse_rules("tss,1st_exon,promoter,tts,intron,gene_body,upstream,downstream");
        assert!(!result);
    }

    #[test]
    fn test_set_distance_kb() {
        let mut config = Config::new();
        config.set_distance_kb(20);
        assert_eq!(config.distance, 20000);

        config.set_distance_kb(-1);
        assert_eq!(config.distance, 20000); // Should not change
    }

    #[test]
    fn test_set_distance_kb_zero() {
        let mut config = Config::new();
        config.set_distance_kb(0);
        assert_eq!(config.distance, 0);
    }

    #[test]
    fn test_max_lookback_distance_default() {
        let config = Config::default();
        // Default: distance = 10000, tss = 200, tts = 0, promoter = 1300
        // max_lookback = max(10000, max(200, 0, 1300)) = max(10000, 1300) = 10000
        assert_eq!(config.max_lookback_distance(), 10000);
    }

    #[test]
    fn test_max_lookback_distance_large_tss() {
        let mut config = Config::default();
        config.tss = 15000.0;
        // max_lookback = max(10000, max(15000, 0, 1300)) = max(10000, 15000) = 15000
        assert_eq!(config.max_lookback_distance(), 15000);
    }

    #[test]
    fn test_max_lookback_distance_large_promoter() {
        let mut config = Config::default();
        config.promoter = 20000.0;
        // max_lookback = max(10000, max(200, 0, 20000)) = max(10000, 20000) = 20000
        assert_eq!(config.max_lookback_distance(), 20000);
    }

    #[test]
    fn test_max_lookback_distance_large_tts() {
        let mut config = Config::default();
        config.tts = 12000.0;
        // max_lookback = max(10000, max(200, 12000, 1300)) = max(10000, 12000) = 12000
        assert_eq!(config.max_lookback_distance(), 12000);
    }

    #[test]
    fn test_config_new_equals_default() {
        let config_new = Config::new();
        let config_default = Config::default();
        assert_eq!(config_new.rules.len(), config_default.rules.len());
        assert_eq!(config_new.perc_area, config_default.perc_area);
        assert_eq!(config_new.perc_region, config_default.perc_region);
        assert_eq!(config_new.tss, config_default.tss);
        assert_eq!(config_new.tts, config_default.tts);
        assert_eq!(config_new.promoter, config_default.promoter);
        assert_eq!(config_new.distance, config_default.distance);
    }

    #[test]
    fn test_config_default_level() {
        let config = Config::default();
        assert_eq!(config.level, ReportLevel::Exon);
    }

    #[test]
    fn test_config_default_id_tags() {
        let config = Config::default();
        assert_eq!(config.gene_id_tag, "gene_id");
        assert_eq!(config.transcript_id_tag, "transcript_id");
    }

    #[test]
    fn test_parse_rules_preserves_order() {
        let mut config = Config::new();
        let result =
            config.parse_rules("DOWNSTREAM,UPSTREAM,GENE_BODY,INTRON,TTS,PROMOTER,1st_EXON,TSS");
        assert!(result);
        // Verify order is preserved
        assert_eq!(config.rules[0], Area::Downstream);
        assert_eq!(config.rules[1], Area::Upstream);
        assert_eq!(config.rules[2], Area::GeneBody);
        assert_eq!(config.rules[3], Area::Intron);
        assert_eq!(config.rules[4], Area::Tts);
        assert_eq!(config.rules[5], Area::Promoter);
        assert_eq!(config.rules[6], Area::FirstExon);
        assert_eq!(config.rules[7], Area::Tss);
    }

    #[test]
    fn test_parse_rules_with_unknown_tag() {
        let mut config = Config::new();
        // One valid, rest invalid
        let result = config.parse_rules("TSS,UNKNOWN1,UNKNOWN2,UNKNOWN3");
        assert!(!result);
    }

    #[test]
    fn test_parse_rules_empty_string() {
        let mut config = Config::new();
        let result = config.parse_rules("");
        assert!(!result);
    }
}

// -------------------------------------------------------------------------
// 6. Types Module Tests - Strand, Area, Gene, Region, ReportLevel
// -------------------------------------------------------------------------

mod test_types_strand {
    use super::*;

    #[test]
    fn test_strand_from_str_positive() {
        let strand: Strand = "+".parse().unwrap();
        assert_eq!(strand, Strand::Positive);
    }

    #[test]
    fn test_strand_from_str_negative() {
        let strand: Strand = "-".parse().unwrap();
        assert_eq!(strand, Strand::Negative);
    }

    #[test]
    fn test_strand_from_str_invalid() {
        let result = ".".parse::<Strand>();
        assert!(result.is_err());
        let result = "".parse::<Strand>();
        assert!(result.is_err());
        let result = "positive".parse::<Strand>();
        assert!(result.is_err());
    }

    #[test]
    fn test_strand_as_str() {
        assert_eq!(Strand::Positive.as_str(), "+");
        assert_eq!(Strand::Negative.as_str(), "-");
    }

    #[test]
    fn test_strand_display() {
        assert_eq!(format!("{}", Strand::Positive), "+");
        assert_eq!(format!("{}", Strand::Negative), "-");
    }
}

mod test_types_area {
    use super::*;

    #[test]
    fn test_area_from_str_all_variants() {
        assert_eq!("TSS".parse::<Area>().unwrap(), Area::Tss);
        assert_eq!("1st_EXON".parse::<Area>().unwrap(), Area::FirstExon);
        assert_eq!("PROMOTER".parse::<Area>().unwrap(), Area::Promoter);
        assert_eq!("TTS".parse::<Area>().unwrap(), Area::Tts);
        assert_eq!("INTRON".parse::<Area>().unwrap(), Area::Intron);
        assert_eq!("GENE_BODY".parse::<Area>().unwrap(), Area::GeneBody);
        assert_eq!("UPSTREAM".parse::<Area>().unwrap(), Area::Upstream);
        assert_eq!("DOWNSTREAM".parse::<Area>().unwrap(), Area::Downstream);
    }

    #[test]
    fn test_area_from_str_invalid() {
        assert!("INVALID".parse::<Area>().is_err());
        assert!("tss".parse::<Area>().is_err()); // lowercase
        assert!("".parse::<Area>().is_err());
    }

    #[test]
    fn test_area_as_str_all_variants() {
        assert_eq!(Area::Tss.as_str(), "TSS");
        assert_eq!(Area::FirstExon.as_str(), "1st_EXON");
        assert_eq!(Area::Promoter.as_str(), "PROMOTER");
        assert_eq!(Area::Tts.as_str(), "TTS");
        assert_eq!(Area::Intron.as_str(), "INTRON");
        assert_eq!(Area::GeneBody.as_str(), "GENE_BODY");
        assert_eq!(Area::Upstream.as_str(), "UPSTREAM");
        assert_eq!(Area::Downstream.as_str(), "DOWNSTREAM");
    }

    #[test]
    fn test_area_display() {
        assert_eq!(format!("{}", Area::Tss), "TSS");
        assert_eq!(format!("{}", Area::FirstExon), "1st_EXON");
        assert_eq!(format!("{}", Area::Promoter), "PROMOTER");
        assert_eq!(format!("{}", Area::Tts), "TTS");
        assert_eq!(format!("{}", Area::Intron), "INTRON");
        assert_eq!(format!("{}", Area::GeneBody), "GENE_BODY");
        assert_eq!(format!("{}", Area::Upstream), "UPSTREAM");
        assert_eq!(format!("{}", Area::Downstream), "DOWNSTREAM");
    }
}

mod test_types_gene {
    use super::*;
    use rgmatch::Gene;
    use rgmatch::types::Exon;

    #[test]
    fn test_gene_new() {
        let gene = Gene::new("GENE001".to_string(), Strand::Positive);
        assert_eq!(gene.gene_id, "GENE001");
        assert_eq!(gene.strand, Strand::Positive);
        assert!(gene.transcripts.is_empty());
        assert_eq!(gene.start, i64::MAX);
        assert_eq!(gene.end, 0);
    }

    #[test]
    fn test_gene_add_transcript() {
        let mut gene = Gene::new("GENE001".to_string(), Strand::Positive);
        let transcript = Transcript::new("T1".to_string());
        gene.add_transcript(transcript);
        assert_eq!(gene.transcripts.len(), 1);
        assert_eq!(gene.transcripts[0].transcript_id, "T1");
    }

    #[test]
    fn test_gene_set_length() {
        let mut gene = Gene::new("GENE001".to_string(), Strand::Positive);
        gene.set_length(1000, 2000);
        assert_eq!(gene.start, 1000);
        assert_eq!(gene.end, 2000);
    }

    #[test]
    fn test_gene_calculate_size() {
        let mut gene = Gene::new("GENE001".to_string(), Strand::Positive);

        let mut t1 = Transcript::new("T1".to_string());
        t1.add_exon(Exon::new(100, 200));
        t1.calculate_size();

        let mut t2 = Transcript::new("T2".to_string());
        t2.add_exon(Exon::new(50, 300));
        t2.calculate_size();

        gene.add_transcript(t1);
        gene.add_transcript(t2);
        gene.calculate_size();

        assert_eq!(gene.start, 50);
        assert_eq!(gene.end, 300);
    }

    #[test]
    fn test_gene_calculate_size_empty() {
        let mut gene = Gene::new("GENE001".to_string(), Strand::Positive);
        gene.calculate_size();
        // Empty gene keeps initial values
        assert_eq!(gene.start, i64::MAX);
        assert_eq!(gene.end, 0);
    }
}

mod test_types_report_level {
    use super::*;

    #[test]
    fn test_report_level_from_str() {
        assert_eq!("exon".parse::<ReportLevel>().unwrap(), ReportLevel::Exon);
        assert_eq!(
            "transcript".parse::<ReportLevel>().unwrap(),
            ReportLevel::Transcript
        );
        assert_eq!("gene".parse::<ReportLevel>().unwrap(), ReportLevel::Gene);
    }

    #[test]
    fn test_report_level_case_insensitive() {
        assert_eq!("EXON".parse::<ReportLevel>().unwrap(), ReportLevel::Exon);
        assert_eq!("Exon".parse::<ReportLevel>().unwrap(), ReportLevel::Exon);
        assert_eq!(
            "TRANSCRIPT".parse::<ReportLevel>().unwrap(),
            ReportLevel::Transcript
        );
        assert_eq!("GENE".parse::<ReportLevel>().unwrap(), ReportLevel::Gene);
    }

    #[test]
    fn test_report_level_from_str_invalid() {
        assert!("invalid".parse::<ReportLevel>().is_err());
        assert!("".parse::<ReportLevel>().is_err());
        assert!("region".parse::<ReportLevel>().is_err());
    }
}

mod test_types_exon {
    use rgmatch::types::Exon;

    #[test]
    fn test_exon_new() {
        let exon = Exon::new(100, 200);
        assert_eq!(exon.start, 100);
        assert_eq!(exon.end, 200);
        assert!(exon.exon_number.is_none());
    }

    #[test]
    fn test_exon_length() {
        let exon = Exon::new(100, 200);
        // length = end - start + 1 = 200 - 100 + 1 = 101
        assert_eq!(exon.length(), 101);
    }

    #[test]
    fn test_exon_length_single_base() {
        let exon = Exon::new(100, 100);
        assert_eq!(exon.length(), 1);
    }
}

mod test_types_transcript {
    use super::*;
    use rgmatch::types::Exon;

    #[test]
    fn test_transcript_new() {
        let transcript = Transcript::new("T1".to_string());
        assert_eq!(transcript.transcript_id, "T1");
        assert!(transcript.exons.is_empty());
        assert_eq!(transcript.start, i64::MAX);
        assert_eq!(transcript.end, 0);
    }

    #[test]
    fn test_transcript_add_exon() {
        let mut transcript = Transcript::new("T1".to_string());
        transcript.add_exon(Exon::new(100, 200));
        assert_eq!(transcript.exons.len(), 1);
    }

    #[test]
    fn test_transcript_set_length() {
        let mut transcript = Transcript::new("T1".to_string());
        transcript.set_length(50, 500);
        assert_eq!(transcript.start, 50);
        assert_eq!(transcript.end, 500);
    }

    #[test]
    fn test_transcript_calculate_size() {
        let mut transcript = Transcript::new("T1".to_string());
        transcript.add_exon(Exon::new(100, 200));
        transcript.add_exon(Exon::new(300, 400));
        transcript.add_exon(Exon::new(50, 150));
        transcript.calculate_size();

        assert_eq!(transcript.start, 50);
        assert_eq!(transcript.end, 400);
    }

    #[test]
    fn test_transcript_calculate_size_empty() {
        let mut transcript = Transcript::new("T1".to_string());
        transcript.calculate_size();
        // Empty transcript keeps initial values
        assert_eq!(transcript.start, i64::MAX);
        assert_eq!(transcript.end, 0);
    }

    #[test]
    fn test_transcript_renumber_exons_positive_unsorted() {
        let mut transcript = Transcript::new("T1".to_string());
        // Add exons in random order
        transcript.add_exon(Exon::new(500, 600));
        transcript.add_exon(Exon::new(100, 200));
        transcript.add_exon(Exon::new(300, 400));

        transcript.renumber_exons(Strand::Positive);

        // Should be sorted by start: [100-200, 300-400, 500-600]
        assert_eq!(transcript.exons[0].start, 100);
        assert_eq!(transcript.exons[0].exon_number, Some("1".to_string()));
        assert_eq!(transcript.exons[1].start, 300);
        assert_eq!(transcript.exons[1].exon_number, Some("2".to_string()));
        assert_eq!(transcript.exons[2].start, 500);
        assert_eq!(transcript.exons[2].exon_number, Some("3".to_string()));
    }

    #[test]
    fn test_transcript_renumber_exons_negative_unsorted() {
        let mut transcript = Transcript::new("T1".to_string());
        transcript.add_exon(Exon::new(500, 600));
        transcript.add_exon(Exon::new(100, 200));
        transcript.add_exon(Exon::new(300, 400));

        transcript.renumber_exons(Strand::Negative);

        // After sorting: [100-200, 300-400, 500-600]
        // For negative strand: first (lowest) gets N, last (highest) gets 1
        assert_eq!(transcript.exons[0].start, 100);
        assert_eq!(transcript.exons[0].exon_number, Some("3".to_string()));
        assert_eq!(transcript.exons[1].start, 300);
        assert_eq!(transcript.exons[1].exon_number, Some("2".to_string()));
        assert_eq!(transcript.exons[2].start, 500);
        assert_eq!(transcript.exons[2].exon_number, Some("1".to_string()));
    }

    #[test]
    fn test_transcript_renumber_single_exon() {
        let mut transcript = Transcript::new("T1".to_string());
        transcript.add_exon(Exon::new(100, 200));

        transcript.renumber_exons(Strand::Positive);
        assert_eq!(transcript.exons[0].exon_number, Some("1".to_string()));

        transcript.renumber_exons(Strand::Negative);
        assert_eq!(transcript.exons[0].exon_number, Some("1".to_string()));
    }
}

mod test_types_region {
    use rgmatch::Region;

    #[test]
    fn test_region_new() {
        let region = Region::new(
            "chr1".to_string(),
            100,
            200,
            vec!["name1".to_string(), "500".to_string()],
        );
        assert_eq!(region.chrom, "chr1");
        assert_eq!(region.start, 100);
        assert_eq!(region.end, 200);
        assert_eq!(region.metadata.len(), 2);
    }

    #[test]
    fn test_region_length() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        // length = end - start + 1 = 200 - 100 + 1 = 101
        assert_eq!(region.length(), 101);
    }

    #[test]
    fn test_region_midpoint() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        // midpoint = (100 + 200) / 2 = 150
        assert_eq!(region.midpoint(), 150);
    }

    #[test]
    fn test_region_midpoint_integer_division() {
        let region = Region::new("chr1".to_string(), 100, 201, vec![]);
        // midpoint = (100 + 201) / 2 = 150 (integer division)
        assert_eq!(region.midpoint(), 150);
    }

    #[test]
    fn test_region_id() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        assert_eq!(region.id(), "chr1_100_200");
    }

    #[test]
    fn test_region_id_with_special_chrom() {
        let region = Region::new("chrUn_gl000220".to_string(), 1000, 2000, vec![]);
        assert_eq!(region.id(), "chrUn_gl000220_1000_2000");
    }
}

mod test_types_candidate {
    use super::*;

    #[test]
    fn test_candidate_new() {
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            50,
            80.5,
            90.5,
            100,
        );
        assert_eq!(candidate.start, 100);
        assert_eq!(candidate.end, 200);
        assert_eq!(candidate.strand, Strand::Positive);
        assert_eq!(candidate.exon_number, "1");
        assert_eq!(candidate.area, Area::Tss);
        assert_eq!(candidate.transcript, "T1");
        assert_eq!(candidate.gene, "G1");
        assert_eq!(candidate.distance, 50);
        assert_eq!(candidate.pctg_region, 80.5);
        assert_eq!(candidate.pctg_area, 90.5);
        assert_eq!(candidate.tss_distance, 100);
    }

    #[test]
    fn test_candidate_negative_distance() {
        let candidate = Candidate::new(
            100,
            200,
            Strand::Negative,
            "2".to_string(),
            Area::Upstream,
            "T1".to_string(),
            "G1".to_string(),
            -500,
            100.0,
            -1.0,
            -1000,
        );
        assert_eq!(candidate.distance, -500);
        assert_eq!(candidate.pctg_area, -1.0);
        assert_eq!(candidate.tss_distance, -1000);
    }
}

// -------------------------------------------------------------------------
// 7. Bug Regression Tests
// -------------------------------------------------------------------------

mod test_bug_regression {
    use super::*;
    use rgmatch::matcher::overlap::match_region_to_genes;
    use rgmatch::types::Exon;
    use rgmatch::{Gene, Region};
    use std::collections::HashSet;

    fn make_test_gene(
        gene_id: &str,
        start: i64,
        end: i64,
        strand: Strand,
        exons: Vec<(i64, i64)>,
    ) -> Gene {
        let mut gene = Gene::new(gene_id.to_string(), strand);
        gene.set_length(start, end);
        let mut transcript = Transcript::new(format!("TRANS_{}", gene_id.replace("GENE", "")));
        for (i, (exon_start, exon_end)) in exons.iter().enumerate() {
            let mut exon = Exon::new(*exon_start, *exon_end);
            exon.exon_number = Some((i + 1).to_string());
            transcript.add_exon(exon);
        }
        transcript.calculate_size();
        transcript.renumber_exons(strand);
        gene.transcripts.push(transcript);
        gene
    }

    /// Bug #1: Test that Case 2 (partial overlap left) doesn't produce duplicate DOWNSTREAM
    #[test]
    fn test_no_duplicate_downstream_case2() {
        // Region [100, 200) partially overlaps exon [51, 150] on the left
        let config = Config::default();
        let region = Region::new("chr1".into(), 100, 200, vec!["region1".into()]);

        // Single-exon gene - triggers Case 2 (partial overlap left on last exon)
        let genes = vec![make_test_gene(
            "GENE001",
            51,
            150,
            Strand::Positive,
            vec![(51, 150)],
        )];

        let last_index = 0;
        let candidates = match_region_to_genes(&region, &genes, &config, last_index);

        // Count DOWNSTREAM candidates for GENE001
        let downstream_count = candidates
            .iter()
            .filter(|c| c.gene == "GENE001" && c.area == Area::Downstream)
            .count();

        assert_eq!(
            downstream_count, 1,
            "GENE001 DOWNSTREAM should appear exactly once, not {}",
            downstream_count
        );
    }

    /// Bug #1: Test that Case 3 (exon inside region) doesn't produce duplicate DOWNSTREAM
    #[test]
    fn test_no_duplicate_downstream_case3() {
        // Region [1000, 1300) completely contains exon [1050, 1200]
        let config = Config::default();
        let region = Region::new("chr1".into(), 1000, 1300, vec!["region2".into()]);

        let genes = vec![make_test_gene(
            "GENE002",
            1050,
            1200,
            Strand::Positive,
            vec![(1050, 1200)],
        )];

        let last_index = 0;
        let candidates = match_region_to_genes(&region, &genes, &config, last_index);

        let downstream_count = candidates
            .iter()
            .filter(|c| c.gene == "GENE002" && c.area == Area::Downstream)
            .count();

        assert_eq!(
            downstream_count, 1,
            "GENE002 DOWNSTREAM should appear exactly once, not {}",
            downstream_count
        );
    }

    /// Bug #2: Test that proximity candidates are preserved when overlapping gene comes later
    #[test]
    fn test_proximity_candidate_preserved() {
        // Region [5000, 5100)
        // GENE003: multi-exon gene ending at 4900 (proximity candidate)
        // GENE004: multi-exon gene with exon overlapping region (overlapping candidate)
        let config = Config::default();
        let region = Region::new("chr1".into(), 5000, 5100, vec!["region3".into()]);

        let genes = vec![
            // GENE003: ends at 4900, 100bp before region - should be DOWNSTREAM proximity
            make_test_gene(
                "GENE003",
                4700,
                4900,
                Strand::Positive,
                vec![(4700, 4750), (4800, 4900)],
            ),
            // GENE004: exon 2 overlaps region at [4950, 5050]
            make_test_gene(
                "GENE004",
                4850,
                5200,
                Strand::Positive,
                vec![(4850, 4900), (4950, 5050)],
            ),
        ];

        let last_index = 0;
        let candidates = match_region_to_genes(&region, &genes, &config, last_index);

        // GENE003 DOWNSTREAM should be preserved (proximity candidate)
        let gene003_downstream = candidates
            .iter()
            .any(|c| c.gene == "GENE003" && c.area == Area::Downstream);

        // GENE004 should also have candidates (overlapping)
        let gene004_present = candidates.iter().any(|c| c.gene == "GENE004");

        assert!(
            gene003_downstream,
            "GENE003 DOWNSTREAM proximity candidate should be preserved"
        );
        assert!(
            gene004_present,
            "GENE004 overlapping candidate should be present"
        );
    }

    /// Combined test: no duplicates across all test regions
    #[test]
    fn test_no_duplicate_lines_overall() {
        let config = Config::default();

        let regions = vec![
            Region::new("chr1".into(), 100, 200, vec!["region1".into()]),
            Region::new("chr1".into(), 1000, 1300, vec!["region2".into()]),
            Region::new("chr1".into(), 5000, 5100, vec!["region3".into()]),
        ];

        let genes = vec![
            make_test_gene("GENE001", 51, 150, Strand::Positive, vec![(51, 150)]),
            make_test_gene("GENE002", 1050, 1200, Strand::Positive, vec![(1050, 1200)]),
            make_test_gene(
                "GENE003",
                4700,
                4900,
                Strand::Positive,
                vec![(4700, 4750), (4800, 4900)],
            ),
            make_test_gene(
                "GENE004",
                4850,
                5200,
                Strand::Positive,
                vec![(4850, 4900), (4950, 5050)],
            ),
        ];

        for region in &regions {
            let last_index = 0;
            let candidates = match_region_to_genes(region, &genes, &config, last_index);

            // Create unique key for each candidate
            let keys: Vec<String> = candidates
                .iter()
                .map(|c| format!("{}_{}_{}", c.gene, c.transcript, c.area))
                .collect();
            let unique_keys: HashSet<_> = keys.iter().collect();

            assert_eq!(
                keys.len(),
                unique_keys.len(),
                "Duplicate candidates found for region {:?}",
                region.id()
            );
        }
    }
}

// -------------------------------------------------------------------------
// 8. Overlap Module Tests - process_candidates_for_output, match_regions_to_genes, find_search_start_index
// -------------------------------------------------------------------------

mod test_parser_gtf_basic {
    // Placeholder - actual GTF tests are in test_parser_gtf module at the end of this file

    #[test]
    fn test_gtf_extract_attribute_edge_cases() {
        // Test that attributes with various formats are handled correctly
        // This tests the public interface behavior

        // Note: The actual extract_attribute function is private, but we can
        // test the parsing behavior through the main parse functions
        // For now, we verify the parsing works through integration-style tests
        // that check the public API behaves correctly
    }
}

mod test_overlap_functions {
    use super::*;
    use rgmatch::types::Exon;
    use rgmatch::{Gene, Region};

    fn make_test_gene(
        gene_id: &str,
        start: i64,
        end: i64,
        strand: Strand,
        exons: Vec<(i64, i64)>,
    ) -> Gene {
        let mut gene = Gene::new(gene_id.to_string(), strand);
        gene.set_length(start, end);
        let mut transcript = Transcript::new(format!("TRANS_{}", gene_id.replace("GENE", "")));
        for (i, (exon_start, exon_end)) in exons.iter().enumerate() {
            let mut exon = Exon::new(*exon_start, *exon_end);
            exon.exon_number = Some((i + 1).to_string());
            transcript.add_exon(exon);
        }
        transcript.calculate_size();
        transcript.renumber_exons(strand);
        gene.transcripts.push(transcript);
        gene
    }

    #[test]
    fn test_find_search_start_index_empty() {
        let genes: Vec<Gene> = vec![];
        assert_eq!(find_search_start_index(&genes, 1000), 0);
    }

    #[test]
    fn test_find_search_start_index_all_before() {
        let genes = vec![
            make_test_gene("G1", 100, 200, Strand::Positive, vec![(100, 200)]),
            make_test_gene("G2", 300, 400, Strand::Positive, vec![(300, 400)]),
        ];
        // All genes start before 1000
        assert_eq!(find_search_start_index(&genes, 1000), 2);
    }

    #[test]
    fn test_find_search_start_index_all_after() {
        let genes = vec![
            make_test_gene("G1", 1000, 2000, Strand::Positive, vec![(1000, 2000)]),
            make_test_gene("G2", 3000, 4000, Strand::Positive, vec![(3000, 4000)]),
        ];
        // All genes start at or after 100
        assert_eq!(find_search_start_index(&genes, 100), 0);
    }

    #[test]
    fn test_find_search_start_index_middle() {
        let genes = vec![
            make_test_gene("G1", 100, 200, Strand::Positive, vec![(100, 200)]),
            make_test_gene("G2", 500, 600, Strand::Positive, vec![(500, 600)]),
            make_test_gene("G3", 1000, 1100, Strand::Positive, vec![(1000, 1100)]),
        ];
        // First gene with start >= 500
        assert_eq!(find_search_start_index(&genes, 500), 1);
        // First gene with start >= 600
        assert_eq!(find_search_start_index(&genes, 600), 2);
    }

    #[test]
    fn test_process_candidates_for_output_empty() {
        let config = Config::default();
        let candidates: Vec<Candidate> = vec![];
        let result = process_candidates_for_output(candidates, &config);
        assert!(result.is_empty());
    }

    #[test]
    fn test_process_candidates_for_output_exon_level() {
        let mut config = Config::default();
        config.level = ReportLevel::Exon;

        let candidates = vec![
            Candidate::new(
                100,
                200,
                Strand::Positive,
                "1".to_string(),
                Area::Tss,
                "T1".to_string(),
                "G1".to_string(),
                0,
                80.0,
                90.0,
                100,
            ),
            Candidate::new(
                300,
                400,
                Strand::Positive,
                "2".to_string(),
                Area::Intron,
                "T1".to_string(),
                "G1".to_string(),
                0,
                60.0,
                70.0,
                200,
            ),
        ];

        let result = process_candidates_for_output(candidates.clone(), &config);
        // Exon level returns all candidates
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_process_candidates_for_output_transcript_level() {
        let mut config = Config::default();
        config.level = ReportLevel::Transcript;

        let candidates = vec![
            Candidate::new(
                100,
                200,
                Strand::Positive,
                "1".to_string(),
                Area::Tss,
                "T1".to_string(),
                "G1".to_string(),
                0,
                80.0,
                95.0,
                100,
            ),
            Candidate::new(
                300,
                400,
                Strand::Positive,
                "2".to_string(),
                Area::Intron,
                "T1".to_string(),
                "G1".to_string(),
                0,
                60.0,
                95.0,
                200,
            ),
        ];

        let result = process_candidates_for_output(candidates, &config);
        // Transcript level returns best per transcript (TSS has higher priority)
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].area, Area::Tss);
    }

    #[test]
    fn test_process_candidates_for_output_gene_level() {
        let mut config = Config::default();
        config.level = ReportLevel::Gene;

        let candidates = vec![
            Candidate::new(
                100,
                200,
                Strand::Positive,
                "1".to_string(),
                Area::Tss,
                "T1".to_string(),
                "G1".to_string(),
                0,
                80.0,
                95.0,
                100,
            ),
            Candidate::new(
                300,
                400,
                Strand::Positive,
                "2".to_string(),
                Area::Tss,
                "T2".to_string(),
                "G1".to_string(),
                0,
                90.0,
                95.0,
                200,
            ),
        ];

        let result = process_candidates_for_output(candidates, &config);
        // Gene level merges candidates with same area
        assert_eq!(result.len(), 1);
        assert!(result[0].transcript.contains("T1"));
        assert!(result[0].transcript.contains("T2"));
    }

    #[test]
    fn test_match_regions_to_genes_empty_regions() {
        let config = Config::default();
        let regions: Vec<Region> = vec![];
        let genes = vec![make_test_gene(
            "G1",
            100,
            200,
            Strand::Positive,
            vec![(100, 200)],
        )];

        let results = match_regions_to_genes(&regions, &genes, &config, 500);
        assert!(results.is_empty());
    }

    #[test]
    fn test_match_regions_to_genes_empty_genes() {
        let config = Config::default();
        let regions = vec![Region::new("chr1".into(), 100, 200, vec![])];
        let genes: Vec<Gene> = vec![];

        let results = match_regions_to_genes(&regions, &genes, &config, 0);
        assert_eq!(results.len(), 1);
        assert!(results[0].1.is_empty()); // No candidates
    }

    #[test]
    fn test_match_regions_to_genes_single_match() {
        let config = Config::default();
        let regions = vec![Region::new("chr1".into(), 100, 200, vec![])];
        let genes = vec![make_test_gene(
            "G1",
            50,
            250,
            Strand::Positive,
            vec![(50, 250)],
        )];

        let results = match_regions_to_genes(&regions, &genes, &config, 500);
        assert_eq!(results.len(), 1);
        assert!(!results[0].1.is_empty());
    }

    #[test]
    fn test_match_regions_to_genes_preserves_region() {
        let config = Config::default();
        let regions = vec![Region::new(
            "chr1".into(),
            100,
            200,
            vec!["metadata1".into()],
        )];
        let genes = vec![make_test_gene(
            "G1",
            50,
            250,
            Strand::Positive,
            vec![(50, 250)],
        )];

        let results = match_regions_to_genes(&regions, &genes, &config, 500);
        assert_eq!(results[0].0.chrom, "chr1");
        assert_eq!(results[0].0.start, 100);
        assert_eq!(results[0].0.end, 200);
        assert_eq!(results[0].0.metadata[0], "metadata1");
    }

    #[test]
    fn test_match_region_to_genes_intron() {
        let config = Config::default();
        // Region inside intron
        let region = Region::new("chr1".into(), 250, 350, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            100,
            500,
            Strand::Positive,
            vec![(100, 200), (400, 500)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        let has_intron = candidates.iter().any(|c| c.area == Area::Intron);
        assert!(has_intron, "Should find intron overlap");
    }

    #[test]
    fn test_match_region_to_genes_first_exon() {
        let config = Config::default();
        // Region overlapping first exon
        let region = Region::new("chr1".into(), 90, 150, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            100,
            500,
            Strand::Positive,
            vec![(100, 200), (400, 500)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        let has_first_exon = candidates.iter().any(|c| c.area == Area::FirstExon);
        assert!(has_first_exon, "Should find first exon overlap");
    }

    #[test]
    fn test_match_region_to_genes_gene_body() {
        let config = Config::default();
        // Region overlapping non-first exon
        let region = Region::new("chr1".into(), 380, 450, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            100,
            500,
            Strand::Positive,
            vec![(100, 200), (400, 500)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        let has_gene_body = candidates.iter().any(|c| c.area == Area::GeneBody);
        assert!(has_gene_body, "Should find gene body overlap");
    }

    #[test]
    fn test_match_region_to_genes_negative_strand() {
        let config = Config::default();
        // Region overlapping first exon on negative strand (last in genomic order)
        let region = Region::new("chr1".into(), 380, 450, vec![]);
        let genes = vec![make_test_gene(
            "G1",
            100,
            500,
            Strand::Negative,
            vec![(100, 200), (400, 500)],
        )];

        let candidates = match_region_to_genes(&region, &genes, &config, 0);
        let has_first_exon = candidates.iter().any(|c| c.area == Area::FirstExon);
        assert!(
            has_first_exon,
            "Should find first exon (last genomic) for negative strand"
        );
    }
}

// -------------------------------------------------------------------------
// 9. Output Module Tests - write_header, format_output_line edge cases
// -------------------------------------------------------------------------

mod test_output {
    use rgmatch::output::{format_output_line, write_header};
    use rgmatch::types::{Area, Candidate, Region, Strand};

    #[test]
    fn test_format_output_line_empty_metadata() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            50,
            80.00,
            90.00,
            500,
        );

        let line = format_output_line(&region, &candidate);

        // Should NOT end with a tab
        assert!(!line.ends_with('\t'));
        // Should end with the pctg_area value
        assert!(line.ends_with("90.00"));
        // Count tabs - should be exactly 9 (10 fields = 9 separators)
        assert_eq!(line.matches('\t').count(), 9);
    }

    #[test]
    fn test_format_output_line_large_coordinates() {
        // Human chr1 is ~249 million bases
        let region = Region::new(
            "chr1".to_string(),
            248956422, // Large start
            248956500, // Large end
            vec![],
        );
        let candidate = Candidate::new(
            248956000,
            248957000,
            Strand::Positive,
            "1".to_string(),
            Area::GeneBody,
            "ENST00000000001".to_string(),
            "ENSG00000000001".to_string(),
            1000000,  // Large distance
            50.0,
            25.0,
            5000000, // Large TSS distance
        );

        let line = format_output_line(&region, &candidate);

        // Verify large coordinates appear correctly
        assert!(line.contains("chr1_248956422_248956500"));
        // Midpoint: (248956422 + 248956500) / 2 = 248956461
        assert!(line.contains("248956461"));
        assert!(line.contains("1000000"));
        assert!(line.contains("5000000"));
    }

    #[test]
    fn test_format_output_line_special_chars_in_metadata() {
        let region = Region::new(
            "chr1".to_string(),
            100,
            200,
            vec![
                "peak_1;gene=ABC".to_string(),
                "score=0.95".to_string(),
                "name with spaces".to_string(),
            ],
        );
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            0,
            100.0,
            100.0,
            0,
        );

        let line = format_output_line(&region, &candidate);

        // Verify special characters are preserved
        assert!(line.contains("peak_1;gene=ABC"));
        assert!(line.contains("score=0.95"));
        assert!(line.contains("name with spaces"));
    }

    #[test]
    fn test_format_output_line_multiple_metadata_columns() {
        let region = Region::new(
            "chr1".to_string(),
            100,
            200,
            vec![
                "name1".to_string(),
                "500".to_string(),
                "+".to_string(),
                "100".to_string(),
                "200".to_string(),
            ],
        );
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            0,
            100.0,
            100.0,
            0,
        );

        let line = format_output_line(&region, &candidate);

        // 10 base fields + 5 metadata = 15 fields = 14 tabs
        assert_eq!(line.matches('\t').count(), 14);

        // Verify metadata appears in correct order
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 15);
        assert_eq!(fields[10], "name1");
        assert_eq!(fields[11], "500");
        assert_eq!(fields[12], "+");
        assert_eq!(fields[13], "100");
        assert_eq!(fields[14], "200");
    }

    #[test]
    fn test_format_output_line_zero_length_region() {
        let region = Region::new("chr1".to_string(), 100, 100, vec![]);
        let candidate = Candidate::new(
            100,
            100,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            0,
            100.0,
            100.0,
            0,
        );

        let line = format_output_line(&region, &candidate);

        // Region ID should be chr1_100_100
        assert!(line.contains("chr1_100_100"));
        // Midpoint of 100-100 is 100
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[0], "chr1_100_100");
        assert_eq!(fields[1], "100"); // midpoint
    }

    #[test]
    fn test_format_output_line_edge_percentages() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);

        // Test 0% values
        let candidate_zero = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Upstream,
            "T1".to_string(),
            "G1".to_string(),
            1000,
            0.0,
            0.0,
            500,
        );
        let line_zero = format_output_line(&region, &candidate_zero);
        let fields: Vec<&str> = line_zero.split('\t').collect();
        assert_eq!(fields[8], "0.00");
        assert_eq!(fields[9], "0.00");

        // Test 100% values
        let candidate_full = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            0,
            100.0,
            100.0,
            0,
        );
        let line_full = format_output_line(&region, &candidate_full);
        assert!(line_full.contains("100.00"));
    }

    #[test]
    fn test_format_output_line_all_area_types() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);

        let area_types = [
            (Area::Tss, "TSS"),
            (Area::FirstExon, "1st_EXON"),
            (Area::Promoter, "PROMOTER"),
            (Area::Tts, "TTS"),
            (Area::Intron, "INTRON"),
            (Area::GeneBody, "GENE_BODY"),
            (Area::Upstream, "UPSTREAM"),
            (Area::Downstream, "DOWNSTREAM"),
        ];

        for (area, expected_str) in area_types {
            let candidate = Candidate::new(
                100,
                200,
                Strand::Positive,
                "1".to_string(),
                area,
                "T1".to_string(),
                "G1".to_string(),
                0,
                50.0,
                50.0,
                0,
            );

            let line = format_output_line(&region, &candidate);
            assert!(
                line.contains(&format!("\t{}\t", expected_str)),
                "Area {:?} should produce '{}' in output, got: {}",
                area,
                expected_str,
                line
            );
        }
    }

    #[test]
    fn test_format_output_line_negative_distance() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Negative,
            "1".to_string(),
            Area::Downstream,
            "T1".to_string(),
            "G1".to_string(),
            -500,  // Negative distance
            50.0,
            50.0,
            -1000, // Negative TSS distance
        );

        let line = format_output_line(&region, &candidate);

        assert!(line.contains("\t-500\t"));
        assert!(line.contains("\t-1000\t"));
    }

    #[test]
    fn test_write_header_single_meta_column() {
        let mut output = Vec::new();
        write_header(&mut output, 1).unwrap();
        let header = String::from_utf8(output).unwrap();

        assert!(header.contains("PercArea\tname"));
        assert!(!header.contains("score")); // Should NOT have second column
    }

    #[test]
    fn test_format_output_line_very_small_percentages() {
        let region = Region::new("chr1".to_string(), 100, 200, vec![]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Upstream,
            "T1".to_string(),
            "G1".to_string(),
            1000,
            0.001, // Very small, should round to 0.00
            0.009, // Very small, should round to 0.01
            500,
        );

        let line = format_output_line(&region, &candidate);

        // 0.001 rounds to 0.00
        // 0.009 rounds to 0.01
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields[8], "0.00"); // pctg_region
        assert_eq!(fields[9], "0.01"); // pctg_area
    }

    #[test]
    fn test_format_output_line_field_count() {
        // Test with 0 metadata
        let region0 = Region::new("chr1".to_string(), 100, 200, vec![]);
        let candidate = Candidate::new(
            100,
            200,
            Strand::Positive,
            "1".to_string(),
            Area::Tss,
            "T1".to_string(),
            "G1".to_string(),
            0,
            100.0,
            100.0,
            0,
        );
        let line0 = format_output_line(&region0, &candidate);
        assert_eq!(line0.split('\t').count(), 10);

        // Test with 3 metadata
        let region3 = Region::new(
            "chr1".to_string(),
            100,
            200,
            vec!["a".to_string(), "b".to_string(), "c".to_string()],
        );
        let line3 = format_output_line(&region3, &candidate);
        assert_eq!(line3.split('\t').count(), 13);

        // Test with 1 metadata
        let region1 = Region::new("chr1".to_string(), 100, 200, vec!["single".to_string()]);
        let line1 = format_output_line(&region1, &candidate);
        assert_eq!(line1.split('\t').count(), 11);
    }

    #[test]
    fn test_write_header_ends_with_newline() {
        let mut output = Vec::new();
        write_header(&mut output, 2).unwrap();
        let header = String::from_utf8(output).unwrap();

        assert!(header.ends_with('\n'));
        assert!(!header.ends_with("\n\n"));

        // Count newlines - should be exactly 1
        assert_eq!(header.matches('\n').count(), 1);
    }
}

// -------------------------------------------------------------------------
// 10. Parser BED Module Tests - BedReader, parse_bed edge cases
// -------------------------------------------------------------------------

mod test_parser_bed {
    use rgmatch::parser::bed::get_bed_headers;

    #[test]
    fn test_parse_bed_insufficient_columns() {
        // Lines with < 3 columns should be silently skipped
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1\t100").unwrap(); // 2 columns - skip
        writeln!(temp_file, "chr1").unwrap(); // 1 column - skip
        writeln!(temp_file, "chr1\t100\t200").unwrap(); // valid
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(100).unwrap().unwrap();

        assert_eq!(chunk.len(), 1);
        assert_eq!(chunk[0].chrom, "chr1");
        assert_eq!(chunk[0].start, 100);
    }

    #[test]
    fn test_parse_bed_metadata_truncation() {
        // BED line with 15 columns (3 required + 12 metadata)
        // Should truncate metadata to 9 columns
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(
            temp_file,
            "chr1\t0\t100\tm1\tm2\tm3\tm4\tm5\tm6\tm7\tm8\tm9\tm10\tm11\tm12"
        )
        .unwrap();
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap().unwrap();

        assert_eq!(chunk[0].metadata.len(), 9);
        assert_eq!(chunk[0].metadata[8], "m9"); // Last kept column
        assert_eq!(reader.num_meta_columns(), 9);
    }

    #[test]
    fn test_parse_bed_large_coordinates() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Test with large but valid i64 coordinates
        let large_start: i64 = 9_000_000_000;
        let large_end: i64 = 9_000_000_100;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1\t{}\t{}", large_start, large_end).unwrap();
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap().unwrap();

        assert_eq!(chunk[0].start, large_start);
        assert_eq!(chunk[0].end, large_end);
    }

    #[test]
    fn test_parse_bed_negative_coordinates() {
        // The parser uses i64, so negative coordinates should parse successfully
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1\t-100\t200").unwrap();
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap().unwrap();

        assert_eq!(chunk[0].start, -100);
        assert_eq!(chunk[0].end, 200);
    }

    #[test]
    fn test_parse_bed_scientific_notation_rejected() {
        // Scientific notation should fail to parse as i64 and line should be skipped
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1\t1e6\t2e6\tscientific").unwrap();
        writeln!(temp_file, "chr1\t1000000\t2000000\tvalid").unwrap();
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap().unwrap();

        // Only the valid integer line should be parsed
        assert_eq!(chunk.len(), 1);
        assert_eq!(chunk[0].start, 1000000);
        assert_eq!(chunk[0].metadata[0], "valid");
    }

    #[test]
    fn test_parse_bed_space_not_delimiter() {
        // The parser uses split('\t') - spaces should NOT be treated as delimiters
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1 with spaces\t100\t200\tname with spaces").unwrap();
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap().unwrap();

        // "chr1 with spaces" should be the chromosome name (spaces preserved)
        assert_eq!(chunk[0].chrom, "chr1 with spaces");
        assert_eq!(chunk[0].metadata[0], "name with spaces");
    }

    #[test]
    fn test_parse_bed_empty_file() {
        use tempfile::NamedTempFile;

        let temp_file = NamedTempFile::new().unwrap();
        // Don't write anything - empty file

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap();

        // Should return None for empty file (EOF immediately)
        assert!(chunk.is_none());
        assert_eq!(reader.num_meta_columns(), 0);
    }

    #[test]
    fn test_parse_bed_all_invalid_lines() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "header\tstart\tend").unwrap();
        writeln!(temp_file, "chr1").unwrap();
        writeln!(temp_file, "chr1\tabc\t200").unwrap(); // non-numeric start
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap();

        // No valid regions, should return None
        assert!(chunk.is_none());
    }

    #[test]
    fn test_get_bed_headers_beyond_nine() {
        // Requesting more than 9 headers should return exactly 9 (all available)
        let headers = get_bed_headers(15);
        assert_eq!(headers.len(), 9);
        assert_eq!(headers[8], "blockStarts");
    }

    #[test]
    fn test_get_bed_headers_exactly_nine() {
        let headers = get_bed_headers(9);
        assert_eq!(headers.len(), 9);
        assert_eq!(
            headers,
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
    fn test_bed_reader_num_meta_columns_tracking() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1\t0\t100\tm1").unwrap(); // 1 meta
        writeln!(temp_file, "chr1\t100\t200\tm1\tm2\tm3").unwrap(); // 3 meta
        writeln!(temp_file, "chr1\t200\t300\tm1\tm2").unwrap(); // 2 meta
        writeln!(temp_file, "chr1\t300\t400\tm1\tm2\tm3\tm4\tm5").unwrap(); // 5 meta
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();

        // Read first chunk of 2
        let _ = reader.read_chunk(2).unwrap();
        assert_eq!(reader.num_meta_columns(), 3); // Max so far is 3

        // Read remaining
        let _ = reader.read_chunk(10).unwrap();
        assert_eq!(reader.num_meta_columns(), 5); // Max is now 5
    }

    #[test]
    fn test_parse_bed_zero_length_region() {
        // In BED format, start < end typically, but parser doesn't validate
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1\t100\t100\tpoint").unwrap();
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap().unwrap();

        // Point region (start == end)
        assert_eq!(chunk[0].start, 100);
        assert_eq!(chunk[0].end, 100);
        assert_eq!(chunk[0].length(), 1); // 100 - 100 + 1 = 1
    }

    #[test]
    fn test_parse_bed_similar_chromosome_names() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1\t100\t200").unwrap();
        writeln!(temp_file, "chr10\t100\t200").unwrap();
        writeln!(temp_file, "chr1_random\t100\t200").unwrap();
        writeln!(temp_file, "chr1\t300\t400").unwrap();
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();
        let chunk = reader.read_chunk(10).unwrap().unwrap();

        // Should have 4 regions with 3 different chromosomes
        assert_eq!(chunk.len(), 4);
        assert_eq!(chunk[0].chrom, "chr1");
        assert_eq!(chunk[1].chrom, "chr10");
        assert_eq!(chunk[2].chrom, "chr1_random");
        assert_eq!(chunk[3].chrom, "chr1");
    }

    #[test]
    fn test_bed_reader_chunk_size_zero() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "chr1\t100\t200").unwrap();
        temp_file.flush().unwrap();

        let mut reader = rgmatch::parser::BedReader::new(temp_file.path()).unwrap();

        // Chunk size 0 should return None since regions vec is empty
        let chunk = reader.read_chunk(0).unwrap();
        assert!(chunk.is_none());

        // But file position hasn't advanced, so next read should work
        let chunk2 = reader.read_chunk(10).unwrap().unwrap();
        assert_eq!(chunk2.len(), 1);
    }
}

// -------------------------------------------------------------------------
// 11. Parser GTF Module Tests - parse_gtf, extract_attribute edge cases
// -------------------------------------------------------------------------

mod test_parser_gtf {
    use rgmatch::types::Strand;

    // Note: extract_attribute is private, so we test through parsing behavior

    #[test]
    fn test_parse_gtf_multiple_transcripts_per_gene() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"chr1	TEST	gene	1000	5000	.	+	.	gene_id "G1";
chr1	TEST	transcript	1000	2500	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	transcript	2000	5000	.	+	.	gene_id "G1"; transcript_id "T2";
chr1	TEST	exon	1000	1200	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	1500	2500	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	2000	2300	.	+	.	gene_id "G1"; transcript_id "T2";
chr1	TEST	exon	3000	4000	.	+	.	gene_id "G1"; transcript_id "T2";
chr1	TEST	exon	4500	5000	.	+	.	gene_id "G1"; transcript_id "T2";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        let genes = &result.genes_by_chrom["chr1"];
        assert_eq!(genes.len(), 1);

        let gene = &genes[0];
        assert_eq!(gene.gene_id, "G1");
        assert_eq!(gene.transcripts.len(), 2);

        // Verify transcript T1
        let t1 = gene
            .transcripts
            .iter()
            .find(|t| t.transcript_id == "T1")
            .unwrap();
        assert_eq!(t1.exons.len(), 2);
        assert_eq!(t1.start, 1000);
        assert_eq!(t1.end, 2500);

        // Verify transcript T2
        let t2 = gene
            .transcripts
            .iter()
            .find(|t| t.transcript_id == "T2")
            .unwrap();
        assert_eq!(t2.exons.len(), 3);
        assert_eq!(t2.start, 2000);
        assert_eq!(t2.end, 5000);
    }

    #[test]
    fn test_parse_gtf_multiple_genes_same_chromosome() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"chr1	TEST	gene	1000	2000	.	+	.	gene_id "G1";
chr1	TEST	gene	5000	6000	.	-	.	gene_id "G2";
chr1	TEST	gene	10000	15000	.	+	.	gene_id "G3";
chr1	TEST	exon	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	5000	5500	.	-	.	gene_id "G2"; transcript_id "T2";
chr1	TEST	exon	10000	12000	.	+	.	gene_id "G3"; transcript_id "T3";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        let genes = &result.genes_by_chrom["chr1"];
        assert_eq!(genes.len(), 3);

        // Verify genes are present (order preserved from GTF)
        assert_eq!(genes[0].gene_id, "G1");
        assert_eq!(genes[0].strand, Strand::Positive);

        assert_eq!(genes[1].gene_id, "G2");
        assert_eq!(genes[1].strand, Strand::Negative);

        assert_eq!(genes[2].gene_id, "G3");
        assert_eq!(genes[2].strand, Strand::Positive);
    }

    #[test]
    fn test_parse_gtf_exon_only_no_gene_transcript_entries() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"chr1	TEST	exon	1000	1200	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	1500	1800	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	2000	2500	.	+	.	gene_id "G1"; transcript_id "T1";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        let gene = &result.genes_by_chrom["chr1"][0];

        // Gene boundaries should be calculated from exons
        assert_eq!(gene.start, 1000);
        assert_eq!(gene.end, 2500);

        // Transcript boundaries should be calculated from exons
        let transcript = &gene.transcripts[0];
        assert_eq!(transcript.start, 1000);
        assert_eq!(transcript.end, 2500);

        // Exon numbering should still work
        assert_eq!(transcript.exons.len(), 3);
        assert_eq!(transcript.exons[0].exon_number, Some("1".to_string()));
        assert_eq!(transcript.exons[1].exon_number, Some("2".to_string()));
        assert_eq!(transcript.exons[2].exon_number, Some("3".to_string()));
    }

    #[test]
    fn test_parse_gtf_invalid_strand_characters() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"chr1	TEST	gene	1000	2000	.	.	.	gene_id "G1";
chr1	TEST	exon	1000	1200	.	.	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	gene	3000	4000	.	+	.	gene_id "G2";
chr1	TEST	exon	3000	3500	.	+	.	gene_id "G2"; transcript_id "T2";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        // Only G2 should be present (G1 has "." strand which is invalid)
        let genes = &result.genes_by_chrom["chr1"];
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].gene_id, "G2");
    }

    #[test]
    fn test_parse_gtf_comment_lines() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"##description: test GTF file
##provider: TEST
#this is a comment
chr1	TEST	gene	1000	2000	.	+	.	gene_id "G1";
# another comment in the middle
chr1	TEST	exon	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1";
##trailing comment
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        assert!(result.genes_by_chrom.contains_key("chr1"));
        let genes = &result.genes_by_chrom["chr1"];
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].gene_id, "G1");
    }

    #[test]
    fn test_parse_gtf_malformed_lines_wrong_column_count() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"chr1	TEST	gene	1000	2000	.	+	.	gene_id "G1";
chr1	TEST	exon	1000	1200
chr1	only_three_columns
chr1	TEST	exon	1500	2000	.	+	.	gene_id "G1"; transcript_id "T1";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        // Should parse successfully, skipping malformed lines
        let genes = &result.genes_by_chrom["chr1"];
        assert_eq!(genes.len(), 1);

        let transcript = &genes[0].transcripts[0];
        // Only one valid exon line
        assert_eq!(transcript.exons.len(), 1);
        assert_eq!(transcript.exons[0].start, 1500);
    }

    #[test]
    fn test_parse_gtf_empty_lines() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"

chr1	TEST	gene	1000	2000	.	+	.	gene_id "G1";

chr1	TEST	exon	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1";


"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        assert!(result.genes_by_chrom.contains_key("chr1"));
        let genes = &result.genes_by_chrom["chr1"];
        assert_eq!(genes.len(), 1);
    }

    #[test]
    fn test_parse_gtf_gene_boundary_calculation_from_exons() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // No gene entry, only exons from two transcripts
        let gtf_content = r#"chr1	TEST	exon	1000	1200	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	1500	1800	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	800	1000	.	+	.	gene_id "G1"; transcript_id "T2";
chr1	TEST	exon	2000	2500	.	+	.	gene_id "G1"; transcript_id "T2";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        let gene = &result.genes_by_chrom["chr1"][0];

        // Gene should span from min(800) to max(2500)
        assert_eq!(gene.start, 800);
        assert_eq!(gene.end, 2500);
    }

    #[test]
    fn test_parse_gtf_max_lengths_calculation() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"chr1	TEST	gene	1000	2000	.	+	.	gene_id "G1";
chr1	TEST	gene	5000	8000	.	+	.	gene_id "G2";
chr2	TEST	gene	100	500	.	+	.	gene_id "G3";
chr1	TEST	exon	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	5000	6000	.	+	.	gene_id "G2"; transcript_id "T2";
chr2	TEST	exon	100	300	.	+	.	gene_id "G3"; transcript_id "T3";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        // chr1: G1 is 1000 (2000-1000), G2 is 3000 (8000-5000), max is 3000
        assert_eq!(result.max_lengths["chr1"], 3000);

        // chr2: G3 is 400 (500-100)
        assert_eq!(result.max_lengths["chr2"], 400);
    }

    #[test]
    fn test_parse_gtf_empty_file() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = "";

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        assert!(result.genes_by_chrom.is_empty());
        assert!(result.max_lengths.is_empty());
    }

    #[test]
    fn test_parse_gtf_comments_only() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"##description: empty GTF
##provider: TEST
#all comments, no data
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        assert!(result.genes_by_chrom.is_empty());
        assert!(result.max_lengths.is_empty());
    }

    #[test]
    fn test_parse_gtf_multiple_chromosomes() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"chr1	TEST	gene	1000	2000	.	+	.	gene_id "G1";
chr2	TEST	gene	3000	4000	.	-	.	gene_id "G2";
chrX	TEST	gene	5000	6000	.	+	.	gene_id "G3";
chr1	TEST	exon	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1";
chr2	TEST	exon	3000	3500	.	-	.	gene_id "G2"; transcript_id "T2";
chrX	TEST	exon	5000	5500	.	+	.	gene_id "G3"; transcript_id "T3";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        assert_eq!(result.genes_by_chrom.len(), 3);
        assert!(result.genes_by_chrom.contains_key("chr1"));
        assert!(result.genes_by_chrom.contains_key("chr2"));
        assert!(result.genes_by_chrom.contains_key("chrX"));

        assert_eq!(result.genes_by_chrom["chr1"][0].gene_id, "G1");
        assert_eq!(result.genes_by_chrom["chr2"][0].gene_id, "G2");
        assert_eq!(result.genes_by_chrom["chrX"][0].gene_id, "G3");
    }

    #[test]
    fn test_parse_gtf_unsorted_exons() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Exons in random order
        let gtf_content = r#"chr1	TEST	exon	3000	3500	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	1000	1200	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	2000	2500	.	+	.	gene_id "G1"; transcript_id "T1";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        let transcript = &result.genes_by_chrom["chr1"][0].transcripts[0];

        // After renumber_exons, should be sorted by start
        assert_eq!(transcript.exons[0].start, 1000);
        assert_eq!(transcript.exons[0].exon_number, Some("1".to_string()));

        assert_eq!(transcript.exons[1].start, 2000);
        assert_eq!(transcript.exons[1].exon_number, Some("2".to_string()));

        assert_eq!(transcript.exons[2].start, 3000);
        assert_eq!(transcript.exons[2].exon_number, Some("3".to_string()));
    }

    #[test]
    fn test_parse_gtf_other_feature_types_ignored() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let gtf_content = r#"chr1	TEST	gene	1000	2000	.	+	.	gene_id "G1";
chr1	TEST	transcript	1000	2000	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	exon	1000	1500	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	CDS	1100	1400	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	start_codon	1100	1102	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	stop_codon	1398	1400	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	five_prime_utr	1000	1099	.	+	.	gene_id "G1"; transcript_id "T1";
chr1	TEST	three_prime_utr	1401	1500	.	+	.	gene_id "G1"; transcript_id "T1";
"#;

        let mut temp_file = NamedTempFile::with_suffix(".gtf").unwrap();
        temp_file.write_all(gtf_content.as_bytes()).unwrap();
        temp_file.flush().unwrap();

        let result =
            rgmatch::parser::parse_gtf(temp_file.path(), "gene_id", "transcript_id").unwrap();

        let transcript = &result.genes_by_chrom["chr1"][0].transcripts[0];
        // Only exon entries should be counted
        assert_eq!(transcript.exons.len(), 1);
    }
}

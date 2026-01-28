//! Matching logic for genomic regions to gene annotations.

pub mod overlap;
pub mod rules;
pub mod tss;
pub mod tts;

pub use overlap::match_regions_to_genes;
pub use rules::{apply_rules, select_transcript};
pub use tss::check_tss;
pub use tts::check_tts;

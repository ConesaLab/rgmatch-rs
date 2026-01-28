//! Parsers for genomic file formats.

pub mod bed;
pub mod gtf;

pub use bed::parse_bed;
pub use gtf::parse_gtf;

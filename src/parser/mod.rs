//! Molecular structure parsers and loaders.
//!
//! Provides parsing for standard PDB and mmCIF file formats, transparent `.gz`
//! decompression, and RCSB HTTPS structure fetching.

pub mod cif;
pub mod pdb;
pub mod rcsb;

use std::io::Read;
use std::path::Path;

use flate2::read::GzDecoder;

pub use cif::parse_cif;
pub use pdb::parse_pdb;
pub use rcsb::fetch_pdb;

use crate::error::{Result, TermPdbError};
use crate::model::Structure;

/// Loads a macromolecular structure from a file path or a 4-letter RCSB PDB ID.
///
/// Supports:
/// - Local uncompressed files (`.pdb`, `.cif`, `.ent`)
/// - Local gzip-compressed files (`.pdb.gz`, `.cif.gz`, `.ent.gz`, or magic bytes `0x1f, 0x8b`)
/// - 4-character RCSB PDB codes (e.g. `"1CRN"`, `"1ubq"`, `"7V67"`)
pub fn load_structure(source: &str) -> Result<Structure> {
    let trimmed = source.trim();
    if trimmed.is_empty() {
        return Err(TermPdbError::InvalidStructure(
            "Empty structure source provided".to_string(),
        ));
    }

    let path = Path::new(trimmed);
    if path.exists() && path.is_file() {
        let bytes = std::fs::read(path)?;
        let content = if bytes.starts_with(&[0x1f, 0x8b])
            || trimmed.ends_with(".gz")
            || trimmed.ends_with(".GZ")
        {
            let mut decoder = GzDecoder::new(&bytes[..]);
            let mut s = String::new();
            decoder.read_to_string(&mut s)?;
            s
        } else {
            String::from_utf8(bytes).map_err(|e| TermPdbError::ParseError(e.to_string()))?
        };

        if is_cif_format(trimmed, &content) {
            parse_cif(&content)
        } else {
            parse_pdb(&content)
        }
    } else if is_potential_pdb_id(trimmed) {
        let data = fetch_pdb(trimmed)?;
        if is_cif_format(trimmed, &data) {
            parse_cif(&data)
        } else {
            parse_pdb(&data)
        }
    } else {
        Err(TermPdbError::Io(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!(
                "File not found or invalid structure identifier: '{}'",
                source
            ),
        )))
    }
}

fn is_potential_pdb_id(s: &str) -> bool {
    let trimmed = s.trim();
    (trimmed.len() == 4
        || (trimmed.len() >= 4 && trimmed.len() <= 12 && trimmed.starts_with("pdb_")))
        && trimmed
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_')
}

fn is_cif_format(path_or_id: &str, content: &str) -> bool {
    let lower_path = path_or_id.to_ascii_lowercase();
    if lower_path.ends_with(".cif") || lower_path.ends_with(".cif.gz") {
        return true;
    }
    if lower_path.ends_with(".pdb")
        || lower_path.ends_with(".pdb.gz")
        || lower_path.ends_with(".ent")
    {
        return false;
    }

    // Check content signatures
    content.contains("_atom_site.")
        || content.contains("data_")
        || content.contains("_struct.")
        || content.contains("_entry.id")
}

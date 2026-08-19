//! RCSB PDB online structure fetcher.
//!
//! Fetches structure files directly from RCSB PDB over HTTPS with transparent
//! Gzip decompression and format fallbacks (.pdb.gz -> .cif.gz -> .pdb -> .cif).

use std::io::Read;
use std::time::Duration;

use flate2::read::GzDecoder;

use crate::error::{Result, TermPdbError};

/// Fetches PDB or mmCIF structure content from RCSB for the given ID.
///
/// Downloads compressed (.gz) or uncompressed files from `https://files.rcsb.org/download/{id}`.
pub fn fetch_pdb(pdb_id: &str) -> Result<String> {
    let trimmed = pdb_id.trim();
    if trimmed.is_empty()
        || !trimmed
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
    {
        return Err(TermPdbError::InvalidStructure(format!(
            "Invalid structure ID: '{}'",
            pdb_id
        )));
    }

    let id_upper = trimmed.to_ascii_uppercase();
    let id_lower = trimmed.to_ascii_lowercase();

    let urls = [
        format!("https://files.rcsb.org/download/{}.pdb.gz", id_upper),
        format!("https://files.rcsb.org/download/{}.cif.gz", id_upper),
        format!("https://files.rcsb.org/download/{}.pdb", id_upper),
        format!("https://files.rcsb.org/download/{}.cif", id_upper),
        format!("https://files.rcsb.org/download/{}.pdb.gz", id_lower),
        format!("https://files.rcsb.org/download/{}.cif.gz", id_lower),
    ];

    let mut last_err = None;

    for url in &urls {
        match try_fetch_url(url) {
            Ok(text) => return Ok(text),
            Err(e) => {
                last_err = Some(e);
            }
        }
    }

    Err(last_err.unwrap_or_else(|| {
        TermPdbError::NetworkError(format!("Failed to fetch structure '{}' from RCSB", pdb_id))
    }))
}

fn try_fetch_url(url: &str) -> Result<String> {
    let resp = ureq::get(url)
        .timeout(Duration::from_secs(10))
        .call()
        .map_err(|e| TermPdbError::NetworkError(format!("HTTP error for {}: {}", url, e)))?;

    let mut reader = resp.into_reader();
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes).map_err(TermPdbError::Io)?;

    if bytes.starts_with(&[0x1f, 0x8b]) || url.ends_with(".gz") {
        let mut decoder = GzDecoder::new(&bytes[..]);
        let mut decompressed = String::new();
        decoder
            .read_to_string(&mut decompressed)
            .map_err(TermPdbError::Io)?;
        Ok(decompressed)
    } else {
        String::from_utf8(bytes).map_err(|e| TermPdbError::ParseError(e.to_string()))
    }
}

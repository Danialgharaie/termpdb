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
/// Downloads compressed (.gz) or uncompressed files from `https://files.rcsb.org/download/{id}`,
/// trying the default endpoint list in order until one succeeds.
pub fn fetch_pdb(pdb_id: &str) -> Result<String> {
    let trimmed = pdb_id.trim();
    let urls = default_download_urls(trimmed);
    fetch_pdb_from(trimmed, &urls)
}

/// Default RCSB download endpoints, tried in order: legacy PDB (compressed,
/// then raw), mmCIF (compressed, then raw), and finally lower-case ID
/// spellings of the compressed variants.
fn default_download_urls(pdb_id: &str) -> Vec<String> {
    let id_upper = pdb_id.to_ascii_uppercase();
    let id_lower = pdb_id.to_ascii_lowercase();
    vec![
        format!("https://files.rcsb.org/download/{id_upper}.pdb.gz"),
        format!("https://files.rcsb.org/download/{id_upper}.cif.gz"),
        format!("https://files.rcsb.org/download/{id_upper}.pdb"),
        format!("https://files.rcsb.org/download/{id_upper}.cif"),
        format!("https://files.rcsb.org/download/{id_lower}.pdb.gz"),
        format!("https://files.rcsb.org/download/{id_lower}.cif.gz"),
    ]
}

/// Downloads structure text by trying each URL in order; the first success wins.
///
/// This is the injectable seam behind [`fetch_pdb`]: callers may supply their
/// own complete endpoint URLs (alternative mirrors, or a local test server) so
/// fallback behavior can be exercised offline. The ID is validated before any
/// request is made, and every candidate URL is tried verbatim.
pub fn fetch_pdb_from(id: &str, urls: &[String]) -> Result<String> {
    let trimmed = id.trim();
    if trimmed.is_empty()
        || !trimmed
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
    {
        return Err(TermPdbError::InvalidStructure(format!(
            "Invalid structure ID: '{}'",
            id
        )));
    }

    let mut last_err = None;

    for url in urls {
        match try_fetch_url(url) {
            Ok(text) => return Ok(text),
            Err(e) => {
                last_err = Some(e);
            }
        }
    }

    Err(last_err.unwrap_or_else(|| {
        TermPdbError::NetworkError(format!(
            "Failed to fetch structure '{}' from {} source(s)",
            trimmed,
            urls.len()
        ))
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

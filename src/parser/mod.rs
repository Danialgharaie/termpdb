//! Molecular structure parsers and loaders.
//!
//! Provides parsing for standard PDB and mmCIF file formats, transparent `.gz`
//! decompression, and RCSB HTTPS structure fetching.

pub mod cif;
pub mod pdb;
pub mod rcsb;

use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

use flate2::read::GzDecoder;

pub use cif::parse_cif;
pub use pdb::parse_pdb;
pub use rcsb::fetch_pdb;

use crate::error::{Result, TermPdbError};
use crate::model::{Bond, BondDetector, Chain, Model, Residue, SecondaryStructure, Structure};

/// Per-model parse accumulator used by PDB and mmCIF loaders.
#[derive(Default)]
pub(crate) struct ModelAccum {
    pub atoms: Vec<crate::model::Atom>,
    pub chain_order: Vec<String>,
    pub chain_residues: HashMap<String, Vec<Residue>>,
    pub serial_to_idx: HashMap<i32, usize>,
}

/// Builds a [`Model`] from parsed atoms/residues and secondary-structure ranges.
pub(crate) fn assemble_model(
    serial: i32,
    accum: ModelAccum,
    helices: &[(String, i32, i32)],
    sheets: &[(String, i32, i32)],
) -> Model {
    let ModelAccum {
        atoms,
        chain_order,
        mut chain_residues,
        ..
    } = accum;

    let mut model = Model::new(serial);
    model.atoms = atoms;

    for chain_id in chain_order {
        let mut chain = Chain::new(&chain_id);
        if let Some(mut residues) = chain_residues.remove(&chain_id) {
            for res in &mut residues {
                for (h_chain, h_init, h_end) in helices {
                    if h_chain == &chain_id && res.seq >= *h_init && res.seq <= *h_end {
                        res.secondary_structure = SecondaryStructure::Helix;
                        break;
                    }
                }
                if res.secondary_structure == SecondaryStructure::Coil {
                    for (s_chain, s_init, s_end) in sheets {
                        if s_chain == &chain_id && res.seq >= *s_init && res.seq <= *s_end {
                            res.secondary_structure = SecondaryStructure::Sheet;
                            break;
                        }
                    }
                }
            }
            chain.residues = residues;
        }
        model.chains.push(chain);
    }

    model.bonds = BondDetector::detect_bonds(&model.atoms);
    model
}

/// Adds CONECT pairs to a model using that model's atom-serial map. Never crosses models.
pub(crate) fn apply_conect(
    model: &mut Model,
    conect_pairs: &[(i32, i32)],
    serial_to_idx: &HashMap<i32, usize>,
) {
    for &(src_serial, dst_serial) in conect_pairs {
        let Some(&idx1) = serial_to_idx.get(&src_serial) else {
            continue;
        };
        let Some(&idx2) = serial_to_idx.get(&dst_serial) else {
            continue;
        };
        if idx1 == idx2 {
            continue;
        }
        let (a, b) = if idx1 < idx2 {
            (idx1, idx2)
        } else {
            (idx2, idx1)
        };
        let exists = model.bonds.iter().any(|bond| {
            let (ba, bb) = if bond.atom1_idx < bond.atom2_idx {
                (bond.atom1_idx, bond.atom2_idx)
            } else {
                (bond.atom2_idx, bond.atom1_idx)
            };
            ba == a && bb == b
        });
        if !exists {
            model.bonds.push(Bond::single(a, b));
        }
    }
}

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

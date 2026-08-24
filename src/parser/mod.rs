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
use crate::model::{
    Bond, BondDetector, Chain, Element, Model, Residue, SecondaryStructure, Structure,
    element_by_symbol,
};

/// Per-model parse accumulator used by PDB and mmCIF loaders.
#[derive(Default)]
pub(crate) struct ModelAccum {
    pub atoms: Vec<crate::model::Atom>,
    pub chain_order: Vec<String>,
    pub chain_residues: HashMap<String, Vec<Residue>>,
    pub serial_to_idx: HashMap<i32, usize>,
}

/// Ranks alternate-location conformers competing for one atom site.
///
/// Lower wins: the blank conformer (`alt_loc == None`, PDB column 17 blank /
/// mmCIF `label_alt_id` `.` or `?`) is the default deposition choice, `A` is
/// the conventional primary alternate, and any other identifier comes last.
fn alt_loc_rank(alt_loc: Option<char>) -> u8 {
    match alt_loc {
        None => 0,
        Some('A') => 1,
        Some(_) => 2,
    }
}

/// Atom-site key used to group alternate-location conformers:
/// `(chain_id, res_seq, insertion code, atom name)`.
type SiteKey = (String, i32, Option<char>, String);

/// Keeps exactly one alternate-location conformer per atom site, dropping the rest.
///
/// Policy: at most one conformer survives per site key
/// `(chain_id, res_seq, ins_code, atom name)` — preference order is the blank
/// conformer, then `'A'`, then the first-encountered alternate. Atoms recorded
/// without an altloc always win their site, so structures with no alternates
/// are untouched. Dropped atoms are removed from `atoms` (with indices
/// compacted and `Atom::index` renumbered), from every residue's
/// `atom_indices`, and from `serial_to_idx`, so bond detection, rendering, and
/// atom/residue counts never see duplicate conformers. `Atom::alt_loc` is kept
/// as-is for the surviving atoms.
///
/// Applied identically to the PDB and mmCIF loaders via [`assemble_model`].
pub(crate) fn keep_one_alt_loc_conformer(accum: &mut ModelAccum) {
    // Fast path: files without any altloc records need no rework at all.
    if !accum.atoms.iter().any(|atom| atom.alt_loc.is_some()) {
        return;
    }

    // Atom carries no insertion code, so recover the per-atom site insertion
    // code from the residue grouping built while parsing.
    let mut ins_codes = vec![None; accum.atoms.len()];
    for residues in accum.chain_residues.values() {
        for res in residues {
            for &idx in &res.atom_indices {
                if let Some(slot) = ins_codes.get_mut(idx) {
                    *slot = res.ins_code;
                }
            }
        }
    }

    // Winning atom position per site: lowest rank, ties (e.g. two `B`
    // records) broken by first occurrence in file order.
    let mut winners: HashMap<SiteKey, (u8, usize)> = HashMap::new();
    for (pos, atom) in accum.atoms.iter().enumerate() {
        let rank = alt_loc_rank(atom.alt_loc);
        let key = (
            atom.chain_id.clone(),
            atom.res_seq,
            ins_codes[pos],
            atom.name.clone(),
        );
        match winners.get(&key) {
            Some(&(best_rank, _)) if best_rank <= rank => {}
            _ => {
                winners.insert(key, (rank, pos));
            }
        }
    }

    let atoms = std::mem::take(&mut accum.atoms);
    let mut keep = vec![false; atoms.len()];
    for &(_, pos) in winners.values() {
        keep[pos] = true;
    }

    let mut remap: HashMap<usize, usize> = HashMap::new();
    for (pos, mut atom) in atoms.into_iter().enumerate() {
        if !keep[pos] {
            continue;
        }
        let new_idx = accum.atoms.len();
        remap.insert(pos, new_idx);
        atom.index = new_idx;
        accum.atoms.push(atom);
    }

    // Compact residue -> atom references; dropped conformers fall out here.
    for residues in accum.chain_residues.values_mut() {
        for res in residues {
            res.atom_indices.retain_mut(|idx| match remap.get(idx) {
                Some(&new_idx) => {
                    *idx = new_idx;
                    true
                }
                None => false,
            });
        }
    }

    // Rebuild the serial map from the surviving atoms so CONECT resolution
    // skips dropped serials and points at compacted indices. Iterating in file
    // order preserves the parser's last-insert-wins behaviour for duplicates.
    accum.serial_to_idx = accum
        .atoms
        .iter()
        .map(|atom| (atom.serial, atom.index))
        .collect();
}

/// Builds a [`Model`] from parsed atoms/residues and secondary-structure ranges.
pub(crate) fn assemble_model(
    serial: i32,
    mut accum: ModelAccum,
    helices: &[(String, i32, i32)],
    sheets: &[(String, i32, i32)],
) -> Model {
    // Collapse alternate-location conformers before anything (bond detection,
    // rendering, counts) can observe duplicate atoms.
    keep_one_alt_loc_conformer(&mut accum);

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

/// Infers an element for an atom whose record carries no usable element symbol.
///
/// Shared by the PDB and mmCIF parsers. Lookup order matters for scientific
/// correctness:
///
/// 1. The periodic table is consulted for the leading two characters first, so
///    heteroatom/ion names such as `HG`, `CD`, or `SE` resolve to mercury,
///    cadmium, and selenium instead of being misread by naming conventions as
///    hydrogen/carbon/sulfur.
/// 2. Standard amino-acid and nucleotide residues are exempt from step 1: their
///    atom names put the element in the first character followed by a locator
///    letter or digit (`CD1`, `NE2`, `OG`, `SD`, `HG1`), and several locators
///    coincide with valid element symbols (Nd, Ne, Ce, He, Hg, Og), so there a
///    naming-convention hit takes precedence.
/// 3. A single-character lookup (`N`, `C`, `O`, `S`, ...) applies last.
pub(crate) fn infer_element(atom_name: &str, res_name: &str, is_hetatm: bool) -> Element {
    let trimmed = atom_name.trim();
    if trimmed.is_empty() {
        return Element::unknown();
    }

    let cleaned = trimmed.trim_start_matches(|c: char| c.is_ascii_digit());
    if cleaned.is_empty() {
        return Element::unknown();
    }

    let res_trimmed = res_name.trim();

    // Special case for protein C-alpha: in non-HETATM or amino acid residues, "CA" is Carbon
    if cleaned.eq_ignore_ascii_case("CA")
        && !is_hetatm
        && !res_trimmed.eq_ignore_ascii_case("CA")
        && !res_trimmed.eq_ignore_ascii_case("CAL")
    {
        return element_by_symbol("C");
    }

    // Index by chars, not bytes: atom names may contain non-ASCII characters
    // and byte slicing would panic on non-char boundaries.
    let mut chars = cleaned.chars();
    let first_char = chars.next();
    let second_char = chars.next().map(|c| c.to_ascii_uppercase());

    if let (Some(first_char), Some(second_char)) = (first_char, second_char) {
        // Periodic table before naming conventions (see doc comment).
        let elem2_end = cleaned
            .char_indices()
            .nth(2)
            .map(|(i, _)| i)
            .unwrap_or(cleaned.len());
        let elem2 = element_by_symbol(&cleaned[..elem2_end]);

        let convention_elem = convention_element(first_char, second_char);
        let symbol_lookup_first =
            convention_elem.is_none() || !is_standard_polymer_residue(res_trimmed);

        if symbol_lookup_first && elem2.atomic_number != 0 {
            return elem2;
        }
        if let Some(elem) = convention_elem {
            return elem;
        }
        if elem2.atomic_number != 0 {
            return elem2;
        }
    }

    // Check 1-letter element
    let elem1_end = cleaned
        .char_indices()
        .nth(1)
        .map(|(i, _)| i)
        .unwrap_or(cleaned.len());
    let candidate1 = &cleaned[..elem1_end];
    let elem1 = element_by_symbol(candidate1);
    if elem1.atomic_number != 0 {
        return elem1;
    }

    Element::unknown()
}

/// Standard IUPAC amino acid / nucleotide atom prefix conventions: returns the
/// element implied by treating `second_char` as a residue-atom locator letter
/// behind the element `first_char` (e.g. `HB2` → hydrogen, `SD` → sulfur).
fn convention_element(first_char: char, second_char: char) -> Option<Element> {
    let by_first = |c: char| first_char == c || first_char == c.to_ascii_lowercase();
    let matches_second = |candidates: &[char]| candidates.contains(&second_char);

    if by_first('H') && matches_second(&['B', 'G', 'D', 'E', 'Z', 'H', 'Q', '1', '2', '3']) {
        return Some(element_by_symbol("H"));
    }
    if by_first('C')
        && matches_second(&['A', 'B', 'G', 'D', 'E', 'Z', 'H', '1', '2', '3', '\'', '*'])
    {
        return Some(element_by_symbol("C"));
    }
    if by_first('N') && matches_second(&['E', 'D', 'H', 'Z', '1', '2', '3']) {
        return Some(element_by_symbol("N"));
    }
    if by_first('O') && matches_second(&['G', 'D', 'E', 'H', 'X', 'P', '1', '2', '3', '\'', '*']) {
        return Some(element_by_symbol("O"));
    }
    if by_first('S') && matches_second(&['D', 'G', 'E', '1', '2', '3']) {
        return Some(element_by_symbol("S"));
    }

    None
}

/// Whether `res_name` is a standard amino acid or nucleotide residue name whose
/// atoms follow the locator naming conventions handled by [`convention_element`].
fn is_standard_polymer_residue(res_name: &str) -> bool {
    matches!(
        res_name.to_ascii_uppercase().as_str(),
        "ALA"
            | "ARG"
            | "ASN"
            | "ASP"
            | "CYS"
            | "GLN"
            | "GLU"
            | "GLY"
            | "HIS"
            | "ILE"
            | "LEU"
            | "LYS"
            | "MET"
            | "PHE"
            | "PRO"
            | "SER"
            | "THR"
            | "TRP"
            | "TYR"
            | "VAL"
            | "DA"
            | "DC"
            | "DG"
            | "DT"
            | "DU"
            | "A"
            | "C"
            | "G"
            | "U"
    )
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

//! PDB format parser.
//!
//! Parses standard Protein Data Bank (PDB) formatted text into a [`Structure`].
//! Supports HEADER, TITLE, HELIX, SHEET, ATOM, HETATM, and CONECT records.

use std::collections::HashMap;

use crate::error::{Result, TermPdbError};
use crate::math::Vec3;
use crate::model::{
    Atom, Bond, Chain, Element, Residue, SecondaryStructure, Structure, element_by_symbol,
};

/// Parses a PDB-formatted string into a [`Structure`].
pub fn parse_pdb(input: &str) -> Result<Structure> {
    let mut structure = Structure::new("");
    let mut title_lines: Vec<String> = Vec::new();
    let mut serial_to_idx: HashMap<i32, usize> = HashMap::new();
    let mut conect_pairs: Vec<(i32, i32)> = Vec::new();

    let mut helices: Vec<(String, i32, i32)> = Vec::new();
    let mut sheets: Vec<(String, i32, i32)> = Vec::new();

    // Map: chain_id -> Vec<Residue>
    // We maintain ordered list of chain IDs and residues for each chain.
    let mut chain_order: Vec<String> = Vec::new();
    let mut chain_residues: HashMap<String, Vec<Residue>> = HashMap::new();

    for line in input.lines() {
        if line.is_empty() {
            continue;
        }

        let record_type = safe_slice(line, 0, 6).trim();

        match record_type {
            "HEADER" => {
                if structure.id_code.is_none() {
                    let id_code = safe_slice(line, 62, 66).trim();
                    if !id_code.is_empty() {
                        structure.id_code = Some(id_code.to_string());
                    }
                }
                let classification = safe_slice(line, 10, 50).trim();
                if !classification.is_empty() {
                    structure
                        .metadata
                        .insert("classification".to_string(), classification.to_string());
                }
                let dep_date = safe_slice(line, 50, 59).trim();
                if !dep_date.is_empty() {
                    structure
                        .metadata
                        .insert("deposition_date".to_string(), dep_date.to_string());
                }
            }
            "TITLE" => {
                let text = safe_slice(line, 10, 80).trim();
                if !text.is_empty() {
                    title_lines.push(text.to_string());
                }
            }
            "HELIX" => {
                let mut chain_id = safe_slice(line, 19, 20).trim().to_string();
                let mut init_seq = safe_slice(line, 21, 25).trim().parse::<i32>().ok();
                let mut end_seq = safe_slice(line, 33, 37).trim().parse::<i32>().ok();

                if init_seq.is_none() || end_seq.is_none() {
                    if end_seq.is_none() {
                        end_seq = safe_slice(line, 32, 39)
                            .split_whitespace()
                            .next()
                            .and_then(|s| s.parse::<i32>().ok());
                    }
                    let tokens: Vec<&str> = line.split_whitespace().collect();
                    if tokens.len() >= 9 {
                        if chain_id.is_empty() {
                            chain_id = tokens[4].to_string();
                        }
                        if init_seq.is_none() {
                            init_seq = tokens[5].parse::<i32>().ok();
                        }
                        if end_seq.is_none() {
                            end_seq = tokens[8].parse::<i32>().ok();
                        }
                    }
                }

                if let (Some(init), Some(end)) = (init_seq, end_seq) {
                    let cid = if chain_id.is_empty() {
                        "A".to_string()
                    } else {
                        chain_id
                    };
                    helices.push((cid, init, end));
                }
            }
            "SHEET" => {
                let mut chain_id = safe_slice(line, 21, 22).trim().to_string();
                let mut init_seq = safe_slice(line, 22, 26).trim().parse::<i32>().ok();
                let mut end_seq = safe_slice(line, 33, 37).trim().parse::<i32>().ok();

                if init_seq.is_none() || end_seq.is_none() {
                    if end_seq.is_none() {
                        end_seq = safe_slice(line, 32, 39)
                            .split_whitespace()
                            .next()
                            .and_then(|s| s.parse::<i32>().ok());
                    }
                    let tokens: Vec<&str> = line.split_whitespace().collect();
                    if tokens.len() >= 10 {
                        if chain_id.is_empty() {
                            chain_id = tokens[5].to_string();
                        }
                        if init_seq.is_none() {
                            init_seq = tokens[6].parse::<i32>().ok();
                        }
                        if end_seq.is_none() {
                            end_seq = tokens[9].parse::<i32>().ok();
                        }
                    }
                }

                if let (Some(init), Some(end)) = (init_seq, end_seq) {
                    let cid = if chain_id.is_empty() {
                        "A".to_string()
                    } else {
                        chain_id
                    };
                    sheets.push((cid, init, end));
                }
            }
            "ATOM" | "HETATM" => {
                let is_hetatm = record_type == "HETATM";
                let serial = safe_slice(line, 6, 11)
                    .trim()
                    .parse::<i32>()
                    .unwrap_or(structure.atoms.len() as i32 + 1);
                let atom_name = safe_slice(line, 12, 16).trim();
                let alt_loc_char = safe_slice(line, 16, 17).chars().next();
                let alt_loc = alt_loc_char.filter(|&c| c != ' ');

                let res_name = safe_slice(line, 17, 20).trim();
                let chain_id_str = safe_slice(line, 21, 22).trim();
                let chain_id = if chain_id_str.is_empty() {
                    "A".to_string()
                } else {
                    chain_id_str.to_string()
                };

                let res_seq = safe_slice(line, 22, 26).trim().parse::<i32>().unwrap_or(1);
                let ins_code_char = safe_slice(line, 26, 27).chars().next();
                let ins_code = ins_code_char.filter(|&c| c != ' ');

                let x = safe_slice(line, 30, 38)
                    .trim()
                    .parse::<f32>()
                    .map_err(|_| {
                        TermPdbError::ParseError(format!("Invalid X coordinate on line: {}", line))
                    })?;
                let y = safe_slice(line, 38, 46)
                    .trim()
                    .parse::<f32>()
                    .map_err(|_| {
                        TermPdbError::ParseError(format!("Invalid Y coordinate on line: {}", line))
                    })?;
                let z = safe_slice(line, 46, 54)
                    .trim()
                    .parse::<f32>()
                    .map_err(|_| {
                        TermPdbError::ParseError(format!("Invalid Z coordinate on line: {}", line))
                    })?;

                let occupancy = safe_slice(line, 54, 60)
                    .trim()
                    .parse::<f32>()
                    .unwrap_or(1.0);
                let b_factor = safe_slice(line, 60, 66)
                    .trim()
                    .parse::<f32>()
                    .unwrap_or(0.0);

                let elem_str = safe_slice(line, 76, 78).trim();
                let element = if !elem_str.is_empty() {
                    let elem = element_by_symbol(elem_str);
                    if elem.atomic_number != 0 {
                        elem
                    } else {
                        infer_element(atom_name, res_name, is_hetatm)
                    }
                } else {
                    infer_element(atom_name, res_name, is_hetatm)
                };

                let charge_str = safe_slice(line, 78, 80).trim();
                let charge = parse_charge(charge_str);

                let atom_idx = structure.atoms.len();
                let mut atom = Atom::new(
                    atom_idx,
                    serial,
                    atom_name,
                    element,
                    Vec3::new(x, y, z),
                    b_factor,
                    res_name,
                    res_seq,
                    &chain_id,
                    is_hetatm,
                );
                atom.occupancy = occupancy;
                atom.alt_loc = alt_loc;
                atom.charge = charge;

                structure.atoms.push(atom);
                serial_to_idx.insert(serial, atom_idx);

                // Add to chain/residue structure
                if !chain_order.contains(&chain_id) {
                    chain_order.push(chain_id.clone());
                    chain_residues.insert(chain_id.clone(), Vec::new());
                }

                let residues = chain_residues.get_mut(&chain_id).unwrap();
                let is_same_as_last = if let Some(last_res) = residues.last() {
                    last_res.seq == res_seq
                        && last_res.ins_code == ins_code
                        && last_res.name == res_name
                } else {
                    false
                };

                if is_same_as_last {
                    residues.last_mut().unwrap().atom_indices.push(atom_idx);
                } else {
                    let mut new_res = Residue::new(res_seq, res_name, &chain_id);
                    new_res.ins_code = ins_code;
                    new_res.atom_indices.push(atom_idx);
                    residues.push(new_res);
                }
            }
            "CONECT" => {
                let src_serial = safe_slice(line, 6, 11).trim().parse::<i32>().ok();
                if let Some(src) = src_serial {
                    let offsets = [11, 16, 21, 26, 31, 36, 41, 46];
                    for &offset in &offsets {
                        let dst_slice = safe_slice(line, offset, offset + 5).trim();
                        if let Ok(dst) = dst_slice.parse::<i32>() {
                            conect_pairs.push((src, dst));
                        }
                    }
                }
            }
            _ => {}
        }
    }

    if !title_lines.is_empty() {
        structure.title = title_lines.join(" ");
    }

    // Build chains
    for chain_id in chain_order {
        let mut chain = Chain::new(&chain_id);
        if let Some(mut residues) = chain_residues.remove(&chain_id) {
            // Apply secondary structure
            for res in &mut residues {
                for (h_chain, h_init, h_end) in &helices {
                    if h_chain == &chain_id && res.seq >= *h_init && res.seq <= *h_end {
                        res.secondary_structure = SecondaryStructure::Helix;
                        break;
                    }
                }
                if res.secondary_structure == SecondaryStructure::Coil {
                    for (s_chain, s_init, s_end) in &sheets {
                        if s_chain == &chain_id && res.seq >= *s_init && res.seq <= *s_end {
                            res.secondary_structure = SecondaryStructure::Sheet;
                            break;
                        }
                    }
                }
            }
            chain.residues = residues;
        }
        structure.add_chain(chain);
    }

    // Auto-detect covalent bonds
    structure.build_bonds();

    // Incorporate explicit CONECT bonds
    for (src_serial, dst_serial) in conect_pairs {
        if let (Some(&idx1), Some(&idx2)) = (
            serial_to_idx.get(&src_serial),
            serial_to_idx.get(&dst_serial),
        ) {
            if idx1 == idx2 {
                continue;
            }
            let (a, b) = if idx1 < idx2 {
                (idx1, idx2)
            } else {
                (idx2, idx1)
            };
            let exists = structure.bonds.iter().any(|bond| {
                let (ba, bb) = if bond.atom1_idx < bond.atom2_idx {
                    (bond.atom1_idx, bond.atom2_idx)
                } else {
                    (bond.atom2_idx, bond.atom1_idx)
                };
                ba == a && bb == b
            });
            if !exists {
                structure.add_bond(Bond::single(a, b));
            }
        }
    }

    Ok(structure)
}

fn safe_slice(s: &str, start: usize, end: usize) -> &str {
    if start >= s.len() {
        return "";
    }
    let end = end.min(s.len());
    &s[start..end]
}

fn parse_charge(charge_str: &str) -> Option<i8> {
    let trimmed = charge_str.trim();
    if trimmed.is_empty() {
        return None;
    }

    if let Some(stripped) = trimmed.strip_suffix('+') {
        let mag = stripped.parse::<i8>().unwrap_or(1);
        Some(mag)
    } else if let Some(stripped) = trimmed.strip_suffix('-') {
        let mag = stripped.parse::<i8>().unwrap_or(1);
        Some(-mag)
    } else if let Some(stripped) = trimmed.strip_prefix('+') {
        let mag = stripped.parse::<i8>().unwrap_or(1);
        Some(mag)
    } else if let Some(stripped) = trimmed.strip_prefix('-') {
        let mag = stripped.parse::<i8>().unwrap_or(1);
        Some(-mag)
    } else {
        trimmed.parse::<i8>().ok()
    }
}

fn infer_element(atom_name: &str, res_name: &str, is_hetatm: bool) -> Element {
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

    // Check 2-letter element matches if cleaned >= 2
    if cleaned.len() >= 2 {
        let first_char = cleaned.chars().next().unwrap();
        let second_char = cleaned.chars().nth(1).unwrap().to_ascii_uppercase();

        // Check standard IUPAC amino acid / nucleotide atom prefix conventions
        if (first_char == 'H' || first_char == 'h')
            && matches!(
                second_char,
                'B' | 'G' | 'D' | 'E' | 'Z' | 'H' | 'Q' | '1' | '2' | '3'
            )
        {
            return element_by_symbol("H");
        }
        if (first_char == 'C' || first_char == 'c')
            && matches!(
                second_char,
                'A' | 'B' | 'G' | 'D' | 'E' | 'Z' | 'H' | '1' | '2' | '3' | '\'' | '*'
            )
        {
            return element_by_symbol("C");
        }
        if (first_char == 'N' || first_char == 'n')
            && matches!(second_char, 'E' | 'D' | 'H' | 'Z' | '1' | '2' | '3')
        {
            return element_by_symbol("N");
        }
        if (first_char == 'O' || first_char == 'o')
            && matches!(
                second_char,
                'G' | 'D' | 'E' | 'H' | 'X' | 'P' | '1' | '2' | '3' | '\'' | '*'
            )
        {
            return element_by_symbol("O");
        }
        if (first_char == 'S' || first_char == 's')
            && matches!(second_char, 'D' | 'G' | 'E' | '1' | '2' | '3')
        {
            return element_by_symbol("S");
        }

        let candidate2 = &cleaned[..2];
        let elem2 = element_by_symbol(candidate2);
        if elem2.atomic_number != 0 {
            return elem2;
        }
    }

    // Check 1-letter element
    let candidate1 = &cleaned[..1];
    let elem1 = element_by_symbol(candidate1);
    if elem1.atomic_number != 0 {
        return elem1;
    }

    Element::unknown()
}

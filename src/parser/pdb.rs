//! PDB format parser.
//!
//! Parses standard Protein Data Bank (PDB) formatted text into a [`Structure`].
//! Supports HEADER, TITLE, HELIX, SHEET, ATOM, HETATM, and CONECT records.

use std::collections::{BTreeMap, HashMap};

use crate::error::{Result, TermPdbError};
use crate::math::Vec3;
use crate::model::assembly::affine_from_rows;
use crate::model::{Assembly, AssemblyGen, Atom, Residue, Structure, element_by_symbol};
use crate::parser::{ModelAccum, apply_conect, assemble_model, infer_element};

/// Parses a PDB-formatted string into a [`Structure`].
pub fn parse_pdb(input: &str) -> Result<Structure> {
    let mut structure = Structure::new("");
    let mut title_lines: Vec<String> = Vec::new();
    let mut conect_pairs: Vec<(i32, i32)> = Vec::new();

    let mut helices: Vec<(String, i32, i32)> = Vec::new();
    let mut sheets: Vec<(String, i32, i32)> = Vec::new();

    let mut accums: BTreeMap<i32, ModelAccum> = BTreeMap::new();
    let mut current_serial: i32 = 1;
    let mut remark350 = Remark350Parser::default();

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
            "MODEL" => {
                let serial = safe_slice(line, 10, 14)
                    .trim()
                    .parse::<i32>()
                    .ok()
                    .or_else(|| {
                        line.split_whitespace()
                            .nth(1)
                            .and_then(|s| s.parse::<i32>().ok())
                    })
                    .unwrap_or(current_serial);
                current_serial = serial;
                accums.entry(current_serial).or_default();
            }
            "ENDMDL" => {}
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
                let accum = accums.entry(current_serial).or_default();
                let serial = safe_slice(line, 6, 11)
                    .trim()
                    .parse::<i32>()
                    .unwrap_or(accum.atoms.len() as i32 + 1);
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

                let atom_idx = accum.atoms.len();
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

                accum.atoms.push(atom);
                accum.serial_to_idx.insert(serial, atom_idx);

                if !accum.chain_order.contains(&chain_id) {
                    accum.chain_order.push(chain_id.clone());
                    accum.chain_residues.insert(chain_id.clone(), Vec::new());
                }

                let residues = accum.chain_residues.get_mut(&chain_id).unwrap();
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
            "REMARK" => {
                if safe_slice(line, 7, 10).trim() == "350" {
                    remark350.feed(line);
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

    let has_ss = !helices.is_empty() || !sheets.is_empty();
    let mut models = Vec::with_capacity(accums.len());
    for (serial, accum) in accums {
        let mut model = assemble_model(serial, accum, &helices, &sheets);
        // Derived after assembly: the altloc policy may have dropped
        // duplicate conformers, shifting atom indices.
        let serial_to_idx: HashMap<i32, usize> = model
            .atoms
            .iter()
            .map(|atom| (atom.serial, atom.index))
            .collect();
        apply_conect(&mut model, &conect_pairs, &serial_to_idx);
        models.push(model);
    }
    structure.set_models(models);
    structure.set_assemblies(remark350.finish());

    if !has_ss {
        crate::model::dssp::assign_dssp(&mut structure);
    }

    Ok(structure)
}

#[derive(Default)]
struct Remark350Parser {
    assemblies: Vec<Assembly>,
    current: Option<Assembly>,
    pending_chains: Vec<String>,
    biomt_rows: HashMap<String, [[f32; 4]; 3]>,
}

impl Remark350Parser {
    fn feed(&mut self, line: &str) {
        let rest = safe_slice(line, 11, line.len()).trim();
        if rest.is_empty() {
            return;
        }
        if let Some(id) = rest.strip_prefix("BIOMOLECULE:") {
            self.flush_current();
            let id = id.trim().to_string();
            if !id.is_empty() {
                self.current = Some(Assembly::new(id));
            }
            return;
        }
        if self.current.is_none() {
            return;
        }
        if let Some(details) = rest.strip_prefix("AUTHOR DETERMINED BIOLOGICAL UNIT:") {
            if let Some(asm) = self.current.as_mut() {
                asm.details = details.trim().to_string();
            }
            return;
        }
        if let Some(details) = rest.strip_prefix("SOFTWARE DETERMINED QUATERNARY STRUCTURE:") {
            if let Some(asm) = self.current.as_mut()
                && asm.details.is_empty()
            {
                asm.details = details.trim().to_string();
            }
            return;
        }
        if let Some(chains) = rest.strip_prefix("APPLY THE FOLLOWING TO CHAINS:") {
            self.flush_gen();
            self.pending_chains = parse_chain_list(chains);
            return;
        }
        if let Some(chains) = rest.strip_prefix("AND CHAINS:") {
            self.pending_chains.extend(parse_chain_list(chains));
            return;
        }
        if rest.contains("BIOMT") {
            parse_biomt_line(rest, &mut self.biomt_rows);
        }
    }

    fn flush_gen(&mut self) {
        let Some(asm) = self.current.as_mut() else {
            return;
        };
        if self.pending_chains.is_empty() || self.biomt_rows.is_empty() {
            return;
        }
        let mut ids: Vec<String> = self.biomt_rows.keys().cloned().collect();
        ids.sort_by(|a, b| match (a.parse::<i32>(), b.parse::<i32>()) {
            (Ok(x), Ok(y)) => x.cmp(&y),
            _ => a.cmp(b),
        });
        for id in &ids {
            if let Some(rows) = self.biomt_rows.get(id) {
                asm.operators.insert(id.clone(), affine_from_rows(*rows));
            }
        }
        asm.gens.push(AssemblyGen {
            oper_expression: ids.join(","),
            chain_ids: self.pending_chains.clone(),
        });
        self.pending_chains.clear();
        self.biomt_rows.clear();
    }

    fn flush_current(&mut self) {
        self.flush_gen();
        if let Some(asm) = self.current.take()
            && !asm.gens.is_empty()
        {
            self.assemblies.push(asm);
        }
        self.pending_chains.clear();
        self.biomt_rows.clear();
    }

    fn finish(mut self) -> Vec<Assembly> {
        self.flush_current();
        self.assemblies
    }
}

fn parse_chain_list(s: &str) -> Vec<String> {
    s.split(',')
        .map(|p| p.trim().trim_end_matches(':').to_string())
        .filter(|p| !p.is_empty())
        .collect()
}

fn parse_biomt_line(rest: &str, rows: &mut HashMap<String, [[f32; 4]; 3]>) {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    let Some(tag) = tokens.first() else {
        return;
    };
    let row_idx = if tag.ends_with('1') {
        0
    } else if tag.ends_with('2') {
        1
    } else if tag.ends_with('3') {
        2
    } else {
        return;
    };
    if tokens.len() < 6 {
        return;
    }
    let id = tokens[1].to_string();
    let Ok(a) = tokens[2].parse::<f32>() else {
        return;
    };
    let Ok(b) = tokens[3].parse::<f32>() else {
        return;
    };
    let Ok(c) = tokens[4].parse::<f32>() else {
        return;
    };
    let Ok(t) = tokens[5].parse::<f32>() else {
        return;
    };
    let entry = rows.entry(id).or_insert([[0.0; 4]; 3]);
    entry[row_idx] = [a, b, c, t];
}

/// Extracts the fixed-column substring `s[start..end]`.
///
/// PDB columns are byte offsets into an ASCII line. A range that falls outside
/// the line or lands inside a multi-byte UTF-8 character means the record is
/// malformed; treat that field as empty rather than panicking on a non-char
/// boundary slice.
fn safe_slice(s: &str, start: usize, end: usize) -> &str {
    if start >= s.len() || end <= start {
        return "";
    }
    s.get(start..end.min(s.len())).unwrap_or("")
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

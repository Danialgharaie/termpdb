//! mmCIF format parser.
//!
//! Parses Macromolecular Crystallographic Information Framework (mmCIF) text
//! into a [`Structure`]. Supports `_atom_site`, `_struct_conf`, `_struct_sheet_range`,
//! and top-level metadata tags.

use std::collections::{BTreeMap, HashMap};

use crate::error::{Result, TermPdbError};
use crate::math::{Mat4, Vec3};
use crate::model::assembly::affine_from_rows;
use crate::model::{Assembly, AssemblyGen, Atom, Residue, Structure, element_by_symbol};
use crate::parser::{ModelAccum, assemble_model, infer_element};

#[derive(Debug, PartialEq, Clone)]
enum CifToken {
    Data(String),
    Loop,
    Tag(String),
    Value(String),
}

/// Parses an mmCIF-formatted string into a [`Structure`].
pub fn parse_cif(input: &str) -> Result<Structure> {
    let tokens = tokenize_cif(input);
    let mut structure = Structure::new("");

    let mut helices: Vec<(String, i32, i32)> = Vec::new();
    let mut sheets: Vec<(String, i32, i32)> = Vec::new();

    let mut accums: BTreeMap<i32, ModelAccum> = BTreeMap::new();
    let mut assembly_meta: Vec<(String, String)> = Vec::new();
    let mut assembly_gens: Vec<(String, String, String)> = Vec::new();
    let mut oper_list: HashMap<String, Mat4> = HashMap::new();

    let mut token_idx = 0;
    while token_idx < tokens.len() {
        match &tokens[token_idx] {
            CifToken::Data(id) => {
                if structure.id_code.is_none() && !id.is_empty() {
                    structure.id_code = Some(id.clone());
                }
                token_idx += 1;
            }
            CifToken::Tag(tag) => {
                token_idx += 1;
                if let Some(CifToken::Value(val)) = tokens.get(token_idx) {
                    if tag == "_entry.id" && structure.id_code.is_none() && val != "?" && val != "."
                    {
                        structure.id_code = Some(val.clone());
                    } else if tag == "_struct.title" && val != "?" && val != "." {
                        structure.title = val.trim().to_string();
                    }
                    token_idx += 1;
                }
            }
            CifToken::Loop => {
                token_idx += 1;
                let mut headers = Vec::new();
                while token_idx < tokens.len() {
                    if let CifToken::Tag(tag) = &tokens[token_idx] {
                        headers.push(tag.clone());
                        token_idx += 1;
                    } else {
                        break;
                    }
                }

                if headers.is_empty() {
                    continue;
                }

                let num_cols = headers.len();
                let mut values = Vec::new();
                while token_idx < tokens.len() {
                    match &tokens[token_idx] {
                        CifToken::Value(val) => {
                            values.push(val.clone());
                            token_idx += 1;
                        }
                        _ => break,
                    }
                }

                // mmCIF loops must be rectangular: a ragged tail would be
                // silently dropped by integer division, corrupting every
                // downstream row. Reject instead of truncating. This single
                // check covers all loop types (rows are sliced here before
                // dispatching to the parse_*_loop helpers).
                if values.len() % num_cols != 0 {
                    return Err(TermPdbError::ParseError(format!(
                        "Malformed CIF loop: tag '{}' defines {} columns but {} values \
                         were collected; mmCIF loops must be rectangular",
                        headers[0],
                        num_cols,
                        values.len()
                    )));
                }

                let num_rows = values.len() / num_cols;
                let is_atom_site = headers.iter().any(|h| h.starts_with("_atom_site."));
                let is_struct_conf = headers.iter().any(|h| h.starts_with("_struct_conf."));
                let is_struct_sheet = headers
                    .iter()
                    .any(|h| h.starts_with("_struct_sheet_range."));
                let is_assembly = headers
                    .iter()
                    .any(|h| h.starts_with("_pdbx_struct_assembly."));
                let is_assembly_gen = headers
                    .iter()
                    .any(|h| h.starts_with("_pdbx_struct_assembly_gen."));
                let is_oper_list = headers
                    .iter()
                    .any(|h| h.starts_with("_pdbx_struct_oper_list."));

                if is_atom_site {
                    parse_atom_site_loop(&headers, &values, num_cols, num_rows, &mut accums)?;
                } else if is_struct_conf {
                    parse_struct_conf_loop(&headers, &values, num_cols, num_rows, &mut helices);
                } else if is_struct_sheet {
                    parse_struct_sheet_loop(&headers, &values, num_cols, num_rows, &mut sheets);
                } else if is_assembly_gen {
                    parse_assembly_gen_loop(
                        &headers,
                        &values,
                        num_cols,
                        num_rows,
                        &mut assembly_gens,
                    );
                } else if is_oper_list {
                    parse_oper_list_loop(&headers, &values, num_cols, num_rows, &mut oper_list);
                } else if is_assembly {
                    parse_assembly_loop(&headers, &values, num_cols, num_rows, &mut assembly_meta);
                }
            }
            CifToken::Value(_) => {
                token_idx += 1;
            }
        }
    }

    let has_ss = !helices.is_empty() || !sheets.is_empty();
    let mut models = Vec::with_capacity(accums.len());
    for (serial, accum) in accums {
        models.push(assemble_model(serial, accum, &helices, &sheets));
    }
    structure.set_models(models);
    structure.set_assemblies(build_cif_assemblies(
        assembly_meta,
        assembly_gens,
        oper_list,
    ));

    if !has_ss {
        crate::model::dssp::assign_dssp(&mut structure);
    }

    Ok(structure)
}

fn build_cif_assemblies(
    meta: Vec<(String, String)>,
    gens: Vec<(String, String, String)>,
    operators: HashMap<String, Mat4>,
) -> Vec<Assembly> {
    let mut by_id: HashMap<String, Assembly> = HashMap::new();
    for (id, details) in meta {
        let mut asm = Assembly::new(&id);
        asm.details = details;
        asm.operators = operators.clone();
        by_id.insert(id, asm);
    }
    for (assembly_id, expr, chains) in gens {
        let asm = by_id.entry(assembly_id.clone()).or_insert_with(|| {
            let mut a = Assembly::new(&assembly_id);
            a.operators = operators.clone();
            a
        });
        let chain_ids = chains
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty() && s != "." && s != "?")
            .collect();
        asm.gens.push(AssemblyGen {
            oper_expression: expr,
            chain_ids,
        });
    }
    let mut assemblies: Vec<Assembly> = by_id.into_values().collect();
    assemblies.sort_by(|a, b| match (a.id.parse::<i32>(), b.id.parse::<i32>()) {
        (Ok(x), Ok(y)) => x.cmp(&y),
        _ => a.id.cmp(&b.id),
    });
    assemblies
        .into_iter()
        .filter(|a| !a.gens.is_empty())
        .collect()
}

fn col_index(headers: &[String], name: &str) -> Option<usize> {
    headers.iter().position(|h| h.eq_ignore_ascii_case(name))
}

fn get_val(row: &[String], idx: Option<usize>) -> Option<&str> {
    idx.and_then(|i| row.get(i).map(|s| s.as_str()))
}

fn parse_assembly_loop(
    headers: &[String],
    values: &[String],
    num_cols: usize,
    num_rows: usize,
    out: &mut Vec<(String, String)>,
) {
    let id_col = col_index(headers, "_pdbx_struct_assembly.id");
    let details_col = col_index(headers, "_pdbx_struct_assembly.details")
        .or_else(|| col_index(headers, "_pdbx_struct_assembly.method_details"));
    for r in 0..num_rows {
        let row = &values[r * num_cols..(r + 1) * num_cols];
        let id = get_val(row, id_col).unwrap_or("");
        if id.is_empty() || id == "." || id == "?" {
            continue;
        }
        let details = get_val(row, details_col)
            .filter(|s| *s != "." && *s != "?")
            .unwrap_or("")
            .to_string();
        out.push((id.to_string(), details));
    }
}

fn parse_assembly_gen_loop(
    headers: &[String],
    values: &[String],
    num_cols: usize,
    num_rows: usize,
    out: &mut Vec<(String, String, String)>,
) {
    let id_col = col_index(headers, "_pdbx_struct_assembly_gen.assembly_id");
    let expr_col = col_index(headers, "_pdbx_struct_assembly_gen.oper_expression");
    let chains_col = col_index(headers, "_pdbx_struct_assembly_gen.asym_id_list");
    for r in 0..num_rows {
        let row = &values[r * num_cols..(r + 1) * num_cols];
        let id = get_val(row, id_col).unwrap_or("");
        let expr = get_val(row, expr_col).unwrap_or("");
        let chains = get_val(row, chains_col).unwrap_or("");
        if id.is_empty() || id == "." || expr.is_empty() || expr == "." {
            continue;
        }
        out.push((id.to_string(), expr.to_string(), chains.to_string()));
    }
}

fn parse_oper_list_loop(
    headers: &[String],
    values: &[String],
    num_cols: usize,
    num_rows: usize,
    out: &mut HashMap<String, Mat4>,
) {
    let id_col = col_index(headers, "_pdbx_struct_oper_list.id");
    let m11 = col_index(headers, "_pdbx_struct_oper_list.matrix[1][1]");
    let m12 = col_index(headers, "_pdbx_struct_oper_list.matrix[1][2]");
    let m13 = col_index(headers, "_pdbx_struct_oper_list.matrix[1][3]");
    let m21 = col_index(headers, "_pdbx_struct_oper_list.matrix[2][1]");
    let m22 = col_index(headers, "_pdbx_struct_oper_list.matrix[2][2]");
    let m23 = col_index(headers, "_pdbx_struct_oper_list.matrix[2][3]");
    let m31 = col_index(headers, "_pdbx_struct_oper_list.matrix[3][1]");
    let m32 = col_index(headers, "_pdbx_struct_oper_list.matrix[3][2]");
    let m33 = col_index(headers, "_pdbx_struct_oper_list.matrix[3][3]");
    let v1 = col_index(headers, "_pdbx_struct_oper_list.vector[1]");
    let v2 = col_index(headers, "_pdbx_struct_oper_list.vector[2]");
    let v3 = col_index(headers, "_pdbx_struct_oper_list.vector[3]");

    let f = |row: &[String], col: Option<usize>, default: f32| -> f32 {
        get_val(row, col)
            .and_then(|s| s.parse::<f32>().ok())
            .unwrap_or(default)
    };

    for r in 0..num_rows {
        let row = &values[r * num_cols..(r + 1) * num_cols];
        let id = get_val(row, id_col).unwrap_or("");
        if id.is_empty() || id == "." || id == "?" {
            continue;
        }
        let rows = [
            [
                f(row, m11, 1.0),
                f(row, m12, 0.0),
                f(row, m13, 0.0),
                f(row, v1, 0.0),
            ],
            [
                f(row, m21, 0.0),
                f(row, m22, 1.0),
                f(row, m23, 0.0),
                f(row, v2, 0.0),
            ],
            [
                f(row, m31, 0.0),
                f(row, m32, 0.0),
                f(row, m33, 1.0),
                f(row, v3, 0.0),
            ],
        ];
        out.insert(id.to_string(), affine_from_rows(rows));
    }
}

fn parse_atom_site_loop(
    headers: &[String],
    values: &[String],
    num_cols: usize,
    num_rows: usize,
    accums: &mut BTreeMap<i32, ModelAccum>,
) -> Result<()> {
    let group_pdb_col = col_index(headers, "_atom_site.group_PDB");
    let id_col = col_index(headers, "_atom_site.id");
    let type_symbol_col = col_index(headers, "_atom_site.type_symbol");
    let label_atom_id_col = col_index(headers, "_atom_site.label_atom_id");
    let auth_atom_id_col = col_index(headers, "_atom_site.auth_atom_id");
    let label_alt_id_col = col_index(headers, "_atom_site.label_alt_id");
    let label_comp_id_col = col_index(headers, "_atom_site.label_comp_id");
    let auth_comp_id_col = col_index(headers, "_atom_site.auth_comp_id");
    let label_asym_id_col = col_index(headers, "_atom_site.label_asym_id");
    let auth_asym_id_col = col_index(headers, "_atom_site.auth_asym_id");
    let label_seq_id_col = col_index(headers, "_atom_site.label_seq_id");
    let auth_seq_id_col = col_index(headers, "_atom_site.auth_seq_id");
    let cartn_x_col = col_index(headers, "_atom_site.Cartn_x");
    let cartn_y_col = col_index(headers, "_atom_site.Cartn_y");
    let cartn_z_col = col_index(headers, "_atom_site.Cartn_z");
    let occupancy_col = col_index(headers, "_atom_site.occupancy");
    let b_iso_col = col_index(headers, "_atom_site.B_iso_or_equiv");
    let formal_charge_col = col_index(headers, "_atom_site.pdbx_formal_charge");
    let ins_code_col = col_index(headers, "_atom_site.pdbx_PDB_ins_code");
    let model_num_col = col_index(headers, "_atom_site.pdbx_PDB_model_num");

    for r in 0..num_rows {
        let row = &values[r * num_cols..(r + 1) * num_cols];

        let group_pdb = get_val(row, group_pdb_col).unwrap_or("ATOM");
        let is_hetatm = group_pdb.eq_ignore_ascii_case("HETATM");

        let model_serial = get_val(row, model_num_col)
            .and_then(|s| {
                if s == "?" || s == "." || s.is_empty() {
                    None
                } else {
                    s.parse::<i32>().ok()
                }
            })
            .unwrap_or(1);

        let accum = accums.entry(model_serial).or_default();

        let serial = get_val(row, id_col)
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or((accum.atoms.len() + 1) as i32);

        let atom_name = get_val(row, auth_atom_id_col)
            .filter(|s| *s != "?" && *s != ".")
            .or_else(|| get_val(row, label_atom_id_col))
            .unwrap_or("X");

        let res_name = get_val(row, auth_comp_id_col)
            .filter(|s| *s != "?" && *s != ".")
            .or_else(|| get_val(row, label_comp_id_col))
            .unwrap_or("UNK");

        let chain_id_str = get_val(row, auth_asym_id_col)
            .filter(|s| *s != "?" && *s != ".")
            .or_else(|| get_val(row, label_asym_id_col))
            .unwrap_or("A");
        let chain_id = if chain_id_str.is_empty() {
            "A".to_string()
        } else {
            chain_id_str.to_string()
        };

        let res_seq = get_val(row, auth_seq_id_col)
            .and_then(|s| s.parse::<i32>().ok())
            .or_else(|| get_val(row, label_seq_id_col).and_then(|s| s.parse::<i32>().ok()))
            .unwrap_or(1);

        let ins_code = get_val(row, ins_code_col)
            .filter(|s| *s != "?" && *s != ".")
            .and_then(|s| s.chars().next());

        let x = get_val(row, cartn_x_col)
            .and_then(|s| s.parse::<f32>().ok())
            .ok_or_else(|| {
                TermPdbError::ParseError(format!("Missing or invalid Cartn_x at row {}", r))
            })?;
        let y = get_val(row, cartn_y_col)
            .and_then(|s| s.parse::<f32>().ok())
            .ok_or_else(|| {
                TermPdbError::ParseError(format!("Missing or invalid Cartn_y at row {}", r))
            })?;
        let z = get_val(row, cartn_z_col)
            .and_then(|s| s.parse::<f32>().ok())
            .ok_or_else(|| {
                TermPdbError::ParseError(format!("Missing or invalid Cartn_z at row {}", r))
            })?;

        let occupancy = get_val(row, occupancy_col)
            .and_then(|s| s.parse::<f32>().ok())
            .unwrap_or(1.0);
        let b_factor = get_val(row, b_iso_col)
            .and_then(|s| s.parse::<f32>().ok())
            .unwrap_or(0.0);

        let elem_sym = get_val(row, type_symbol_col)
            .filter(|s| *s != "?" && *s != ".")
            .unwrap_or("");
        let element = if !elem_sym.is_empty() {
            let elem = element_by_symbol(elem_sym);
            if elem.atomic_number != 0 {
                elem
            } else {
                infer_element(atom_name, res_name, is_hetatm)
            }
        } else {
            infer_element(atom_name, res_name, is_hetatm)
        };

        let alt_loc = get_val(row, label_alt_id_col)
            .filter(|s| *s != "?" && *s != ".")
            .and_then(|s| s.chars().next());

        let charge = get_val(row, formal_charge_col).and_then(|s| s.parse::<i8>().ok());

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
            last_res.seq == res_seq && last_res.ins_code == ins_code && last_res.name == res_name
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

    Ok(())
}

fn parse_struct_conf_loop(
    headers: &[String],
    values: &[String],
    num_cols: usize,
    num_rows: usize,
    helices: &mut Vec<(String, i32, i32)>,
) {
    let conf_type_col = col_index(headers, "_struct_conf.conf_type_id");
    let beg_auth_asym_col = col_index(headers, "_struct_conf.beg_auth_asym_id");
    let beg_label_asym_col = col_index(headers, "_struct_conf.beg_label_asym_id");
    let beg_auth_seq_col = col_index(headers, "_struct_conf.beg_auth_seq_id");
    let beg_label_seq_col = col_index(headers, "_struct_conf.beg_label_seq_id");

    let end_auth_seq_col = col_index(headers, "_struct_conf.end_auth_seq_id");
    let end_label_seq_col = col_index(headers, "_struct_conf.end_label_seq_id");

    for r in 0..num_rows {
        let row = &values[r * num_cols..(r + 1) * num_cols];
        let conf_type = get_val(row, conf_type_col).unwrap_or("HELX_P");
        if conf_type.starts_with("HELX") {
            let beg_chain = get_val(row, beg_auth_asym_col)
                .filter(|s| *s != "?" && *s != ".")
                .or_else(|| get_val(row, beg_label_asym_col))
                .unwrap_or("A");
            let beg_seq = get_val(row, beg_auth_seq_col)
                .and_then(|s| s.parse::<i32>().ok())
                .or_else(|| get_val(row, beg_label_seq_col).and_then(|s| s.parse::<i32>().ok()))
                .unwrap_or(0);
            let end_seq = get_val(row, end_auth_seq_col)
                .and_then(|s| s.parse::<i32>().ok())
                .or_else(|| get_val(row, end_label_seq_col).and_then(|s| s.parse::<i32>().ok()))
                .unwrap_or(0);

            helices.push((beg_chain.to_string(), beg_seq, end_seq));
        }
    }
}

fn parse_struct_sheet_loop(
    headers: &[String],
    values: &[String],
    num_cols: usize,
    num_rows: usize,
    sheets: &mut Vec<(String, i32, i32)>,
) {
    let beg_auth_asym_col = col_index(headers, "_struct_sheet_range.beg_auth_asym_id");
    let beg_label_asym_col = col_index(headers, "_struct_sheet_range.beg_label_asym_id");
    let beg_auth_seq_col = col_index(headers, "_struct_sheet_range.beg_auth_seq_id");
    let beg_label_seq_col = col_index(headers, "_struct_sheet_range.beg_label_seq_id");

    let end_auth_seq_col = col_index(headers, "_struct_sheet_range.end_auth_seq_id");
    let end_label_seq_col = col_index(headers, "_struct_sheet_range.end_label_seq_id");

    for r in 0..num_rows {
        let row = &values[r * num_cols..(r + 1) * num_cols];
        let beg_chain = get_val(row, beg_auth_asym_col)
            .filter(|s| *s != "?" && *s != ".")
            .or_else(|| get_val(row, beg_label_asym_col))
            .unwrap_or("A");
        let beg_seq = get_val(row, beg_auth_seq_col)
            .and_then(|s| s.parse::<i32>().ok())
            .or_else(|| get_val(row, beg_label_seq_col).and_then(|s| s.parse::<i32>().ok()))
            .unwrap_or(0);
        let end_seq = get_val(row, end_auth_seq_col)
            .and_then(|s| s.parse::<i32>().ok())
            .or_else(|| get_val(row, end_label_seq_col).and_then(|s| s.parse::<i32>().ok()))
            .unwrap_or(0);

        sheets.push((beg_chain.to_string(), beg_seq, end_seq));
    }
}

fn tokenize_cif(input: &str) -> Vec<CifToken> {
    let chars: Vec<char> = input.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;

    while i < chars.len() {
        let c = chars[i];

        if c.is_whitespace() {
            i += 1;
            continue;
        }

        if c == '#' {
            // Comment: skip until newline
            while i < chars.len() && chars[i] != '\n' {
                i += 1;
            }
            continue;
        }

        // Check for semicolon-delimited multiline string at start of line
        if c == ';'
            && (i == 0
                || chars[i - 1] == '\n'
                || (i > 1 && chars[i - 1] == '\r' && chars[i - 2] == '\n'))
        {
            i += 1; // skip leading ';'
            let start_val = i;
            while i < chars.len() {
                if chars[i] == ';'
                    && (chars[i - 1] == '\n'
                        || (i > 1 && chars[i - 1] == '\r' && chars[i - 2] == '\n'))
                {
                    break;
                }
                i += 1;
            }
            let val_str: String = chars[start_val..i].iter().collect();
            if i < chars.len() && chars[i] == ';' {
                i += 1; // skip trailing ';'
            }
            tokens.push(CifToken::Value(val_str.trim().to_string()));
            continue;
        }

        if c == '\'' || c == '"' {
            let quote = c;
            i += 1;
            let start_val = i;
            while i < chars.len() && chars[i] != quote {
                i += 1;
            }
            let val_str: String = chars[start_val..i].iter().collect();
            if i < chars.len() && chars[i] == quote {
                i += 1;
            }
            tokens.push(CifToken::Value(val_str));
            continue;
        }

        if c == '_' {
            let start_tag = i;
            while i < chars.len() && !chars[i].is_whitespace() && chars[i] != '#' {
                i += 1;
            }
            let tag_str: String = chars[start_tag..i].iter().collect();
            tokens.push(CifToken::Tag(tag_str));
            continue;
        }

        // Regular word token. `data_...` is only a data-block marker when it
        // begins a line (column 0); mid-line occurrences are ordinary values,
        // otherwise a stray `data_foo` inside loop values would terminate
        // value collection early and silently shift every subsequent row.
        let start_word = i;
        while i < chars.len() && !chars[i].is_whitespace() && chars[i] != '#' {
            i += 1;
        }
        let word: String = chars[start_word..i].iter().collect();
        let at_line_start = start_word == 0
            || chars[start_word - 1] == '\n'
            || (start_word > 1 && chars[start_word - 1] == '\r' && chars[start_word - 2] == '\n');
        if word.eq_ignore_ascii_case("loop_") {
            tokens.push(CifToken::Loop);
        } else if at_line_start && word.to_ascii_lowercase().starts_with("data_") {
            tokens.push(CifToken::Data(word[5..].to_string()));
        } else {
            tokens.push(CifToken::Value(word));
        }
    }

    tokens
}

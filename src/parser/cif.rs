//! mmCIF format parser.
//!
//! Parses Macromolecular Crystallographic Information Framework (mmCIF) text
//! into a [`Structure`]. Supports `_atom_site`, `_struct_conf`, `_struct_sheet_range`,
//! and top-level metadata tags.

use std::collections::HashMap;

use crate::error::{Result, TermPdbError};
use crate::math::Vec3;
use crate::model::{
    Atom, Chain, Element, Residue, SecondaryStructure, Structure, element_by_symbol,
};

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

    let mut chain_order: Vec<String> = Vec::new();
    let mut chain_residues: HashMap<String, Vec<Residue>> = HashMap::new();

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

                let num_rows = values.len() / num_cols;
                let is_atom_site = headers.iter().any(|h| h.starts_with("_atom_site."));
                let is_struct_conf = headers.iter().any(|h| h.starts_with("_struct_conf."));
                let is_struct_sheet = headers
                    .iter()
                    .any(|h| h.starts_with("_struct_sheet_range."));

                if is_atom_site {
                    parse_atom_site_loop(
                        &headers,
                        &values,
                        num_cols,
                        num_rows,
                        &mut structure,
                        &mut chain_order,
                        &mut chain_residues,
                    )?;
                } else if is_struct_conf {
                    parse_struct_conf_loop(&headers, &values, num_cols, num_rows, &mut helices);
                } else if is_struct_sheet {
                    parse_struct_sheet_loop(&headers, &values, num_cols, num_rows, &mut sheets);
                }
            }
            CifToken::Value(_) => {
                token_idx += 1;
            }
        }
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

    Ok(structure)
}

fn col_index(headers: &[String], name: &str) -> Option<usize> {
    headers.iter().position(|h| h.eq_ignore_ascii_case(name))
}

fn get_val(row: &[String], idx: Option<usize>) -> Option<&str> {
    idx.and_then(|i| row.get(i).map(|s| s.as_str()))
}

#[allow(clippy::too_many_arguments)]
fn parse_atom_site_loop(
    headers: &[String],
    values: &[String],
    num_cols: usize,
    num_rows: usize,
    structure: &mut Structure,
    chain_order: &mut Vec<String>,
    chain_residues: &mut HashMap<String, Vec<Residue>>,
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

    for r in 0..num_rows {
        let row = &values[r * num_cols..(r + 1) * num_cols];

        let group_pdb = get_val(row, group_pdb_col).unwrap_or("ATOM");
        let is_hetatm = group_pdb.eq_ignore_ascii_case("HETATM");

        let serial = get_val(row, id_col)
            .and_then(|s| s.parse::<i32>().ok())
            .unwrap_or((structure.atoms.len() + 1) as i32);

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
                infer_cif_element(atom_name, res_name, is_hetatm)
            }
        } else {
            infer_cif_element(atom_name, res_name, is_hetatm)
        };

        let alt_loc = get_val(row, label_alt_id_col)
            .filter(|s| *s != "?" && *s != ".")
            .and_then(|s| s.chars().next());

        let charge = get_val(row, formal_charge_col).and_then(|s| s.parse::<i8>().ok());

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

        // Add to chain/residue structure
        if !chain_order.contains(&chain_id) {
            chain_order.push(chain_id.clone());
            chain_residues.insert(chain_id.clone(), Vec::new());
        }

        let residues = chain_residues.get_mut(&chain_id).unwrap();
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

fn infer_cif_element(atom_name: &str, res_name: &str, is_hetatm: bool) -> Element {
    let trimmed = atom_name.trim();
    if trimmed.is_empty() {
        return Element::unknown();
    }

    let cleaned = trimmed.trim_start_matches(|c: char| c.is_ascii_digit());
    if cleaned.is_empty() {
        return Element::unknown();
    }

    let res_trimmed = res_name.trim();
    if cleaned.eq_ignore_ascii_case("CA")
        && !is_hetatm
        && !res_trimmed.eq_ignore_ascii_case("CA")
        && !res_trimmed.eq_ignore_ascii_case("CAL")
    {
        return element_by_symbol("C");
    }

    if cleaned.len() >= 2 {
        let first_char = cleaned.chars().next().unwrap();
        let second_char = cleaned.chars().nth(1).unwrap().to_ascii_uppercase();

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

    let candidate1 = &cleaned[..1];
    let elem1 = element_by_symbol(candidate1);
    if elem1.atomic_number != 0 {
        return elem1;
    }

    Element::unknown()
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

        // Regular word token
        let start_word = i;
        while i < chars.len() && !chars[i].is_whitespace() && chars[i] != '#' {
            i += 1;
        }
        let word: String = chars[start_word..i].iter().collect();
        if word.eq_ignore_ascii_case("loop_") {
            tokens.push(CifToken::Loop);
        } else if word.to_ascii_lowercase().starts_with("data_") {
            tokens.push(CifToken::Data(word[5..].to_string()));
        } else {
            tokens.push(CifToken::Value(word));
        }
    }

    tokens
}

//! Biological assemblies: symmetry operators and expansion of the asymmetric unit.

use std::collections::{HashMap, HashSet};

use crate::error::{Result, TermPdbError};
use crate::math::{Mat4, Vec3};
use crate::model::Model;
use crate::model::bond::BondDetector;
use crate::model::chain::Chain;

/// One generation rule: apply an operator expression to a set of chains.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyGen {
    pub oper_expression: String,
    pub chain_ids: Vec<String>,
}

/// A biological assembly (biomolecule) from REMARK 350 or `_pdbx_struct_assembly`.
#[derive(Debug, Clone, PartialEq)]
pub struct Assembly {
    /// File identifier (`1`, `2`, …), not a compacted index.
    pub id: String,
    /// Free-text description (e.g. "author_defined_assembly", "DIMERIC").
    pub details: String,
    pub gens: Vec<AssemblyGen>,
    /// Operators keyed by id as written in the file (`1`, `P_1`, …).
    pub operators: HashMap<String, Mat4>,
}

impl Assembly {
    pub fn new(id: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            details: String::new(),
            gens: Vec::new(),
            operators: HashMap::new(),
        }
    }
}

/// Parses an mmCIF/PDB operator expression into groups of operator ids.
///
/// `1,2` and `(1,2)` are one group (two copies). `(1)(2)` is a product
/// (compose the two operators into one copy).
pub fn parse_oper_expression(expr: &str) -> Result<Vec<Vec<String>>> {
    let trimmed = expr.trim();
    if trimmed.is_empty() || trimmed == "." || trimmed == "?" {
        return Err(TermPdbError::ParseError(
            "Empty assembly operator expression".to_string(),
        ));
    }

    let mut groups: Vec<Vec<String>> = Vec::new();
    let mut current = String::new();
    let mut in_paren = false;

    let flush_group = |raw: &str, groups: &mut Vec<Vec<String>>| -> Result<()> {
        let ids = expand_oper_group(raw)?;
        if !ids.is_empty() {
            groups.push(ids);
        }
        Ok(())
    };

    for c in trimmed.chars() {
        match c {
            '(' => {
                if !current.trim().is_empty() {
                    flush_group(&current, &mut groups)?;
                    current.clear();
                }
                in_paren = true;
            }
            ')' => {
                flush_group(&current, &mut groups)?;
                current.clear();
                in_paren = false;
            }
            _ => current.push(c),
        }
    }

    if in_paren {
        return Err(TermPdbError::ParseError(format!(
            "Unbalanced parentheses in operator expression '{expr}'"
        )));
    }
    if !current.trim().is_empty() {
        flush_group(&current, &mut groups)?;
    }

    if groups.is_empty() {
        return Err(TermPdbError::ParseError(format!(
            "No operators in expression '{expr}'"
        )));
    }
    Ok(groups)
}

fn expand_oper_group(raw: &str) -> Result<Vec<String>> {
    let mut ids = Vec::new();
    for token in raw.split(',') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        if let Some((a, b)) = token.split_once('-') {
            let a = a.trim();
            let b = b.trim();
            if let (Ok(start), Ok(end)) = (a.parse::<i32>(), b.parse::<i32>()) {
                if end < start {
                    return Err(TermPdbError::ParseError(format!(
                        "Invalid operator range '{token}'"
                    )));
                }
                for n in start..=end {
                    ids.push(n.to_string());
                }
                continue;
            }
        }
        ids.push(token.to_string());
    }
    Ok(ids)
}

fn cartesian_product(groups: &[Vec<String>]) -> Vec<Vec<String>> {
    let mut acc: Vec<Vec<String>> = vec![Vec::new()];
    for group in groups {
        let mut next = Vec::new();
        for prefix in &acc {
            for id in group {
                let mut row = prefix.clone();
                row.push(id.clone());
                next.push(row);
            }
        }
        acc = next;
    }
    acc
}

fn compose_operators(ids: &[String], operators: &HashMap<String, Mat4>) -> Result<Mat4> {
    let mut m = Mat4::identity();
    for id in ids.iter().rev() {
        let Some(op) = operators.get(id) else {
            return Err(TermPdbError::InvalidStructure(format!(
                "Assembly operator '{id}' not found"
            )));
        };
        m = op.mul(&m);
    }
    Ok(m)
}

/// Matrices to apply (one copy per matrix) for an operator expression.
pub fn assembly_transforms(
    expr: &str,
    operators: &HashMap<String, Mat4>,
) -> Result<Vec<(Vec<String>, Mat4)>> {
    let groups = parse_oper_expression(expr)?;
    let combos = cartesian_product(&groups);
    let mut out = Vec::with_capacity(combos.len());
    for ids in combos {
        let mat = compose_operators(&ids, operators)?;
        out.push((ids, mat));
    }
    Ok(out)
}

/// Builds a new model by applying `assembly` to the asymmetric-unit model.
pub fn expand_model(asu: &Model, assembly: &Assembly) -> Result<Model> {
    let mut out = Model::new(asu.serial);
    let mut used_chains: HashSet<String> = HashSet::new();

    for rule in &assembly.gens {
        let transforms = assembly_transforms(&rule.oper_expression, &assembly.operators)?;
        for chain_id in &rule.chain_ids {
            let Some(src_chain) = asu
                .chains
                .iter()
                .find(|c| c.id.eq_ignore_ascii_case(chain_id))
            else {
                continue;
            };
            for (_ids, mat) in &transforms {
                let dest_id = unique_chain_id(&src_chain.id, &mut used_chains);
                copy_chain_transformed(asu, src_chain, &dest_id, mat, &mut out);
            }
        }
    }

    out.bonds = BondDetector::detect_bonds(&out.atoms);
    Ok(out)
}

fn unique_chain_id(base: &str, used: &mut HashSet<String>) -> String {
    if used.insert(base.to_string()) {
        return base.to_string();
    }
    let mut n = 2;
    loop {
        let id = format!("{base}{n}");
        if used.insert(id.clone()) {
            return id;
        }
        n += 1;
    }
}

fn copy_chain_transformed(
    asu: &Model,
    src_chain: &Chain,
    dest_id: &str,
    mat: &Mat4,
    out: &mut Model,
) {
    let mut new_chain = Chain::new(dest_id);
    for res in &src_chain.residues {
        let mut new_res = res.clone();
        new_res.chain_id = dest_id.to_string();
        new_res.atom_indices.clear();
        for &old_idx in &res.atom_indices {
            let Some(src) = asu.atoms.get(old_idx) else {
                continue;
            };
            let mut atom = src.clone();
            atom.pos = mat.transform_point(src.pos);
            atom.chain_id = dest_id.to_string();
            atom.index = out.atoms.len();
            new_res.atom_indices.push(atom.index);
            out.atoms.push(atom);
        }
        if !new_res.atom_indices.is_empty() {
            new_chain.residues.push(new_res);
        }
    }
    if !new_chain.residues.is_empty() {
        out.chains.push(new_chain);
    }
}

/// Builds a 4×4 affine transform from a row-major 3×3 rotation and translation.
pub fn affine_from_rows(rows: [[f32; 4]; 3]) -> Mat4 {
    Mat4::from_rotation_translation(
        [
            [rows[0][0], rows[0][1], rows[0][2]],
            [rows[1][0], rows[1][1], rows[1][2]],
            [rows[2][0], rows[2][1], rows[2][2]],
        ],
        Vec3::new(rows[0][3], rows[1][3], rows[2][3]),
    )
}

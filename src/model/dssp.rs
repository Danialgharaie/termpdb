//! Pure Rust DSSP Secondary Structure Assignment.
//!
//! Implements the Kabsch-Sander (DSSP) electrostatic backbone hydrogen-bonding
//! energy calculation ($E < -0.5\text{ kcal/mol}$) to automatically detect
//! $\alpha$-helices, $3_{10}$-helices, and $\beta$-sheets.

use crate::math::Vec3;
use crate::model::{SecondaryStructure, Structure};

/// Electrostatic constant $q_1 q_2 \cdot 332 = 0.42 \cdot 0.20 \cdot 332.0 = 27.888\text{ kcal}\cdot\text{\AA}/\text{mol}$.
const DSSP_COUPLING_CONSTANT: f32 = 27.888;
const DSSP_HBOND_CUTOFF: f32 = -0.5;

/// Calculates the electrostatic hydrogen bond energy between donor $N-H$ and acceptor $C=O$ in kcal/mol.
pub fn calculate_dssp_hbond_energy(c: Vec3, o: Vec3, n: Vec3, h: Vec3) -> f32 {
    let r_on = o.distance(&n).max(0.5);
    let r_ch = c.distance(&h).max(0.5);
    let r_oh = o.distance(&h).max(0.5);
    let r_cn = c.distance(&n).max(0.5);

    DSSP_COUPLING_CONSTANT * ((1.0 / r_on) + (1.0 / r_ch) - (1.0 / r_oh) - (1.0 / r_cn))
}

struct BackboneResidue {
    chain_idx: usize,
    res_idx: usize,
    c_pos: Option<Vec3>,
    o_pos: Option<Vec3>,
    n_pos: Option<Vec3>,
    h_pos: Option<Vec3>,
}

/// Runs DSSP secondary structure assignment on `structure` and updates `SecondaryStructure` on all residues.
/// Returns the total number of residues assigned to Helix or Sheet.
#[allow(clippy::needless_range_loop)]
pub fn assign_dssp(structure: &mut Structure) -> usize {
    let atoms = structure.atoms().to_vec();
    let mut backbone = Vec::new();

    for (c_i, chain) in structure.chains().iter().enumerate() {
        let mut prev_c: Option<Vec3> = None;

        for (r_i, res) in chain.residues.iter().enumerate() {
            let mut n_pos = None;
            let mut ca_pos = None;
            let mut c_pos = None;
            let mut o_pos = None;
            let mut explicit_h = None;

            for &atom_idx in &res.atom_indices {
                if let Some(atom) = atoms.get(atom_idx) {
                    match atom.name.as_str() {
                        "N" => n_pos = Some(atom.pos),
                        "CA" => ca_pos = Some(atom.pos),
                        "C" => c_pos = Some(atom.pos),
                        "O" => o_pos = Some(atom.pos),
                        "H" | "HN" | "H1" => explicit_h = Some(atom.pos),
                        _ => {}
                    }
                }
            }

            // Approximate backbone H position if missing
            let h_pos = explicit_h.or_else(|| {
                let n = n_pos?;
                let ca = ca_pos?;
                if let Some(c_prev) = prev_c {
                    let u = (n - c_prev).normalize();
                    let v = (n - ca).normalize();
                    let dir = (u + v).normalize();
                    Some(n + dir * 1.01)
                } else {
                    let dir = (n - ca).normalize();
                    Some(n + dir * 1.01)
                }
            });

            if let Some(c) = c_pos {
                prev_c = Some(c);
            }

            backbone.push(BackboneResidue {
                chain_idx: c_i,
                res_idx: r_i,
                c_pos,
                o_pos,
                n_pos,
                h_pos,
            });
        }
    }

    let n = backbone.len();
    if n < 4 {
        return 0;
    }

    // Build H-bond matrix: hb[i][j] is true if residue i (acceptor C=O) forms H-bond to j (donor N-H)
    let mut hb = vec![vec![false; n]; n];

    for i in 0..n {
        let acc_c = backbone[i].c_pos;
        let acc_o = backbone[i].o_pos;
        if let (Some(c), Some(o)) = (acc_c, acc_o) {
            for j in 0..n {
                if (i as isize - j as isize).abs() < 2 {
                    continue; // Skip same or adjacent residue
                }
                let don_n = backbone[j].n_pos;
                let don_h = backbone[j].h_pos;
                if let (Some(n_pos), Some(h_pos)) = (don_n, don_h) {
                    let energy = calculate_dssp_hbond_energy(c, o, n_pos, h_pos);
                    if energy < DSSP_HBOND_CUTOFF {
                        hb[i][j] = true;
                    }
                }
            }
        }
    }

    let mut ss = vec![SecondaryStructure::Coil; n];

    // 1. Detect 4-turn Alpha Helices: hb[i][i+4]
    for i in 0..(n.saturating_sub(4)) {
        if hb[i][i + 4] {
            for k in (i + 1)..=(i + 4) {
                if k < n && ss[k] == SecondaryStructure::Coil {
                    ss[k] = SecondaryStructure::Helix;
                }
            }
        }
    }

    // 2. Detect 3-turn 3_10 Helices: hb[i][i+3]
    for i in 0..(n.saturating_sub(3)) {
        if hb[i][i + 3] {
            for k in (i + 1)..=(i + 3) {
                if k < n && ss[k] == SecondaryStructure::Coil {
                    ss[k] = SecondaryStructure::Helix;
                }
            }
        }
    }

    // 3. Detect Beta Sheets (Antiparallel and Parallel ladders)
    for i in 0..n {
        for j in 0..n {
            if i >= j || (i as isize - j as isize).abs() < 3 {
                continue;
            }
            // Antiparallel bridge: (i->j and j->i) or (i-1->j+1 and j-1->i+1)
            let ap_bridge = (hb[i][j] && hb[j][i])
                || (i > 0
                    && j + 1 < n
                    && j > 0
                    && i + 1 < n
                    && hb[i - 1][j + 1]
                    && hb[j - 1][i + 1]);

            // Parallel bridge: (i-1->j and j->i+1) or (j-1->i and i->j+1)
            let p_bridge = (i > 0 && i + 1 < n && hb[i - 1][j] && hb[j][i + 1])
                || (j > 0 && j + 1 < n && hb[j - 1][i] && hb[i][j + 1]);

            if ap_bridge || p_bridge {
                ss[i] = SecondaryStructure::Sheet;
                ss[j] = SecondaryStructure::Sheet;
            }
        }
    }

    // Apply assigned secondary structure back to structure residues
    let mut assigned_count = 0;
    for (k, state) in ss.into_iter().enumerate() {
        let b = &backbone[k];
        if state != SecondaryStructure::Coil {
            assigned_count += 1;
        }
        if let Some(chain) = structure.chains_mut().get_mut(b.chain_idx)
            && let Some(res) = chain.residues.get_mut(b.res_idx)
        {
            res.secondary_structure = state;
        }
    }

    assigned_count
}

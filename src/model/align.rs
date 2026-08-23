//! Needleman-Wunsch Sequence Alignment & 3D Structural Superposition.

use crate::math::Vec3;
use crate::math::kabsch::{KabschResult, kabsch_align};
use crate::model::Structure;

/// Result of superimposing two macromolecular structures.
#[derive(Debug, Clone, PartialEq)]
pub struct SuperpositionResult {
    /// Kabsch rotation matrix and translation vector
    pub kabsch: KabschResult,
    /// Per-residue coordinate distance after alignment for all matched residue pairs
    pub per_residue_rmsd: Vec<f32>,
    /// Total number of aligned residue pairs
    pub aligned_pairs: usize,
}

/// Converts a 3-letter amino acid residue name to a 1-letter code (or 'X' if unknown).
pub fn residue_to_one_letter(res_name: &str) -> char {
    match res_name.trim().to_ascii_uppercase().as_str() {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        // Nucleic acids
        "DA" | "A" => 'A',
        "DC" | "C" => 'C',
        "DG" | "G" => 'G',
        "DT" | "T" => 'T',
        "DU" | "U" => 'U',
        _ => 'X',
    }
}

/// Needleman-Wunsch global sequence alignment between two strings.
/// Returns a vector of matched 0-based indices `(Option<index_a>, Option<index_b>)`.
#[allow(clippy::needless_range_loop)]
pub fn needleman_wunsch(seq_a: &str, seq_b: &str) -> Vec<(Option<usize>, Option<usize>)> {
    let chars_a: Vec<char> = seq_a.chars().collect();
    let chars_b: Vec<char> = seq_b.chars().collect();

    let n = chars_a.len();
    let m = chars_b.len();

    if n == 0 && m == 0 {
        return Vec::new();
    }
    if n == 0 {
        return (0..m).map(|j| (None, Some(j))).collect();
    }
    if m == 0 {
        return (0..n).map(|i| (Some(i), None)).collect();
    }

    let match_score = 2i32;
    let mismatch_score = -1i32;
    let gap_penalty = -2i32;

    // DP Matrix
    let mut dp = vec![vec![0i32; m + 1]; n + 1];

    for i in 0..=n {
        dp[i][0] = (i as i32) * gap_penalty;
    }
    for j in 0..=m {
        dp[0][j] = (j as i32) * gap_penalty;
    }

    for i in 1..=n {
        for j in 1..=m {
            let score_diag = dp[i - 1][j - 1]
                + if chars_a[i - 1] == chars_b[j - 1] {
                    match_score
                } else {
                    mismatch_score
                };
            let score_up = dp[i - 1][j] + gap_penalty;
            let score_left = dp[i][j - 1] + gap_penalty;
            dp[i][j] = score_diag.max(score_up).max(score_left);
        }
    }

    // Traceback
    let mut alignment = Vec::new();
    let mut i = n;
    let mut j = m;

    while i > 0 || j > 0 {
        if i > 0
            && j > 0
            && dp[i][j]
                == dp[i - 1][j - 1]
                    + if chars_a[i - 1] == chars_b[j - 1] {
                        match_score
                    } else {
                        mismatch_score
                    }
        {
            alignment.push((Some(i - 1), Some(j - 1)));
            i -= 1;
            j -= 1;
        } else if i > 0 && dp[i][j] == dp[i - 1][j] + gap_penalty {
            alignment.push((Some(i - 1), None));
            i -= 1;
        } else {
            alignment.push((None, Some(j - 1)));
            j -= 1;
        }
    }

    alignment.reverse();
    alignment
}

/// Pairs corresponding C-alpha (or nucleic P) coordinates between two structures using sequence alignment.
pub fn pair_ca_coordinates(s1: &Structure, s2: &Structure) -> Vec<(Vec3, Vec3)> {
    let mut pairs = Vec::new();

    for c1 in s1.chains() {
        if let Some(c2) = s2.chains().iter().find(|c| c.id == c1.id).or_else(|| s2.chains().first()) {
            let seq1: String = c1.residues.iter().map(|r| residue_to_one_letter(&r.name)).collect();
            let seq2: String = c2.residues.iter().map(|r| residue_to_one_letter(&r.name)).collect();

            let alignment = needleman_wunsch(&seq1, &seq2);
            let atoms1 = s1.atoms();
            let atoms2 = s2.atoms();

            for (opt_i, opt_j) in alignment {
                if let (Some(i), Some(j)) = (opt_i, opt_j) {
                    let r1 = &c1.residues[i];
                    let r2 = &c2.residues[j];

                    let ca1 = r1.ca_atom(atoms1).or_else(|| r1.atom_indices.first().map(|&idx| &atoms1[idx]));
                    let ca2 = r2.ca_atom(atoms2).or_else(|| r2.atom_indices.first().map(|&idx| &atoms2[idx]));

                    if let (Some(a1), Some(a2)) = (ca1, ca2) {
                        pairs.push((a1.pos, a2.pos));
                    }
                }
            }
        }
    }

    pairs
}

/// Superimposes `mobile` onto `target` in place using the Kabsch algorithm, minimizing RMSD.
pub fn superimpose_structures(
    mobile: &mut Structure,
    target: &Structure,
) -> Option<SuperpositionResult> {
    let pairs = pair_ca_coordinates(mobile, target);
    if pairs.len() < 3 {
        return None;
    }

    let p: Vec<Vec3> = pairs.iter().map(|(p1, _)| *p1).collect();
    let q: Vec<Vec3> = pairs.iter().map(|(_, p2)| *p2).collect();

    let kabsch = kabsch_align(&p, &q)?;

    // Compute per-residue RMSD
    let mut per_residue_rmsd = Vec::with_capacity(pairs.len());
    for &(p1, p2) in &pairs {
        let transformed = kabsch.transform_point(p1);
        per_residue_rmsd.push(transformed.distance(&p2));
    }

    // Apply transformation in-place to all atoms in mobile active model
    for atom in mobile.atoms_mut() {
        atom.pos = kabsch.transform_point(atom.pos);
    }

    Some(SuperpositionResult {
        kabsch,
        per_residue_rmsd,
        aligned_pairs: pairs.len(),
    })
}

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
    /// `None` for proline (no amide hydrogen -> cannot donate).
    h_pos: Option<Vec3>,
}

/// Upper bound on the acceptor-C to donor-N distance worth evaluating. A bond
/// passing the -0.5 kcal/mol cutoff requires an O...H contact well under
/// 3.5 A, which bounds C...N far below this; used only to prune candidates.
const DONOR_SEARCH_RADIUS: f32 = 8.0;

/// Minimal local hash grid mapping backbone-residue indices by their donor-N
/// position, for radius queries. (`SpatialGrid` is atom-typed; H-bond
/// candidates are residues.)
struct DonorGrid {
    cell_size: f32,
    cells: std::collections::HashMap<(i32, i32, i32), Vec<u32>>,
}

impl DonorGrid {
    fn new(cell_size: f32) -> Self {
        Self {
            cell_size,
            cells: std::collections::HashMap::new(),
        }
    }

    fn key(pos: Vec3, inv: f32) -> (i32, i32, i32) {
        (
            (pos.x * inv).floor() as i32,
            (pos.y * inv).floor() as i32,
            (pos.z * inv).floor() as i32,
        )
    }

    fn insert(&mut self, pos: Vec3, idx: u32) {
        let key = Self::key(pos, 1.0 / self.cell_size);
        self.cells.entry(key).or_default().push(idx);
    }

    /// All inserted indices stored in cells overlapping the query sphere.
    /// Exact per-candidate distance checks are the caller's job (it knows
    /// which residue position to measure against).
    fn query(&self, center: Vec3, radius: f32) -> impl Iterator<Item = u32> + '_ {
        let inv = 1.0 / self.cell_size;
        let (cx, cy, cz) = (
            (center.x * inv).floor() as i32,
            (center.y * inv).floor() as i32,
            (center.z * inv).floor() as i32,
        );
        let r_cells = (radius / self.cell_size).ceil() as i32;
        let r_cells = r_cells.max(0) as u32;
        self.cells
            .iter()
            .filter(move |&(&(gx, gy, gz), _)| {
                gx.abs_diff(cx) <= r_cells
                    && gy.abs_diff(cy) <= r_cells
                    && gz.abs_diff(cz) <= r_cells
            })
            .flat_map(|(_, idxs)| idxs.iter().copied())
    }
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

            // Approximate backbone H position if missing. Proline is excluded:
            // its nitrogen is part of the pyrrolidine ring and carries no
            // amide hydrogen, so it can never donate (DSSP rule).
            let is_proline = res.name.trim().eq_ignore_ascii_case("PRO");
            let h_pos = if is_proline {
                None
            } else {
                explicit_h.or_else(|| {
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
                })
            };

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

    // --- H-bond detection (sparse + spatially pruned) -----------------------
    //
    // `hbond_acceptors[i]` holds the sorted donor indices j for which
    // E(C=O_i, N-H_j) < cutoff. Donor candidates come from a local hash grid
    // over N positions, so cost is O(n * neighbors) time and O(n + bonds)
    // memory instead of the previous O(n^2) matrix (which needed ~10 GB at
    // 100k residues). Pairs may span chains -- inter-chain beta-sheets are
    // real chemistry -- but same-chain near-sequence pairs are skipped.
    let mut hbond_acceptors: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut donors = DonorGrid::new(DONOR_SEARCH_RADIUS);
    for (j, bb) in backbone.iter().enumerate() {
        if let (Some(n_pos), Some(_h)) = (bb.n_pos, bb.h_pos) {
            donors.insert(n_pos, j as u32);
        }
    }

    for i in 0..n {
        let (Some(acc_c), Some(acc_o)) = (backbone[i].c_pos, backbone[i].o_pos) else {
            continue;
        };
        for j in donors.query(acc_c, DONOR_SEARCH_RADIUS) {
            let j = j as usize;
            let don_n = match backbone[j].n_pos {
                Some(p) => p,
                None => continue,
            };
            let don_h = match backbone[j].h_pos {
                Some(p) => p,
                None => continue, // proline or missing hydrogen: cannot donate
            };
            if backbone[i].chain_idx == backbone[j].chain_idx && (i as isize - j as isize).abs() < 2
            {
                continue; // Same or adjacent residue in the same chain.
            }
            let energy = calculate_dssp_hbond_energy(acc_c, acc_o, don_n, don_h);
            if energy < DSSP_HBOND_CUTOFF {
                hbond_acceptors[i].push(j as u32);
            }
        }
    }
    for row in &mut hbond_acceptors {
        row.sort_unstable();
    }

    let has_hbond = |acc: usize, don: usize| -> bool {
        acc < n && don < n && hbond_acceptors[acc].binary_search(&(don as u32)).is_ok()
    };
    let same_chain =
        |a: usize, b: usize| -> bool { backbone[a].chain_idx == backbone[b].chain_idx };

    let mut ss = vec![SecondaryStructure::Coil; n];

    // 1+2. Helices require two CONSECUTIVE n-turns within one chain (the DSSP
    // definition of an H/G structure). A lone n-turn is noise and no longer
    // assigns anything; cross-chain i -> i+n pairs are never turns.
    let turn =
        |i: usize, k: usize| -> bool { i + k < n && same_chain(i, i + k) && has_hbond(i, i + k) };

    // 4-turn alpha helices: overlapping turns at i and i+1 mark i+1..=i+4.
    for i in 0..n.saturating_sub(4) {
        if turn(i, 4) && turn(i + 1, 4) {
            for slot in &mut ss[(i + 1)..=(i + 4)] {
                *slot = SecondaryStructure::Helix;
            }
        }
    }

    // 3-turn 3_10 helices: same consecutive-turn rule, fills remaining coil.
    for i in 0..n.saturating_sub(3) {
        if turn(i, 3) && turn(i + 1, 3) {
            for slot in &mut ss[(i + 1)..=(i + 3)] {
                if *slot == SecondaryStructure::Coil {
                    *slot = SecondaryStructure::Helix;
                }
            }
        }
    }

    // 3. Beta bridges. A bridge predicate on outer pair (i, j) references at
    // most eight concrete H-bonds involving i, j and their chain-local
    // sequence neighbors, so instead of scanning all n^2 residue pairs we
    // enumerate candidate outer pairs from the H-bond list itself: every pair
    // that could satisfy a predicate is a shift-by-one of some observed bond.
    let mut candidates: std::collections::BTreeSet<(usize, usize)> =
        std::collections::BTreeSet::new();
    let mut push_candidate = |a: isize, b: isize| {
        if a >= 0 && b >= 0 {
            let (lo, hi) = if a <= b {
                (a as usize, b as usize)
            } else {
                (b as usize, a as usize)
            };
            // Shifts may push an endpoint past the last residue; such outer
            // pairs were unreachable under the original per-clause bounds
            // (e.g. `j + 1 < n`), so drop them here.
            if hi < n {
                candidates.insert((lo, hi));
            }
        }
    };
    for acc in 0..n {
        for &don_u in &hbond_acceptors[acc] {
            let don = don_u as usize;
            for &(x, y) in &[(acc, don), (don, acc)] {
                push_candidate(x as isize, y as isize);
                push_candidate(x as isize + 1, y as isize - 1);
                push_candidate(x as isize, y as isize - 1);
                push_candidate(x as isize + 1, y as isize);
            }
        }
    }

    for &(i, j) in &candidates {
        if j - i < 3 {
            continue;
        }
        // Sequence-neighbor steps must stay inside one chain; the bridging
        // H-bonds themselves may cross chains.
        let prev_ok_i = i > 0 && same_chain(i, i - 1);
        let next_ok_i = i + 1 < n && same_chain(i, i + 1);
        let prev_ok_j = j > 0 && same_chain(j, j - 1);
        let next_ok_j = j + 1 < n && same_chain(j, j + 1);

        // Antiparallel bridge: (i->j and j->i) or (i-1->j+1 and j-1->i+1)
        let ap_bridge = (has_hbond(i, j) && has_hbond(j, i))
            || (prev_ok_i
                && next_ok_j
                && prev_ok_j
                && next_ok_i
                && has_hbond(i - 1, j + 1)
                && has_hbond(j - 1, i + 1));

        // Parallel bridge: (i-1->j and j->i+1) or (j-1->i and i->j+1)
        let p_bridge = (prev_ok_i && next_ok_i && has_hbond(i - 1, j) && has_hbond(j, i + 1))
            || (prev_ok_j && next_ok_j && has_hbond(j - 1, i) && has_hbond(i, j + 1));

        if ap_bridge || p_bridge {
            ss[i] = SecondaryStructure::Sheet;
            ss[j] = SecondaryStructure::Sheet;
        }
    }

    // Apply assigned secondary structure back to structure residues.
    //
    // Precondition: `chains()` (the active view) and `chains_mut()` (the
    // deposited model) must address the same residues -- true because callers
    // run DSSP at parse time, before any biological-assembly activation.
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

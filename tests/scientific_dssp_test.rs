use termpdb::math::Vec3;
use termpdb::model::dssp::{assign_dssp, calculate_dssp_hbond_energy};
use termpdb::model::{Atom, Chain, Element, Residue, SecondaryStructure, Structure};

fn create_helical_structure(num_residues: usize) -> Structure {
    let mut structure = Structure::new("helix_peptide");
    let mut chain = Chain::new("A");
    let mut atom_counter = 0;

    let r_helix = 2.3f32;
    let pitch_per_res = 1.5f32;
    let angle_step = 100.0f32.to_radians();

    for i in 0..num_residues {
        let res_num = (i + 1) as i32;
        let mut res = Residue::new(res_num, "ALA", "A");

        let theta = (i as f32) * angle_step;
        let z = (i as f32) * pitch_per_res;

        // C-alpha
        let ca_pos = Vec3::new(r_helix * theta.cos(), r_helix * theta.sin(), z);
        // N
        let theta_n = theta - 0.35;
        let n_pos = Vec3::new(
            (r_helix - 0.4) * theta_n.cos(),
            (r_helix - 0.4) * theta_n.sin(),
            z - 0.6,
        );
        // C
        let theta_c = theta + 0.35;
        let c_pos = Vec3::new(
            (r_helix - 0.3) * theta_c.cos(),
            (r_helix - 0.3) * theta_c.sin(),
            z + 0.6,
        );
        // O pointing parallel along helix axis toward +Z (residue i+4)
        let o_pos = c_pos + Vec3::new(0.0, 0.0, 1.23);

        let elem_c = Element {
            symbol: "C",
            name: "Carbon",
            atomic_number: 6,
            covalent_radius: 0.76,
            vdw_radius: 1.7,
            cpk_color: (144, 144, 144),
        };
        let elem_n = Element {
            symbol: "N",
            name: "Nitrogen",
            atomic_number: 7,
            covalent_radius: 0.71,
            vdw_radius: 1.55,
            cpk_color: (48, 80, 248),
        };
        let elem_o = Element {
            symbol: "O",
            name: "Oxygen",
            atomic_number: 8,
            covalent_radius: 0.66,
            vdw_radius: 1.52,
            cpk_color: (255, 13, 13),
        };

        let atom_n_idx = atom_counter;
        structure.add_atom(Atom::new(
            atom_n_idx,
            atom_n_idx as i32 + 1,
            "N",
            elem_n,
            n_pos,
            20.0,
            "ALA",
            res_num,
            "A",
            false,
        ));
        res.atom_indices.push(atom_n_idx);
        atom_counter += 1;

        let atom_ca_idx = atom_counter;
        structure.add_atom(Atom::new(
            atom_ca_idx,
            atom_ca_idx as i32 + 1,
            "CA",
            elem_c,
            ca_pos,
            20.0,
            "ALA",
            res_num,
            "A",
            false,
        ));
        res.atom_indices.push(atom_ca_idx);
        atom_counter += 1;

        let atom_c_idx = atom_counter;
        structure.add_atom(Atom::new(
            atom_c_idx,
            atom_c_idx as i32 + 1,
            "C",
            elem_c,
            c_pos,
            20.0,
            "ALA",
            res_num,
            "A",
            false,
        ));
        res.atom_indices.push(atom_c_idx);
        atom_counter += 1;

        let atom_o_idx = atom_counter;
        structure.add_atom(Atom::new(
            atom_o_idx,
            atom_o_idx as i32 + 1,
            "O",
            elem_o,
            o_pos,
            20.0,
            "ALA",
            res_num,
            "A",
            false,
        ));
        res.atom_indices.push(atom_o_idx);
        atom_counter += 1;

        chain.residues.push(res);
    }

    structure.add_chain(chain);
    structure
}

#[test]
fn test_dssp_hbond_energy_ideal_geometry() {
    // Ideal peptide C=O and N-H backbone hydrogen bond (~2.0 A H..O distance)
    let c = Vec3::new(0.0, 0.0, 0.0);
    let o = Vec3::new(1.23, 0.0, 0.0);
    let h = Vec3::new(3.13, 0.0, 0.0);
    let n = Vec3::new(4.13, 0.0, 0.0);

    let energy = calculate_dssp_hbond_energy(c, o, n, h);
    assert!(
        energy < -0.5,
        "Backbone H-bond energy should be < -0.5 kcal/mol, got {energy}"
    );
    assert!(
        energy > -4.0,
        "Backbone H-bond energy should be physically bounded, got {energy}"
    );
}

#[test]
fn test_dssp_assignment_helical_peptide() {
    let mut structure = create_helical_structure(10);
    let count = assign_dssp(&mut structure);
    assert!(
        count > 0,
        "DSSP should detect helical residues in canonical helix, assigned {count}"
    );

    let chain = &structure.chains()[0];
    let helices = chain
        .residues
        .iter()
        .filter(|r| r.secondary_structure == SecondaryStructure::Helix)
        .count();
    assert!(helices > 0,);
}

// ---------------------------------------------------------------------------
// Synthetic-geometry regression tests for the DSSP assignment rules.
//
// Bond unit template (colinear along +x, from the ideal-energy test):
//   acceptor: C = x0,        O = x0 + 1.23  (N/H parked off-axis, inert)
//   donor:    H = x0 + 3.13, N = x0 + 4.13 (C/O parked off-axis, inert)
// which yields E ~= -2.9 kcal/mol -- a solid H-bond from the acceptor
// residue to the donor residue.
// ---------------------------------------------------------------------------

fn elem(symbol: &'static str, z: u8, cov: f32, vdw: f32, cpk: (u8, u8, u8)) -> Element {
    Element {
        symbol,
        name: symbol,
        atomic_number: z,
        covalent_radius: cov,
        vdw_radius: vdw,
        cpk_color: cpk,
    }
}

fn elem_c() -> Element {
    elem("C", 6, 0.76, 1.70, (144, 144, 144))
}
fn elem_n() -> Element {
    elem("N", 7, 0.71, 1.55, (48, 80, 248))
}
fn elem_o() -> Element {
    elem("O", 8, 0.66, 1.52, (255, 13, 13))
}
fn elem_h() -> Element {
    elem("H", 1, 0.31, 1.20, (255, 255, 255))
}

/// Adds a residue with explicit N/C/O/H backbone atoms, in that order.
fn add_backbone_residue(
    structure: &mut Structure,
    chain: &mut Chain,
    chain_id: &str,
    seq: i32,
    name: &str,
    atoms: [Vec3; 4],
) {
    let [n, c, o, h] = atoms;
    let mut res = Residue::new(seq, name, chain_id);
    let next_idx = structure.atom_count();
    for (offset, (atom_name, pos, el)) in [
        ("N", n, elem_n()),
        ("C", c, elem_c()),
        ("O", o, elem_o()),
        ("H", h, elem_h()),
    ]
    .into_iter()
    .enumerate()
    {
        let idx = structure.add_atom(Atom::new(
            next_idx + offset,
            (next_idx + offset) as i32 + 1,
            atom_name,
            el,
            pos,
            0.0,
            name,
            seq,
            chain_id,
            false,
        ));
        res.atom_indices.push(idx);
    }
    chain.residues.push(res);
}

/// Adds an inert CA-only residue far off in space; it can neither donate nor
/// accept, keeping global residue indexing predictable without side effects.
fn add_inert_filler(structure: &mut Structure, chain: &mut Chain, chain_id: &str, seq: i32) {
    let idx = structure.add_atom(Atom::new(
        structure.atom_count(),
        structure.atom_count() as i32 + 1,
        "CA",
        elem_c(),
        Vec3::new(seq as f32 * 3.0, 0.0, 50.0),
        0.0,
        "GLY",
        seq,
        chain_id,
        false,
    ));
    let mut res = Residue::new(seq, "GLY", chain_id);
    res.atom_indices.push(idx);
    chain.residues.push(res);
}

/// One ideal H-bond unit translated along x by `x0`:
/// `(acceptor (n, c, o, h), donor (n, c, o, h))` backbone positions.
/// The unused atoms of each partner are parked ~7 A off-axis so they cannot
/// form stray bonds.
const fn unit(x0: f32) -> [(Vec3, Vec3, Vec3, Vec3); 2] {
    [
        // Acceptor: on-axis C=O; N and H parked at y=z=5.
        (
            Vec3::new(x0, 5.0, 5.0),
            Vec3::new(x0, 0.0, 0.0),
            Vec3::new(x0 + 1.23, 0.0, 0.0),
            Vec3::new(x0 + 1.0, 5.0, 5.0),
        ),
        // Donor: on-axis H<N pointing back at the acceptor O; C/O parked.
        (
            Vec3::new(x0 + 4.13, 0.0, 0.0),
            Vec3::new(x0 + 7.0, 5.0, 5.0),
            Vec3::new(x0 + 7.6, 5.4, 5.0),
            Vec3::new(x0 + 3.13, 0.0, 0.0),
        ),
    ]
}

fn count_ss(structure: &Structure, ss: SecondaryStructure) -> usize {
    structure
        .chains()
        .iter()
        .flat_map(|c| c.residues.iter())
        .filter(|r| r.secondary_structure == ss)
        .count()
}

#[test]
fn test_dssp_single_isolated_nturn_assigns_nothing() {
    // One lone n-turn (gap 4) between residues 2 and 6 with perfect geometry;
    // all other residues are inert. DSSP requires TWO consecutive turns to
    // call a helix, so nothing may be assigned.
    let mut structure = Structure::new("single_turn");
    let mut chain = Chain::new("A");
    add_inert_filler(&mut structure, &mut chain, "A", 1); // idx 0
    add_inert_filler(&mut structure, &mut chain, "A", 2); // idx 1
    let [acc, _don] = unit(0.0);
    add_backbone_residue(
        &mut structure,
        &mut chain,
        "A",
        3,
        "ALA",
        [acc.0, acc.1, acc.2, acc.3],
    ); // idx 2
    add_inert_filler(&mut structure, &mut chain, "A", 4); // idx 3
    add_inert_filler(&mut structure, &mut chain, "A", 5); // idx 4
    add_inert_filler(&mut structure, &mut chain, "A", 6); // idx 5
    let [_acc, don] = unit(0.0);
    add_backbone_residue(
        &mut structure,
        &mut chain,
        "A",
        7,
        "ALA",
        [don.0, don.1, don.2, don.3],
    ); // idx 6
    add_inert_filler(&mut structure, &mut chain, "A", 8); // idx 7
    structure.add_chain(chain);

    assert_eq!(
        assign_dssp(&mut structure),
        0,
        "a lone n-turn must not assign anything"
    );
    assert_eq!(count_ss(&structure, SecondaryStructure::Helix), 0);
    assert_eq!(count_ss(&structure, SecondaryStructure::Sheet), 0);
}

/// Two consecutive ideal turns: acceptors at idx 2/3, donors at idx 6/7.
/// `acceptor_name` / `donor_name` control the residue identities.
fn build_two_turn_structure(acceptor_name: &str, donor_name: &str) -> Structure {
    let mut structure = Structure::new("two_turn");
    let mut chain = Chain::new("A");
    add_inert_filler(&mut structure, &mut chain, "A", 1); // idx 0
    add_inert_filler(&mut structure, &mut chain, "A", 2); // idx 1
    let [acc1, _] = unit(0.0);
    add_backbone_residue(
        &mut structure,
        &mut chain,
        "A",
        3,
        acceptor_name,
        [acc1.0, acc1.1, acc1.2, acc1.3],
    ); // idx 2
    let [acc2, _] = unit(30.0);
    add_backbone_residue(
        &mut structure,
        &mut chain,
        "A",
        4,
        acceptor_name,
        [acc2.0, acc2.1, acc2.2, acc2.3],
    ); // idx 3
    add_inert_filler(&mut structure, &mut chain, "A", 5); // idx 4
    add_inert_filler(&mut structure, &mut chain, "A", 6); // idx 5
    let [_a, don1] = unit(0.0);
    add_backbone_residue(
        &mut structure,
        &mut chain,
        "A",
        7,
        donor_name,
        [don1.0, don1.1, don1.2, don1.3],
    ); // idx 6
    let [_a, don2] = unit(30.0);
    add_backbone_residue(
        &mut structure,
        &mut chain,
        "A",
        8,
        donor_name,
        [don2.0, don2.1, don2.2, don2.3],
    ); // idx 7
    structure.add_chain(chain);
    structure
}

#[test]
fn test_dssp_proline_never_donates() {
    // Two consecutive turns with PROLINE donors given perfectly placed
    // explicit hydrogens: proline's ring nitrogen has no amide H and must not
    // donate no matter what the file claims.
    let mut structure = build_two_turn_structure("ALA", "PRO");
    assert_eq!(
        assign_dssp(&mut structure),
        0,
        "proline donors must not produce H-bonds"
    );
    assert_eq!(count_ss(&structure, SecondaryStructure::Helix), 0);
}

#[test]
fn test_dssp_proline_as_acceptor_still_participates() {
    // Mirror of the donation test: proline ACCEPTORS keep their C=O, so two
    // consecutive ALA-donor turns onto PRO acceptors must still mark a helix.
    let mut structure = build_two_turn_structure("PRO", "ALA");
    // DSSP marks {i+1..=i+4} for the consecutive-turn pair starting at i,
    // i.e. residues 3, 4, 5, 6.
    let assigned = assign_dssp(&mut structure);
    assert_eq!(assigned, 4, "ALA->PRO turns are valid helical bonds");
    assert_eq!(count_ss(&structure, SecondaryStructure::Helix), 4);
}

#[test]
fn test_dssp_helix_never_spans_chain_break() {
    // Turns (g1->g5) and (g2->g6) have perfect geometry but cross the A|B
    // chain boundary and must not fire; a positive control with identical
    // geometry inside chain C must still mark its four helix residues.
    let mut structure = Structure::new("chain_break");

    // Chain A: filler, acceptors for the cross-chain pair (globals 1, 2).
    let mut chain_a = Chain::new("A");
    add_inert_filler(&mut structure, &mut chain_a, "A", 1); // g0
    let [acc1, _] = unit(0.0);
    add_backbone_residue(
        &mut structure,
        &mut chain_a,
        "A",
        2,
        "ALA",
        [acc1.0, acc1.1, acc1.2, acc1.3],
    ); // g1
    let [acc2, _] = unit(30.0);
    add_backbone_residue(
        &mut structure,
        &mut chain_a,
        "A",
        3,
        "ALA",
        [acc2.0, acc2.1, acc2.2, acc2.3],
    ); // g2
    structure.add_chain(chain_a);

    // Chain B: fillers then the matching donors (globals 5, 6).
    let mut chain_b = Chain::new("B");
    add_inert_filler(&mut structure, &mut chain_b, "B", 1); // g3
    add_inert_filler(&mut structure, &mut chain_b, "B", 2); // g4
    let [_a, don1] = unit(0.0);
    add_backbone_residue(
        &mut structure,
        &mut chain_b,
        "B",
        3,
        "ALA",
        [don1.0, don1.1, don1.2, don1.3],
    ); // g5
    let [_a, don2] = unit(30.0);
    add_backbone_residue(
        &mut structure,
        &mut chain_b,
        "B",
        4,
        "ALA",
        [don2.0, don2.1, don2.2, don2.3],
    ); // g6
    structure.add_chain(chain_b);

    // Chain C (control): intra-chain consecutive turns (g8->g12, g9->g13).
    let mut chain_c = Chain::new("C");
    add_inert_filler(&mut structure, &mut chain_c, "C", 1); // g7
    let [c_acc1, _] = unit(100.0);
    add_backbone_residue(
        &mut structure,
        &mut chain_c,
        "C",
        2,
        "ALA",
        [c_acc1.0, c_acc1.1, c_acc1.2, c_acc1.3],
    ); // g8
    let [c_acc2, _] = unit(130.0);
    add_backbone_residue(
        &mut structure,
        &mut chain_c,
        "C",
        3,
        "ALA",
        [c_acc2.0, c_acc2.1, c_acc2.2, c_acc2.3],
    ); // g9
    add_inert_filler(&mut structure, &mut chain_c, "C", 4); // g10
    add_inert_filler(&mut structure, &mut chain_c, "C", 5); // g11
    let [_a, c_don1] = unit(100.0);
    add_backbone_residue(
        &mut structure,
        &mut chain_c,
        "C",
        6,
        "ALA",
        [c_don1.0, c_don1.1, c_don1.2, c_don1.3],
    ); // g12
    let [_a, c_don2] = unit(130.0);
    add_backbone_residue(
        &mut structure,
        &mut chain_c,
        "C",
        7,
        "ALA",
        [c_don2.0, c_don2.1, c_don2.2, c_don2.3],
    ); // g13
    structure.add_chain(chain_c);

    let assigned = assign_dssp(&mut structure);
    // DSSP marks {i+1..=i+4} for the consecutive-turn pair starting at g8/g9:
    // residues g9..=g12 (four).
    assert_eq!(count_ss(&structure, SecondaryStructure::Helix), 4);
    assert_eq!(assigned, 4);
    // Chains A and B (the cross-chain pair) stay coil.
    let ab_coil = structure.chains()[0]
        .residues
        .iter()
        .chain(structure.chains()[1].residues.iter())
        .all(|r| r.secondary_structure == SecondaryStructure::Coil);
    assert!(ab_coil, "cross-chain turn pairs must not label chains A/B");
}

/// Absolute coordinates of a reciprocal (antiparallel) bridge pair:
/// each residue's C=O faces the other's N-H, E ~= -2.58 both directions.
const BRIDGE_X: (Vec3, Vec3, Vec3, Vec3) = (
    Vec3::new(3.0, 0.0, 0.0),   // N
    Vec3::new(-3.23, 0.0, 0.0), // C
    Vec3::new(-2.0, 0.0, 0.0),  // O
    Vec3::new(4.0, 0.0, 0.0),   // H
);
const BRIDGE_Y: (Vec3, Vec3, Vec3, Vec3) = (
    Vec3::new(1.0, 0.0, 0.0),  // N
    Vec3::new(7.23, 0.0, 0.0), // C
    Vec3::new(6.0, 0.0, 0.0),  // O
    Vec3::new(0.0, 0.0, 0.0),  // H
);

#[test]
fn test_dssp_antiparallel_bridge_assigns_sheet() {
    let mut structure = Structure::new("sheet");
    let mut chain = Chain::new("A");
    add_inert_filler(&mut structure, &mut chain, "A", 1); // idx 0
    add_inert_filler(&mut structure, &mut chain, "A", 2); // idx 1
    let (xn, xc, xo, xh) = BRIDGE_X;
    add_backbone_residue(&mut structure, &mut chain, "A", 3, "VAL", [xn, xc, xo, xh]); // idx 2 = X
    for seq in 4..=9 {
        add_inert_filler(&mut structure, &mut chain, "A", seq); // idx 3..8
    }
    let (yn, yc, yo, yh) = BRIDGE_Y;
    add_backbone_residue(&mut structure, &mut chain, "A", 10, "VAL", [yn, yc, yo, yh]); // idx 9 = Y
    add_inert_filler(&mut structure, &mut chain, "A", 11); // idx 10
    structure.add_chain(chain);

    let assigned = assign_dssp(&mut structure);
    assert_eq!(assigned, 2, "exactly the two bridging residues");
    assert_eq!(count_ss(&structure, SecondaryStructure::Sheet), 2);
}

#[test]
fn test_dssp_allows_inter_chain_sheet_bridge() {
    // Same mutual geometry split across chains A and B: inter-chain beta
    // sheets are real chemistry and must still be detected. Fillers in chain
    // A keep the global index gap >= 3.
    let mut structure = Structure::new("interchain_sheet");
    let mut chain_a = Chain::new("A");
    let (xn, xc, xo, xh) = BRIDGE_X;
    add_backbone_residue(
        &mut structure,
        &mut chain_a,
        "A",
        1,
        "VAL",
        [xn, xc, xo, xh],
    ); // g0 = X
    for seq in 2..=4 {
        add_inert_filler(&mut structure, &mut chain_a, "A", seq); // g1..g3
    }
    structure.add_chain(chain_a);

    let mut chain_b = Chain::new("B");
    let (yn, yc, yo, yh) = BRIDGE_Y;
    add_backbone_residue(
        &mut structure,
        &mut chain_b,
        "B",
        1,
        "VAL",
        [yn, yc, yo, yh],
    ); // g4 = Y
    structure.add_chain(chain_b);

    let assigned = assign_dssp(&mut structure);
    assert_eq!(assigned, 2);
    assert_eq!(count_ss(&structure, SecondaryStructure::Sheet), 2);
}

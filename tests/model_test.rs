use termpdb::math::Vec3;
use termpdb::model::{
    Atom, Bond, BondDetector, BondOrder, Chain, Element, Model, Residue, SecondaryStructure,
    Structure, element_by_atomic_number, element_by_symbol,
};

#[test]
fn test_element_table_known_elements() {
    let carbon = element_by_symbol("C");
    assert_eq!(carbon.symbol, "C");
    assert_eq!(carbon.atomic_number, 6);
    assert!(carbon.covalent_radius > 0.7 && carbon.covalent_radius < 0.85);
    assert!(carbon.vdw_radius > 1.5 && carbon.vdw_radius < 1.9);

    let nitrogen = element_by_symbol("N");
    assert_eq!(nitrogen.symbol, "N");
    assert_eq!(nitrogen.atomic_number, 7);

    let oxygen = element_by_symbol("O");
    assert_eq!(oxygen.symbol, "O");
    assert_eq!(oxygen.atomic_number, 8);

    let hydrogen = element_by_symbol("H");
    assert_eq!(hydrogen.symbol, "H");
    assert_eq!(hydrogen.atomic_number, 1);

    let sulfur = element_by_symbol("S");
    assert_eq!(sulfur.symbol, "S");
    assert_eq!(sulfur.atomic_number, 16);

    let phosphorus = element_by_symbol("P");
    assert_eq!(phosphorus.symbol, "P");
    assert_eq!(phosphorus.atomic_number, 15);

    let iron = element_by_symbol("Fe");
    assert_eq!(iron.symbol, "Fe");
    assert_eq!(iron.atomic_number, 26);

    let zinc = element_by_symbol("Zn");
    assert_eq!(zinc.symbol, "Zn");
    assert_eq!(zinc.atomic_number, 30);
}

#[test]
fn test_element_case_insensitivity_and_trim() {
    let ca1 = element_by_symbol("ca");
    let ca2 = element_by_symbol("CA");
    let ca3 = element_by_symbol(" Ca ");
    assert_eq!(ca1.symbol, "Ca");
    assert_eq!(ca2.symbol, "Ca");
    assert_eq!(ca3.symbol, "Ca");
    assert_eq!(ca1.atomic_number, 20);

    let mg = element_by_symbol("mg");
    assert_eq!(mg.symbol, "Mg");
    assert_eq!(mg.atomic_number, 12);
}

#[test]
fn test_element_by_atomic_number() {
    let h = element_by_atomic_number(1);
    assert_eq!(h.symbol, "H");

    let c = element_by_atomic_number(6);
    assert_eq!(c.symbol, "C");

    let fe = element_by_atomic_number(26);
    assert_eq!(fe.symbol, "Fe");

    let unknown = element_by_atomic_number(250);
    assert_eq!(unknown.atomic_number, 0);
}

#[test]
fn test_element_unknown_fallback() {
    let unk = element_by_symbol("XYZ_UNKNOWN");
    assert_eq!(unk.atomic_number, 0);
    assert!(unk.covalent_radius > 0.0);
    assert!(unk.vdw_radius > 0.0);
}

#[test]
fn test_atom_creation_and_properties() {
    let c_elem = element_by_symbol("C");
    let atom = Atom::new(
        0,
        101,
        "CA",
        c_elem,
        Vec3::new(1.0, 2.0, 3.0),
        25.5,
        "ALA",
        1,
        "A",
        false,
    );

    assert_eq!(atom.index, 0);
    assert_eq!(atom.serial, 101);
    assert_eq!(atom.name, "CA");
    assert_eq!(atom.element.symbol, "C");
    assert_eq!(atom.pos, Vec3::new(1.0, 2.0, 3.0));
    assert_eq!(atom.b_factor, 25.5);
    assert_eq!(atom.res_name, "ALA");
    assert_eq!(atom.res_seq, 1);
    assert_eq!(atom.chain_id, "A");
    assert!(!atom.is_hetatm);
    assert!(atom.is_c_alpha());
    assert!(atom.is_backbone());
    assert!(!atom.is_hydrogen());

    let h_elem = element_by_symbol("H");
    let h_atom = Atom::new(
        1,
        102,
        "H",
        h_elem,
        Vec3::new(1.5, 2.5, 3.5),
        10.0,
        "ALA",
        1,
        "A",
        false,
    );
    assert!(h_atom.is_hydrogen());
    assert!(!h_atom.is_c_alpha());

    let het_ca = Atom::new(
        2,
        103,
        "CA",
        element_by_symbol("Ca"),
        Vec3::new(10.0, 10.0, 10.0),
        30.0,
        "CA",
        999,
        "A",
        true,
    );
    // HETATM calcium is not a protein C-alpha
    assert!(!het_ca.is_c_alpha());
}

#[test]
fn test_residue_classification_and_hydrophobicity() {
    let mut ala = Residue::new(1, "ALA", "A");
    ala.atom_indices = vec![0, 1, 2, 3, 4];
    ala.secondary_structure = SecondaryStructure::Helix;

    assert!(ala.is_amino_acid());
    assert!(!ala.is_nucleic());
    assert!(!ala.is_water());
    assert_eq!(ala.one_letter_code(), 'A');
    assert!(ala.hydrophobicity_score() > 0.0);

    let arg = Residue::new(2, "ARG", "A");
    assert!(arg.is_amino_acid());
    assert_eq!(arg.one_letter_code(), 'R');
    assert!(arg.hydrophobicity_score() < 0.0);

    let rna_a = Residue::new(1, "A", "B");
    assert!(!rna_a.is_amino_acid());
    assert!(rna_a.is_nucleic());

    let dna_da = Residue::new(2, "DA", "B");
    assert!(dna_da.is_nucleic());

    let hoh = Residue::new(100, "HOH", "W");
    assert!(hoh.is_water());
    assert!(!hoh.is_amino_acid());
    assert!(!hoh.is_nucleic());
}

#[test]
fn test_residue_ca_atom_lookup() {
    let c_elem = element_by_symbol("C");
    let n_elem = element_by_symbol("N");
    let atoms = vec![
        Atom::new(
            0,
            1,
            "N",
            n_elem,
            Vec3::new(0.0, 0.0, 0.0),
            10.0,
            "ALA",
            1,
            "A",
            false,
        ),
        Atom::new(
            1,
            2,
            "CA",
            c_elem,
            Vec3::new(1.4, 0.0, 0.0),
            10.0,
            "ALA",
            1,
            "A",
            false,
        ),
        Atom::new(
            2,
            3,
            "C",
            c_elem,
            Vec3::new(2.0, 1.2, 0.0),
            10.0,
            "ALA",
            1,
            "A",
            false,
        ),
    ];

    let mut res = Residue::new(1, "ALA", "A");
    res.atom_indices = vec![0, 1, 2];

    let ca = res.ca_atom(&atoms);
    assert!(ca.is_some());
    assert_eq!(ca.unwrap().index, 1);
    assert_eq!(ca.unwrap().name, "CA");
    assert_eq!(res.ca_atom_index(&atoms), Some(1));
}

#[test]
fn test_chain_trace_and_residue_access() {
    let c_elem = element_by_symbol("C");
    let atoms = vec![
        Atom::new(
            0,
            1,
            "CA",
            c_elem,
            Vec3::new(0.0, 0.0, 0.0),
            10.0,
            "ALA",
            1,
            "A",
            false,
        ),
        Atom::new(
            1,
            2,
            "CB",
            c_elem,
            Vec3::new(0.5, 1.0, 0.0),
            10.0,
            "ALA",
            1,
            "A",
            false,
        ),
        Atom::new(
            2,
            3,
            "CA",
            c_elem,
            Vec3::new(3.8, 0.0, 0.0),
            12.0,
            "GLY",
            2,
            "A",
            false,
        ),
        Atom::new(
            3,
            4,
            "CA",
            c_elem,
            Vec3::new(7.6, 0.0, 0.0),
            15.0,
            "VAL",
            3,
            "A",
            false,
        ),
    ];

    let mut chain = Chain::new("A");
    let mut res1 = Residue::new(1, "ALA", "A");
    res1.atom_indices = vec![0, 1];
    let mut res2 = Residue::new(2, "GLY", "A");
    res2.atom_indices = vec![2];
    let mut res3 = Residue::new(3, "VAL", "A");
    res3.atom_indices = vec![3];

    chain.residues.push(res1);
    chain.residues.push(res2);
    chain.residues.push(res3);

    let ca_trace = chain.ca_atoms(&atoms);
    assert_eq!(ca_trace.len(), 3);
    assert_eq!(ca_trace[0].index, 0);
    assert_eq!(ca_trace[1].index, 2);
    assert_eq!(ca_trace[2].index, 3);

    let ca_indices = chain.ca_atom_indices(&atoms);
    assert_eq!(ca_indices, vec![0, 2, 3]);

    assert!(chain.get_residue(2).is_some());
    assert_eq!(chain.get_residue(2).unwrap().name, "GLY");
    assert!(chain.get_residue(999).is_none());
}

#[test]
fn test_bond_detector_water_molecule() {
    let o_elem = element_by_symbol("O");
    let h_elem = element_by_symbol("H");

    // Standard water geometry: O at origin, H1 and H2 at ~0.96 Å, angle 104.5°
    let o = Atom::new(
        0,
        1,
        "O",
        o_elem,
        Vec3::new(0.0, 0.0, 0.0),
        0.0,
        "HOH",
        1,
        "A",
        true,
    );
    let h1 = Atom::new(
        1,
        2,
        "H1",
        h_elem,
        Vec3::new(0.757, 0.586, 0.0),
        0.0,
        "HOH",
        1,
        "A",
        true,
    );
    let h2 = Atom::new(
        2,
        3,
        "H2",
        h_elem,
        Vec3::new(-0.757, 0.586, 0.0),
        0.0,
        "HOH",
        1,
        "A",
        true,
    );

    let atoms = vec![o, h1, h2];
    let bonds = BondDetector::detect_bonds(&atoms);

    // Should detect O-H1 and O-H2, but NOT H1-H2 (distance is ~1.514 Å > 0.31+0.31+0.45=1.07 Å)
    assert_eq!(bonds.len(), 2);
    let has_o_h1 = bonds
        .iter()
        .any(|b| (b.atom1_idx == 0 && b.atom2_idx == 1) || (b.atom1_idx == 1 && b.atom2_idx == 0));
    let has_o_h2 = bonds
        .iter()
        .any(|b| (b.atom1_idx == 0 && b.atom2_idx == 2) || (b.atom1_idx == 2 && b.atom2_idx == 0));
    let has_h1_h2 = bonds
        .iter()
        .any(|b| (b.atom1_idx == 1 && b.atom2_idx == 2) || (b.atom1_idx == 2 && b.atom2_idx == 1));

    assert!(has_o_h1, "Expected O-H1 bond");
    assert!(has_o_h2, "Expected O-H2 bond");
    assert!(!has_h1_h2, "H1-H2 should not be bonded");
}

#[test]
fn test_bond_detector_no_self_or_coincident_bonds() {
    let c_elem = element_by_symbol("C");
    let c1 = Atom::new(
        0,
        1,
        "C1",
        c_elem,
        Vec3::new(0.0, 0.0, 0.0),
        0.0,
        "MOL",
        1,
        "A",
        false,
    );
    let c2 = Atom::new(
        1,
        2,
        "C2",
        c_elem,
        Vec3::new(0.0, 0.0, 0.0),
        0.0,
        "MOL",
        1,
        "A",
        false,
    );

    let bonds = BondDetector::detect_bonds(&[c1, c2]);
    assert_eq!(
        bonds.len(),
        0,
        "Atoms with distance < 0.4 Å should not form a bond"
    );
}

#[test]
fn test_bond_detector_distant_atoms() {
    let c_elem = element_by_symbol("C");
    let c1 = Atom::new(
        0,
        1,
        "C1",
        c_elem,
        Vec3::new(0.0, 0.0, 0.0),
        0.0,
        "MOL",
        1,
        "A",
        false,
    );
    let c2 = Atom::new(
        1,
        2,
        "C2",
        c_elem,
        Vec3::new(10.0, 0.0, 0.0),
        0.0,
        "MOL",
        1,
        "A",
        false,
    );

    let bonds = BondDetector::detect_bonds(&[c1, c2]);
    assert_eq!(bonds.len(), 0);
}

#[test]
fn test_structure_centroid_and_normalization() {
    let c_elem = element_by_symbol("C");
    let mut structure = Structure::new("Test Tetra");

    structure.add_atom(Atom::new(
        0,
        1,
        "C1",
        c_elem,
        Vec3::new(10.0, 10.0, 10.0),
        0.0,
        "ALA",
        1,
        "A",
        false,
    ));
    structure.add_atom(Atom::new(
        1,
        2,
        "C2",
        c_elem,
        Vec3::new(12.0, 10.0, 10.0),
        0.0,
        "ALA",
        1,
        "A",
        false,
    ));
    structure.add_atom(Atom::new(
        2,
        3,
        "C3",
        c_elem,
        Vec3::new(10.0, 12.0, 10.0),
        0.0,
        "ALA",
        1,
        "A",
        false,
    ));
    structure.add_atom(Atom::new(
        3,
        4,
        "C4",
        c_elem,
        Vec3::new(10.0, 10.0, 12.0),
        0.0,
        "ALA",
        1,
        "A",
        false,
    ));

    let com = structure.center_of_mass();
    assert!((com.x - 10.5).abs() < 1e-5);
    assert!((com.y - 10.5).abs() < 1e-5);
    assert!((com.z - 10.5).abs() < 1e-5);

    let initial_radius = structure.bounding_sphere_radius();
    assert!(initial_radius > 0.0);

    structure.center_and_normalize();
    let new_com = structure.center_of_mass();
    assert!(new_com.x.abs() < 1e-5);
    assert!(new_com.y.abs() < 1e-5);
    assert!(new_com.z.abs() < 1e-5);

    let normalized_radius = structure.bounding_sphere_radius();
    assert!((normalized_radius - initial_radius).abs() < 1e-5);
}

#[test]
fn test_structure_build_bonds_and_counts() {
    let c_elem = element_by_symbol("C");
    let n_elem = element_by_symbol("N");

    let mut structure = Structure::new("Dipeptide");
    let mut chain = Chain::new("A");

    let mut res1 = Residue::new(1, "ALA", "A");
    res1.atom_indices = vec![0, 1];
    let mut res2 = Residue::new(2, "GLY", "A");
    res2.atom_indices = vec![2];

    chain.residues.push(res1);
    chain.residues.push(res2);
    structure.add_chain(chain);

    // C1 and C2 are 1.5 Å apart (bonded)
    structure.add_atom(Atom::new(
        0,
        1,
        "CA",
        c_elem,
        Vec3::new(0.0, 0.0, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    ));
    structure.add_atom(Atom::new(
        1,
        2,
        "C",
        c_elem,
        Vec3::new(1.5, 0.0, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    ));
    // N is 1.33 Å from C (peptide bond)
    structure.add_atom(Atom::new(
        2,
        3,
        "N",
        n_elem,
        Vec3::new(2.83, 0.0, 0.0),
        10.0,
        "GLY",
        2,
        "A",
        false,
    ));

    structure.build_bonds();

    assert_eq!(structure.atom_count(), 3);
    assert_eq!(structure.chain_count(), 1);
    assert_eq!(structure.residue_count(), 2);
    assert_eq!(structure.bonds().len(), 2);
    assert_eq!(structure.ca_atoms().len(), 1);
}

#[test]
fn test_bond_order_and_bond_helpers() {
    let b1 = Bond::new(0, 1, 1.0);
    assert_eq!(b1.atom1_idx, 0);
    assert_eq!(b1.atom2_idx, 1);
    assert_eq!(b1.order, 1.0);
    assert_eq!(b1.other(0), Some(1));
    assert_eq!(b1.other(1), Some(0));
    assert_eq!(b1.other(2), None);

    assert_eq!(BondOrder::Single.as_f32(), 1.0);
    assert_eq!(BondOrder::Double.as_f32(), 2.0);
    assert_eq!(BondOrder::Triple.as_f32(), 3.0);
    assert_eq!(BondOrder::Aromatic.as_f32(), 1.5);
    assert_eq!(BondOrder::from_f32(1.0), BondOrder::Single);
    assert_eq!(BondOrder::from_f32(2.0), BondOrder::Double);
    assert_eq!(BondOrder::from_f32(3.0), BondOrder::Triple);
    assert_eq!(BondOrder::from_f32(1.5), BondOrder::Aromatic);
    assert_eq!(BondOrder::from_f32(4.0), BondOrder::Other(4.0));
}

#[test]
fn test_structure_metadata_and_bfactor_range() {
    let mut structure = Structure::with_id("1CRN", "Crambin");
    assert_eq!(structure.id_code.as_deref(), Some("1CRN"));
    assert_eq!(structure.title, "Crambin");

    let c = element_by_symbol("C");
    structure.add_atom(Atom::new(
        0,
        1,
        "CA",
        c,
        Vec3::new(0.0, 0.0, 0.0),
        12.5,
        "THR",
        1,
        "A",
        false,
    ));
    structure.add_atom(Atom::new(
        1,
        2,
        "CB",
        c,
        Vec3::new(1.0, 1.0, 1.0),
        85.0,
        "THR",
        1,
        "A",
        false,
    ));
    structure.add_atom(Atom::new(
        2,
        3,
        "CG2",
        c,
        Vec3::new(2.0, 2.0, 2.0),
        45.0,
        "THR",
        1,
        "A",
        false,
    ));

    let (min_b, max_b) = structure.b_factor_range();
    assert!((min_b - 12.5).abs() < 1e-4);
    assert!((max_b - 85.0).abs() < 1e-4);
}

#[test]
fn test_secondary_structure_predicates() {
    let h = SecondaryStructure::Helix;
    let s = SecondaryStructure::Sheet;
    let c = SecondaryStructure::Coil;
    let def = SecondaryStructure::default();

    assert!(h.is_helix());
    assert!(!h.is_sheet());
    assert!(s.is_sheet());
    assert!(!s.is_coil());
    assert!(c.is_coil());
    assert_eq!(def, SecondaryStructure::Coil);
}

#[test]
fn test_element_default() {
    let def = Element::default();
    assert_eq!(def.atomic_number, 0);
    assert_eq!(def.symbol, "X");
}

#[test]
fn test_structure_set_models_activates_lowest_serial() {
    let mut m1 = Model::new(1);
    m1.atoms.push(Atom::new(
        0,
        1,
        "CA",
        element_by_symbol("C"),
        Vec3::new(0.0, 0.0, 0.0),
        0.0,
        "ALA",
        1,
        "A",
        false,
    ));
    let mut m4 = Model::new(4);
    m4.atoms.push(Atom::new(
        0,
        1,
        "CA",
        element_by_symbol("C"),
        Vec3::new(4.0, 0.0, 0.0),
        0.0,
        "ALA",
        1,
        "A",
        false,
    ));
    let mut m2 = Model::new(2);
    m2.atoms.push(Atom::new(
        0,
        1,
        "CA",
        element_by_symbol("C"),
        Vec3::new(2.0, 0.0, 0.0),
        0.0,
        "ALA",
        1,
        "A",
        false,
    ));

    let mut structure = Structure::new("ensemble");
    structure.set_models(vec![m4, m1, m2]);

    assert_eq!(structure.model_serials(), vec![1, 2, 4]);
    assert_eq!(structure.active_model_serial(), 1);
    assert_eq!(structure.max_model_serial(), 4);
    assert!((structure.atoms()[0].pos.x - 0.0).abs() < 1e-4);

    structure.next_model();
    assert_eq!(structure.active_model_serial(), 2);
    structure.next_model();
    assert_eq!(structure.active_model_serial(), 4);
    structure.next_model();
    assert_eq!(structure.active_model_serial(), 1);
    structure.prev_model();
    assert_eq!(structure.active_model_serial(), 4);
}

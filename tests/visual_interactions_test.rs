use termpdb::math::Vec3;
use termpdb::model::interactions::{InteractionKind, detect_interactions};
use termpdb::model::{Atom, Chain, Residue, Structure, element_by_symbol};

#[test]
fn test_disulfide_bond_detection() {
    let mut structure = Structure::new("disulfide_test");
    let mut chain = Chain::new("A");

    // Two CYS residues with SG atoms 2.05 A apart
    let mut r1 = Residue::new(1, "CYS", "A");
    let a1 = Atom::new(
        0,
        1,
        "SG",
        element_by_symbol("S"),
        Vec3::new(0.0, 0.0, 0.0),
        10.0,
        "CYS",
        1,
        "A",
        false,
    );
    r1.atom_indices.push(0);

    let mut r2 = Residue::new(5, "CYS", "A");
    let a2 = Atom::new(
        1,
        2,
        "SG",
        element_by_symbol("S"),
        Vec3::new(2.05, 0.0, 0.0),
        10.0,
        "CYS",
        5,
        "A",
        false,
    );
    r2.atom_indices.push(1);

    chain.residues.push(r1);
    chain.residues.push(r2);
    structure.add_chain(chain);
    structure.add_atom(a1);
    structure.add_atom(a2);

    let interactions = detect_interactions(&structure);
    assert_eq!(interactions.len(), 1);
    assert_eq!(interactions[0].kind, InteractionKind::Disulfide);
    assert!((interactions[0].distance - 2.05).abs() < 1e-3);
}

#[test]
fn test_hydrogen_bond_detection() {
    let mut structure = Structure::new("hbond_test");
    let mut chain = Chain::new("A");

    // Donor N atom at (0, 0, 0) in res 1, Acceptor O atom at (2.8, 0, 0) in res 5
    let mut r1 = Residue::new(1, "ALA", "A");
    let a1 = Atom::new(
        0,
        1,
        "N",
        element_by_symbol("N"),
        Vec3::new(0.0, 0.0, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    );
    r1.atom_indices.push(0);

    let mut r2 = Residue::new(5, "VAL", "A");
    let a2 = Atom::new(
        1,
        2,
        "O",
        element_by_symbol("O"),
        Vec3::new(2.8, 0.0, 0.0),
        10.0,
        "VAL",
        5,
        "A",
        false,
    );
    r2.atom_indices.push(1);

    chain.residues.push(r1);
    chain.residues.push(r2);
    structure.add_chain(chain);
    structure.add_atom(a1);
    structure.add_atom(a2);

    let interactions = detect_interactions(&structure);
    assert_eq!(interactions.len(), 1);
    assert_eq!(interactions[0].kind, InteractionKind::HydrogenBond);
    assert!((interactions[0].distance - 2.8).abs() < 1e-3);
}

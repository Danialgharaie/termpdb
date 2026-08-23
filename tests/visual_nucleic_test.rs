use termpdb::math::Vec3;
use termpdb::model::{Atom, Residue, element_by_symbol};
use termpdb::render::representations::nucleic::{base_color, build_base_slab, is_nucleic_residue};

#[test]
fn test_nucleic_residue_recognition_and_colors() {
    assert!(is_nucleic_residue("DA"));
    assert!(is_nucleic_residue("DT"));
    assert!(is_nucleic_residue("DG"));
    assert!(is_nucleic_residue("DC"));
    assert!(is_nucleic_residue("A"));
    assert!(is_nucleic_residue("U"));
    assert!(!is_nucleic_residue("ALA"));
    assert!(!is_nucleic_residue("CYS"));

    let col_a = base_color("DA");
    let col_t = base_color("DT");
    let col_g = base_color("DG");
    let col_c = base_color("DC");
    let col_u = base_color("U");

    assert_ne!(col_a, col_t);
    assert_ne!(col_g, col_c);
    assert_ne!(col_u, (0, 0, 0));
    assert!(col_a.1 > col_a.0); // Green dominant
    assert!(col_t.0 > col_t.1); // Red dominant
}

#[test]
fn test_nucleic_base_slab_generation() {
    let mut res = Residue::new(1, "DA", "A");
    let a_p = Atom::new(
        0,
        1,
        "P",
        element_by_symbol("P"),
        Vec3::new(0.0, 0.0, 0.0),
        10.0,
        "DA",
        1,
        "A",
        false,
    );
    let a_n9 = Atom::new(
        1,
        2,
        "N9",
        element_by_symbol("N"),
        Vec3::new(2.0, 0.0, 0.0),
        10.0,
        "DA",
        1,
        "A",
        false,
    );
    let a_c4 = Atom::new(
        2,
        3,
        "C4",
        element_by_symbol("C"),
        Vec3::new(3.0, 1.0, 0.0),
        10.0,
        "DA",
        1,
        "A",
        false,
    );
    res.atom_indices.push(0);
    res.atom_indices.push(1);
    res.atom_indices.push(2);

    let atoms = vec![a_p, a_n9, a_c4];
    let slab = build_base_slab(&res, &atoms, Vec3::new(0.0, 0.0, 0.0));
    assert!(slab.is_some());
}

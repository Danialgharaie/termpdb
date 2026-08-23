use clap::Parser;
use termpdb::cli::Cli;
use termpdb::math::Vec3;
use termpdb::model::align::superimpose_structures;
use termpdb::model::geometry::{calculate_bond_angle, calculate_dihedral_angle, classify_ramachandran};
use termpdb::model::{Atom, Chain, Element, Residue, Structure};

fn make_dummy_atom(idx: usize, pos: Vec3, name: &str) -> Atom {
    Atom::new(
        idx,
        (idx + 1) as i32,
        name,
        Element {
            symbol: "C",
            name: "Carbon",
            atomic_number: 6,
            covalent_radius: 0.76,
            vdw_radius: 1.7,
            cpk_color: (144, 144, 144),
        },
        pos,
        20.0,
        "ALA",
        1,
        "A",
        false,
    )
}

#[test]
fn test_cli_parsing_scientific_flags() {
    let args = vec![
        "termpdb",
        "struct1.pdb",
        "struct2.pdb",
        "--align",
        "--dssp",
        "--angle",
        "A:1:CA,A:2:CA,A:3:CA",
        "--dihedral",
        "A:1:N,A:1:CA,A:1:C,A:2:N",
    ];

    let cli = Cli::try_parse_from(args).expect("Failed to parse CLI with scientific flags");
    assert_eq!(cli.files.len(), 2);
    assert!(cli.align);
    assert!(cli.dssp);
    assert_eq!(cli.angle.as_deref(), Some("A:1:CA,A:2:CA,A:3:CA"));
    assert_eq!(cli.dihedral.as_deref(), Some("A:1:N,A:1:CA,A:1:C,A:2:N"));
}

#[test]
fn test_multi_structure_superposition_workflow() {
    let mut s1 = Structure::new("struct1");
    let mut c1 = Chain::new("A");
    for i in 0..5 {
        let mut r = Residue::new((i + 1) as i32, "ALA", "A");
        let atom = make_dummy_atom(i, Vec3::new(i as f32 * 3.8, 0.0, 0.0), "CA");
        r.atom_indices.push(i);
        s1.add_atom(atom);
        c1.residues.push(r);
    }
    s1.add_chain(c1);

    let mut s2 = Structure::new("struct2");
    let mut c2 = Chain::new("A");
    for i in 0..5 {
        let mut r = Residue::new((i + 1) as i32, "ALA", "A");
        let atom = make_dummy_atom(i, Vec3::new(i as f32 * 3.8 + 5.0, 2.0, -1.0), "CA");
        r.atom_indices.push(i);
        s2.add_atom(atom);
        c2.residues.push(r);
    }
    s2.add_chain(c2);

    let res = superimpose_structures(&mut s1, &s2).expect("Superposition should align structures");
    assert!(res.kabsch.rmsd < 1e-4);
    assert_eq!(res.aligned_pairs, 5);
}

#[test]
fn test_geometry_and_ramachandran_classification() {
    let p1 = Vec3::new(0.0, 1.0, 0.0);
    let p2 = Vec3::new(0.0, 0.0, 0.0);
    let p3 = Vec3::new(1.0, 0.0, 0.0);
    let p4 = Vec3::new(1.0, 1.0, 0.0);

    let angle = calculate_bond_angle(p1, p2, p3);
    assert!((angle - 90.0).abs() < 1e-3);

    let dihedral = calculate_dihedral_angle(p1, p2, p3, p4);
    assert!(dihedral.abs() < 1e-3);

    let rama = classify_ramachandran(-60.0, -45.0);
    assert_eq!(rama.name(), "α-Helix (Favored)");
}

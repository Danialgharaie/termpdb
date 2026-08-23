use termpdb::math::Vec3;
use termpdb::model::{Atom, Chain, Residue, Structure, element_by_symbol};
use termpdb::render::{Camera, Visibility};
use termpdb::select::{
    Selection, atom_distance, atom_label, distance_report, parse_atom_spec, pick_atom_at_screen,
    resolve_atom,
};

fn two_residue_structure() -> Structure {
    let mut structure = Structure::new("dipeptide");
    let mut chain = Chain::new("A");
    let mut r1 = Residue::new(1, "ALA", "A");
    let mut r2 = Residue::new(12, "GLY", "A");

    let c = element_by_symbol("C");
    let n = element_by_symbol("N");
    let h = element_by_symbol("H");
    let o = element_by_symbol("O");

    let i_n = structure.add_atom(Atom::new(
        0,
        1,
        "N",
        n,
        Vec3::new(-1.4, 0.0, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    ));
    let i_ca = structure.add_atom(Atom::new(
        0,
        2,
        "CA",
        c,
        Vec3::new(0.0, 0.0, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    ));
    let i_h = structure.add_atom(Atom::new(
        0,
        3,
        "H",
        h,
        Vec3::new(-1.9, 0.8, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    ));
    let i_ca2 = structure.add_atom(Atom::new(
        0,
        4,
        "CA",
        c,
        Vec3::new(3.8, 0.0, 0.0),
        10.0,
        "GLY",
        12,
        "A",
        false,
    ));

    r1.atom_indices.extend([i_n, i_ca, i_h]);
    r2.atom_indices.push(i_ca2);
    chain.residues.push(r1);

    let mut water = Residue::new(47, "HOH", "A");
    let i_o = structure.add_atom(Atom::new(
        0,
        5,
        "O",
        o,
        Vec3::new(20.0, 0.0, 0.0),
        10.0,
        "HOH",
        47,
        "A",
        true,
    ));
    water.atom_indices.push(i_o);
    chain.residues.push(r2);
    chain.residues.push(water);
    structure.add_chain(chain);
    structure
}

#[test]
fn test_parse_atom_spec_forms() {
    let a = parse_atom_spec("A:12:CA").unwrap();
    assert_eq!(a.chain_id.as_deref(), Some("A"));
    assert_eq!(a.res_seq, 12);
    assert_eq!(a.atom_name.as_deref(), Some("CA"));

    let b = parse_atom_spec("A/12/N").unwrap();
    assert_eq!(b.atom_name.as_deref(), Some("N"));

    let c = parse_atom_spec("12").unwrap();
    assert_eq!(c.chain_id, None);
    assert_eq!(c.res_seq, 12);
    assert_eq!(c.atom_name, None);

    let d = parse_atom_spec("A 12 CA").unwrap();
    assert_eq!(d.chain_id.as_deref(), Some("A"));
    assert_eq!(d.res_seq, 12);
    assert_eq!(d.atom_name.as_deref(), Some("CA"));

    let e = parse_atom_spec("12:N").unwrap();
    assert_eq!(e.chain_id, None);
    assert_eq!(e.res_seq, 12);
    assert_eq!(e.atom_name.as_deref(), Some("N"));

    assert!(parse_atom_spec("").is_err());
    assert!(parse_atom_spec("CA").is_err());
}

#[test]
fn test_resolve_prefers_ca_and_named_atoms() {
    let s = two_residue_structure();
    let ca = resolve_atom(&s, &parse_atom_spec("A:1").unwrap(), None).unwrap();
    assert_eq!(s.atoms()[ca].name, "CA");

    let n = resolve_atom(&s, &parse_atom_spec("A:1:N").unwrap(), None).unwrap();
    assert_eq!(s.atoms()[n].name, "N");

    let g = resolve_atom(&s, &parse_atom_spec("12").unwrap(), None).unwrap();
    assert_eq!(s.atoms()[g].res_seq, 12);
    assert_eq!(s.atoms()[g].name, "CA");

    assert!(resolve_atom(&s, &parse_atom_spec("A:99").unwrap(), None).is_err());
    assert!(resolve_atom(&s, &parse_atom_spec("A:1:ZZ").unwrap(), None).is_err());
}

#[test]
fn test_resolve_skips_hidden_hydrogens_unless_named() {
    let s = two_residue_structure();
    let vis = Visibility {
        show_waters: false,
        show_hydrogens: false,
    };
    let ca = resolve_atom(&s, &parse_atom_spec("1").unwrap(), Some(&vis)).unwrap();
    assert_eq!(s.atoms()[ca].name, "CA");

    let h = resolve_atom(&s, &parse_atom_spec("A:1:H").unwrap(), Some(&vis)).unwrap();
    assert_eq!(s.atoms()[h].name, "H");

    assert!(resolve_atom(&s, &parse_atom_spec("47").unwrap(), Some(&vis)).is_err());
    let o = resolve_atom(&s, &parse_atom_spec("A:47:O").unwrap(), Some(&vis)).unwrap();
    assert_eq!(s.atoms()[o].res_name, "HOH");
}

#[test]
fn test_distance_is_3_8_angstrom() {
    let s = two_residue_structure();
    let i = resolve_atom(&s, &parse_atom_spec("A:1:CA").unwrap(), None).unwrap();
    let j = resolve_atom(&s, &parse_atom_spec("A:12:CA").unwrap(), None).unwrap();
    let d = atom_distance(&s, i, j).unwrap();
    assert!((d - 3.8).abs() < 1e-4);

    let report = distance_report(&s, "A:1:CA,A:12:CA").unwrap();
    assert_eq!(report, "A:1:CA  A:12:CA  3.800");
    assert!(distance_report(&s, "A:1:CA").is_err());
}

#[test]
fn test_selection_fifo_and_status() {
    let s = two_residue_structure();
    let a = resolve_atom(&s, &parse_atom_spec("A:1:CA").unwrap(), None).unwrap();
    let b = resolve_atom(&s, &parse_atom_spec("A:12:CA").unwrap(), None).unwrap();
    let c = resolve_atom(&s, &parse_atom_spec("A:1:N").unwrap(), None).unwrap();

    let mut sel = Selection::default();
    sel.pick(a);
    assert_eq!(sel.status_line(&s).as_deref(), Some("A:1:CA"));
    sel.pick(b);
    let line = sel.status_line(&s).unwrap();
    assert!(line.contains("3.80 Å"), "{line}");
    sel.pick(c);
    assert_eq!(sel.atoms(), &[a, b, c]);
    sel.pick(c);
    assert_eq!(sel.atoms(), &[a, b]);
    sel.clear();
    assert!(sel.is_empty());
}

#[test]
fn test_atom_label() {
    let s = two_residue_structure();
    assert_eq!(atom_label(&s, 1), "A:1:CA");
}

#[test]
fn test_pick_atom_at_screen_hits_centered_atom() {
    let s = two_residue_structure();
    let mut camera = Camera::new();
    camera.fit_structure(s.center_of_mass(), s.bounding_sphere_radius());
    let width = 80;
    let height = 48;
    let ca = resolve_atom(&s, &parse_atom_spec("A:1:CA").unwrap(), None).unwrap();
    let (sx, sy, _) = camera
        .world_to_screen(s.atoms()[ca].pos, width, height)
        .expect("CA on screen");
    let picked = pick_atom_at_screen(&s, &camera, Visibility::ALL, width, height, sx, sy, 6.0);
    assert_eq!(picked, Some(ca));
    assert_eq!(
        pick_atom_at_screen(&s, &camera, Visibility::ALL, width, height, 0.0, 0.0, 1.0),
        None
    );
}

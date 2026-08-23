use termpdb::math::Vec3;
use termpdb::model::geometry::{
    RamachandranRegion, calculate_bond_angle, calculate_dihedral_angle, classify_ramachandran,
};
use termpdb::model::{Atom, Element, Structure};
use termpdb::select::Selection;

fn create_atom_at(idx: usize, pos: Vec3) -> Atom {
    Atom::new(
        idx,
        (idx + 1) as i32,
        "CA",
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
        (idx + 1) as i32,
        "A",
        false,
    )
}

#[test]
fn test_bond_angle_calculation() {
    let p1 = Vec3::new(1.0, 0.0, 0.0);
    let p2 = Vec3::new(0.0, 0.0, 0.0); // Vertex
    let p3 = Vec3::new(0.0, 1.0, 0.0);

    let angle_90 = calculate_bond_angle(p1, p2, p3);
    assert!((angle_90 - 90.0).abs() < 1e-4);

    let p4 = Vec3::new(-1.0, 0.0, 0.0);
    let angle_180 = calculate_bond_angle(p1, p2, p4);
    assert!((angle_180 - 180.0).abs() < 1e-4);
}

#[test]
fn test_dihedral_angle_calculation() {
    // Cis conformation: dihedral angle should be 0 deg
    let p1 = Vec3::new(1.0, 1.0, 0.0);
    let p2 = Vec3::new(0.0, 1.0, 0.0);
    let p3 = Vec3::new(0.0, 0.0, 0.0);
    let p4 = Vec3::new(1.0, 0.0, 0.0);

    let cis_angle = calculate_dihedral_angle(p1, p2, p3, p4);
    assert!(cis_angle.abs() < 1e-4, "Cis dihedral should be ~0 deg, got {cis_angle}");

    // Trans conformation: dihedral angle should be 180 deg
    let p4_trans = Vec3::new(-1.0, 0.0, 0.0);
    let trans_angle = calculate_dihedral_angle(p1, p2, p3, p4_trans);
    assert!((trans_angle.abs() - 180.0).abs() < 1e-4, "Trans dihedral should be 180 deg, got {trans_angle}");
}

#[test]
fn test_ramachandran_classification() {
    assert_eq!(
        classify_ramachandran(-60.0, -45.0),
        RamachandranRegion::AlphaHelix
    );
    assert_eq!(
        classify_ramachandran(-120.0, 135.0),
        RamachandranRegion::BetaSheet
    );
    assert_eq!(
        classify_ramachandran(60.0, 45.0),
        RamachandranRegion::LeftHandedAlpha
    );
    assert_eq!(
        classify_ramachandran(80.0, -120.0),
        RamachandranRegion::Outlier
    );
}

#[test]
fn test_selection_geometry_report() {
    let mut structure = Structure::new("test_geom");
    structure.add_atom(create_atom_at(0, Vec3::new(0.0, 0.0, 0.0)));
    structure.add_atom(create_atom_at(1, Vec3::new(1.5, 0.0, 0.0)));
    structure.add_atom(create_atom_at(2, Vec3::new(1.5, 1.5, 0.0)));
    structure.add_atom(create_atom_at(3, Vec3::new(0.0, 1.5, 0.0)));

    let mut sel = Selection::default();
    sel.pick(0);
    assert_eq!(sel.status(&structure), "Selected: A:1:CA");

    sel.pick(1);
    let stat2 = sel.status(&structure);
    assert!(stat2.contains("1.50 Å"), "2 atoms should show distance: {stat2}");

    sel.pick(2);
    let stat3 = sel.status(&structure);
    assert!(stat3.contains("Angle: 90.0°"), "3 atoms should show angle: {stat3}");

    sel.pick(3);
    let stat4 = sel.status(&structure);
    assert!(stat4.contains("Dihedral:"), "4 atoms should show dihedral: {stat4}");
}

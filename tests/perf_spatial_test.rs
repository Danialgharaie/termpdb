use termpdb::math::Vec3;
use termpdb::model::spatial::SpatialGrid;
use termpdb::model::{Atom, Element};

fn create_atom(index: usize, pos: Vec3) -> Atom {
    Atom::new(
        index,
        (index + 1) as i32,
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
        index as i32 + 1,
        "A",
        false,
    )
}

#[test]
fn test_spatial_grid_neighbor_queries() {
    let atoms = vec![
        create_atom(0, Vec3::new(0.0, 0.0, 0.0)),
        create_atom(1, Vec3::new(1.0, 0.0, 0.0)),
        create_atom(2, Vec3::new(0.0, 2.0, 0.0)),
        create_atom(3, Vec3::new(10.0, 10.0, 10.0)),
    ];

    let grid = SpatialGrid::new(&atoms, 3.0);
    assert_eq!(grid.cell_size(), 3.0);

    // Atom 0 at origin: neighbors within 2.5 A should include 0, 1, 2 but not 3
    let neighbors_0 = grid.neighbors_within(Vec3::new(0.0, 0.0, 0.0), 2.5);
    assert_eq!(neighbors_0.len(), 3);
    assert!(neighbors_0.contains(&0));
    assert!(neighbors_0.contains(&1));
    assert!(neighbors_0.contains(&2));
    assert!(!neighbors_0.contains(&3));

    // Atom 3 at (10, 10, 10): neighbor query should only find atom 3
    let neighbors_3 = grid.neighbors_within(Vec3::new(10.0, 10.0, 10.0), 2.5);
    assert_eq!(neighbors_3.len(), 1);
    assert_eq!(neighbors_3[0], 3);
}

#[test]
fn test_spatial_grid_burial_detection() {
    let mut atoms = Vec::new();
    let center = Vec3::new(0.0, 0.0, 0.0);
    atoms.push(create_atom(0, center)); // Core center atom

    // Surround center atom with 14 atoms in all directions at radius 3.0 A
    let offsets = [
        Vec3::new(3.0, 0.0, 0.0),
        Vec3::new(-3.0, 0.0, 0.0),
        Vec3::new(0.0, 3.0, 0.0),
        Vec3::new(0.0, -3.0, 0.0),
        Vec3::new(0.0, 0.0, 3.0),
        Vec3::new(0.0, 0.0, -3.0),
        Vec3::new(2.1, 2.1, 0.0),
        Vec3::new(-2.1, 2.1, 0.0),
        Vec3::new(2.1, -2.1, 0.0),
        Vec3::new(-2.1, -2.1, 0.0),
        Vec3::new(0.0, 2.1, 2.1),
        Vec3::new(0.0, -2.1, 2.1),
        Vec3::new(2.1, 0.0, 2.1),
        Vec3::new(-2.1, 0.0, 2.1),
    ];

    for (i, offset) in offsets.into_iter().enumerate() {
        atoms.push(create_atom(i + 1, offset));
    }

    let grid = SpatialGrid::new(&atoms, 4.0);
    let buried = grid.compute_buried_atoms(3.2, 10);

    // Atom 0 has 14 neighbors within 3.0 A (>= 10), so it is buried
    assert!(
        buried[0],
        "Center atom with 14 neighbors should be marked buried"
    );
    // Exterior atoms have far fewer surrounding neighbors, so they should NOT be buried
    assert!(!buried[1], "Exterior atom should not be marked buried");
}

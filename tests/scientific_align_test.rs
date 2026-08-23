use termpdb::math::Vec3;
use termpdb::model::align::{needleman_wunsch, pair_ca_coordinates, superimpose_structures};
use termpdb::model::{Atom, Chain, Element, Residue, Structure};

fn create_poly_ala(length: usize, offset: Vec3) -> Structure {
    let mut structure = Structure::new("poly_ala");
    let mut chain = Chain::new("A");

    for i in 0..length {
        let res_num = (i + 1) as i32;
        let mut res = Residue::new(res_num, "ALA", "A");
        let pos = Vec3::new(i as f32 * 3.8, 0.0, 0.0) + offset;
        let atom = Atom::new(
            i,
            (i + 1) as i32,
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
            res_num,
            "A",
            false,
        );
        res.atom_indices.push(i);
        structure.add_atom(atom);
        chain.residues.push(res);
    }
    structure.add_chain(chain);
    structure
}

#[test]
fn test_needleman_wunsch_alignment() {
    // Exact match
    let seq_a = "HEAGAWGHEE";
    let seq_b = "HEAGAWGHEE";
    let alignment = needleman_wunsch(seq_a, seq_b);
    assert_eq!(alignment.len(), 10);
    for (i, &(opt_a, opt_b)) in alignment.iter().enumerate() {
        assert_eq!(opt_a, Some(i));
        assert_eq!(opt_b, Some(i));
    }

    // Sequence with deletion: seq_b missing 'W'
    let seq_c = "HEAGAGHEE";
    let alignment_gap = needleman_wunsch(seq_a, seq_c);
    assert!(!alignment_gap.is_empty());
}

#[test]
fn test_pair_ca_coordinates() {
    let s1 = create_poly_ala(5, Vec3::ZERO);
    let s2 = create_poly_ala(5, Vec3::new(10.0, 0.0, 0.0));

    let pairs = pair_ca_coordinates(&s1, &s2);
    assert_eq!(pairs.len(), 5);

    for (i, &(p1, p2)) in pairs.iter().enumerate() {
        assert_eq!(p1, Vec3::new(i as f32 * 3.8, 0.0, 0.0));
        assert_eq!(p2, Vec3::new(i as f32 * 3.8 + 10.0, 0.0, 0.0));
    }
}

#[test]
fn test_superimpose_structures() {
    let mut s1 = create_poly_ala(6, Vec3::new(0.0, 0.0, 0.0));
    // s2 is translated by (2, 3, 4)
    let s2 = create_poly_ala(6, Vec3::new(2.0, 3.0, 4.0));

    let res = superimpose_structures(&mut s1, &s2).expect("Superposition should succeed");
    assert!(res.kabsch.rmsd < 1e-4, "Superposition of rigid translation should have 0 RMSD");
    assert_eq!(res.per_residue_rmsd.len(), 6);
    for &rmsd in &res.per_residue_rmsd {
        assert!(rmsd < 1e-4);
    }

    // Coordinates of s1 should now match s2
    for (a1, a2) in s1.atoms().iter().zip(s2.atoms().iter()) {
        assert!(a1.pos.distance(&a2.pos) < 1e-3);
    }
}

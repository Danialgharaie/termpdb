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
    assert!(
        helices > 0,
        "Should assign SecondaryStructure::Helix to helical residues"
    );
}

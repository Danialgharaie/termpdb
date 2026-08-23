use termpdb::math::Vec3;
use termpdb::model::{Atom, Chain, Residue, SecondaryStructure, Structure, element_by_symbol};
use termpdb::parser::parse_pdb;
use termpdb::render::representations::{
    LOD_BACKBONE_ATOMS, LOD_CALPHA_ATOMS, LodLevel, atom_passes_lod,
};
use termpdb::render::{
    Camera, ColorScheme, Framebuffer, Lighting, LodMode, RenderMode, Visibility, render_structure,
};

fn create_synthetic_protein() -> Structure {
    let mut structure = Structure::new("Synthetic Multi-Domain Protein");

    let elem_c = element_by_symbol("C");
    let elem_n = element_by_symbol("N");
    let elem_o = element_by_symbol("O");
    let elem_zn = element_by_symbol("Zn");

    let mut chain_a = Chain::new("A");

    // Build 15 residues in Chain A:
    // Residues 1..=5: Alpha-Helix (coiled coordinates)
    // Residues 6..=10: Beta-Sheet (extended zigzag coordinates)
    // Residues 11..=15: Coil / Loop (extended curved coordinates)
    for res_idx in 1..=15 {
        let (res_name, ss) = if res_idx <= 5 {
            ("ALA", SecondaryStructure::Helix)
        } else if res_idx <= 10 {
            ("VAL", SecondaryStructure::Sheet)
        } else {
            ("GLY", SecondaryStructure::Coil)
        };

        let mut residue = Residue::new(res_idx, res_name, "A");
        residue.secondary_structure = ss;

        // Position coordinates based on secondary structure
        let ca_pos = if res_idx <= 5 {
            // Helical spiral
            let angle = (res_idx as f32) * 1.7; // ~100 deg per residue
            Vec3::new(angle.cos() * 2.3, angle.sin() * 2.3, (res_idx as f32) * 1.5)
        } else if res_idx <= 10 {
            // Beta strand zigzag
            let sign = if res_idx % 2 == 0 { 1.0 } else { -1.0 };
            Vec3::new(sign * 1.2, 5.0 + (res_idx as f32) * 2.5, 7.5)
        } else {
            // Coil curve
            let t = (res_idx - 10) as f32;
            Vec3::new(
                t * 2.0,
                18.0 + (t * 0.8).sin() * 3.0,
                7.5 + (t * 0.8).cos() * 3.0,
            )
        };

        let n_pos = ca_pos + Vec3::new(-0.8, -0.6, 0.4);
        let c_pos = ca_pos + Vec3::new(0.8, 0.6, -0.4);
        let o_pos = c_pos + Vec3::new(0.2, 1.2, 0.0);

        let atom_n = Atom::new(
            0,
            (res_idx - 1) * 4 + 1,
            "N",
            elem_n,
            n_pos,
            15.0,
            res_name,
            res_idx,
            "A",
            false,
        );
        let atom_ca = Atom::new(
            0,
            (res_idx - 1) * 4 + 2,
            "CA",
            elem_c,
            ca_pos,
            12.0,
            res_name,
            res_idx,
            "A",
            false,
        );
        let atom_c = Atom::new(
            0,
            (res_idx - 1) * 4 + 3,
            "C",
            elem_c,
            c_pos,
            12.0,
            res_name,
            res_idx,
            "A",
            false,
        );
        let atom_o = Atom::new(
            0,
            (res_idx - 1) * 4 + 4,
            "O",
            elem_o,
            o_pos,
            18.0,
            res_name,
            res_idx,
            "A",
            false,
        );

        let idx_n = structure.add_atom(atom_n);
        let idx_ca = structure.add_atom(atom_ca);
        let idx_c = structure.add_atom(atom_c);
        let idx_o = structure.add_atom(atom_o);

        residue.atom_indices.push(idx_n);
        residue.atom_indices.push(idx_ca);
        residue.atom_indices.push(idx_c);
        residue.atom_indices.push(idx_o);

        chain_a.residues.push(residue);
    }

    structure.add_chain(chain_a);

    // Add a HETATM ligand (Zinc ion)
    let atom_zn = Atom::new(
        0,
        100,
        "ZN",
        elem_zn,
        Vec3::new(5.0, 5.0, 5.0),
        25.0,
        "ZN",
        99,
        "A",
        true,
    );
    structure.add_atom(atom_zn);

    // Automatically build bonds
    structure.build_bonds();

    structure
}

fn count_drawn_pixels(buffer: &Framebuffer, bg: (u8, u8, u8)) -> usize {
    buffer.pixels.iter().filter(|&&p| p != bg).count()
}

#[test]
fn test_render_mode_enum_methods() {
    let modes = RenderMode::all();
    assert_eq!(modes.len(), 5);
    assert!(modes.contains(&RenderMode::Trace));
    assert!(modes.contains(&RenderMode::BallAndStick));
    assert!(modes.contains(&RenderMode::Ribbon));
    assert!(modes.contains(&RenderMode::Vdw));
    assert!(modes.contains(&RenderMode::Wireframe));

    for mode in modes {
        assert!(!mode.name().is_empty());
        assert_eq!(mode.next().prev(), *mode);
        assert_eq!(mode.prev().next(), *mode);
    }

    // Verify cyclic ordering
    let mode = RenderMode::Trace;
    let next = mode.next();
    assert_ne!(mode, next);
}

#[test]
fn test_render_empty_structure() {
    let structure = Structure::new("Empty");
    let mut camera = Camera::new();
    camera.fit_structure(Vec3::ZERO, 10.0);
    let mut buffer = Framebuffer::new(40, 40);
    let bg = (0, 0, 0);
    buffer.clear(bg);
    let lighting = Lighting::default();

    for &mode in RenderMode::all() {
        render_structure(
            &structure,
            mode,
            ColorScheme::Cpk,
            &camera,
            &mut buffer,
            &lighting,
            Visibility::ALL,
            LodMode::Full,
        );
        assert_eq!(count_drawn_pixels(&buffer, bg), 0);
    }
}

#[test]
fn test_render_trace_mode() {
    let structure = create_synthetic_protein();
    let mut camera = Camera::new();
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    camera.fit_structure(com, radius);

    let mut buffer = Framebuffer::new(80, 48);
    let bg = (10, 10, 10);
    buffer.clear(bg);
    let lighting = Lighting::default();

    render_structure(
        &structure,
        RenderMode::Trace,
        ColorScheme::SecondaryStructure,
        &camera,
        &mut buffer,
        &lighting,
        Visibility::ALL,
        LodMode::Full,
    );

    let drawn = count_drawn_pixels(&buffer, bg);
    assert!(
        drawn > 50,
        "Trace mode should render backbone lines and ligand atoms, drawn={}",
        drawn
    );
}

#[test]
fn test_render_ball_and_stick_mode() {
    let structure = create_synthetic_protein();
    let mut camera = Camera::new();
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    camera.fit_structure(com, radius);

    let mut buffer = Framebuffer::new(80, 48);
    let bg = (10, 10, 10);
    buffer.clear(bg);
    let lighting = Lighting::default();

    render_structure(
        &structure,
        RenderMode::BallAndStick,
        ColorScheme::Cpk,
        &camera,
        &mut buffer,
        &lighting,
        Visibility::ALL,
        LodMode::Full,
    );

    let drawn = count_drawn_pixels(&buffer, bg);
    assert!(
        drawn > 100,
        "Ball & Stick mode should render atoms and bonds, drawn={}",
        drawn
    );
}

#[test]
fn test_render_ribbon_mode() {
    let structure = create_synthetic_protein();
    let mut camera = Camera::new();
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    camera.fit_structure(com, radius);

    let mut buffer = Framebuffer::new(80, 48);
    let bg = (10, 10, 10);
    buffer.clear(bg);
    let lighting = Lighting::default();

    render_structure(
        &structure,
        RenderMode::Ribbon,
        ColorScheme::SecondaryStructure,
        &camera,
        &mut buffer,
        &lighting,
        Visibility::ALL,
        LodMode::Full,
    );

    let drawn = count_drawn_pixels(&buffer, bg);
    assert!(
        drawn > 100,
        "Ribbon mode should render cartoon ribbons and ligand atoms, drawn={}",
        drawn
    );
}

#[test]
fn test_render_vdw_mode() {
    let structure = create_synthetic_protein();
    let mut camera = Camera::new();
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    camera.fit_structure(com, radius);

    let mut buffer_vdw = Framebuffer::new(80, 48);
    let mut buffer_bas = Framebuffer::new(80, 48);
    let bg = (10, 10, 10);
    buffer_vdw.clear(bg);
    buffer_bas.clear(bg);
    let lighting = Lighting::default();

    render_structure(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        &camera,
        &mut buffer_vdw,
        &lighting,
        Visibility::ALL,
        LodMode::Full,
    );

    render_structure(
        &structure,
        RenderMode::BallAndStick,
        ColorScheme::Cpk,
        &camera,
        &mut buffer_bas,
        &lighting,
        Visibility::ALL,
        LodMode::Full,
    );

    let drawn_vdw = count_drawn_pixels(&buffer_vdw, bg);
    let drawn_bas = count_drawn_pixels(&buffer_bas, bg);

    assert!(
        drawn_vdw > 100,
        "VDW mode should render full space-filling spheres, drawn={}",
        drawn_vdw
    );
    // VDW spheres have much larger radius than ball-and-stick spheres
    assert!(
        drawn_vdw > drawn_bas,
        "VDW space-filling ({}) should cover more pixels than Ball&Stick ({})",
        drawn_vdw,
        drawn_bas
    );
}

#[test]
fn test_render_all_modes_and_all_color_schemes() {
    let structure = create_synthetic_protein();
    let mut camera = Camera::new();
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    camera.fit_structure(com, radius);

    let lighting = Lighting::default();
    let bg = (0, 0, 0);

    for &mode in RenderMode::all() {
        for &scheme in ColorScheme::all() {
            let mut buffer = Framebuffer::new(60, 36);
            buffer.clear(bg);

            render_structure(
                &structure,
                mode,
                scheme,
                &camera,
                &mut buffer,
                &lighting,
                Visibility::ALL,
                LodMode::Full,
            );

            let drawn = count_drawn_pixels(&buffer, bg);
            assert!(
                drawn > 0,
                "Rendering with mode {:?} and scheme {:?} produced 0 drawn pixels",
                mode,
                scheme
            );
        }
    }
}

#[test]
fn test_render_real_crambin_pdb() {
    let pdb_text = r#"HEADER    PLANT SEED PROTEIN                      30-APR-81   1CRN
TITLE     WATER STRUCTURE OF A HYDROPHOBIC PROTEIN AT ATOMIC RESOLUTION.
HELIX    1  H1 THR A    7  GLY A   17  1                                  11
SHEET    1  S1 2 CYS A   1  CYS A    4  0
SHEET    2  S1 2 CYS A  32  ILE A   35 -1
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O
ATOM      5  N   THR A   2      15.115  11.555   5.265  1.00  7.81           N
ATOM      6  CA  THR A   2      13.856  11.469   6.066  1.00  8.28           C
ATOM      7  C   THR A   2      12.639  12.016   5.297  1.00  7.24           C
ATOM      8  O   THR A   2      12.753  13.064   4.646  1.00  7.56           O
ATOM      9  N   GLY A   7      10.000  10.000  10.000  1.00 15.00           N
ATOM     10  CA  GLY A   7      10.800  10.500  11.000  1.00 15.00           C
ATOM     11  C   GLY A   7      11.500  11.200  11.800  1.00 15.00           C
ATOM     12  O   GLY A   7      11.800  12.200  11.500  1.00 15.00           O
HETATM   13  O   HOH A  47      15.021  -1.423   2.492  1.00 25.00           O
END
"#;

    let mut structure = parse_pdb(pdb_text).expect("Failed to parse PDB");
    structure.build_bonds();

    let mut camera = Camera::new();
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    camera.fit_structure(com, radius);
    let lighting = Lighting::default();

    for &mode in RenderMode::all() {
        let mut buffer = Framebuffer::new(60, 36);
        buffer.clear((0, 0, 0));

        render_structure(
            &structure,
            mode,
            ColorScheme::Cpk,
            &camera,
            &mut buffer,
            &lighting,
            Visibility::ALL,
            LodMode::Full,
        );

        let drawn = count_drawn_pixels(&buffer, (0, 0, 0));
        assert!(
            drawn > 0,
            "1CRN rendering in mode {:?} should produce drawn pixels",
            mode
        );
    }
}

fn visibility_fixture() -> Structure {
    let mut structure = Structure::new("vis");
    let mut chain = Chain::new("A");
    let mut protein = Residue::new(1, "ALA", "A");
    let mut water = Residue::new(2, "HOH", "A");

    let ca = structure.add_atom(Atom::new(
        0,
        1,
        "CA",
        element_by_symbol("C"),
        Vec3::new(0.0, 0.0, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    ));
    let h = structure.add_atom(Atom::new(
        0,
        2,
        "H",
        element_by_symbol("H"),
        Vec3::new(8.0, 0.0, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    ));
    let o = structure.add_atom(Atom::new(
        0,
        3,
        "O",
        element_by_symbol("O"),
        Vec3::new(12.0, 0.0, 0.0),
        10.0,
        "HOH",
        2,
        "A",
        true,
    ));

    protein.atom_indices.extend([ca, h]);
    water.atom_indices.push(o);
    chain.residues.push(protein);
    chain.residues.push(water);
    structure.add_chain(chain);
    structure
}

#[test]
fn test_hiding_waters_draws_fewer_vdw_pixels_than_showing_all() {
    let structure = visibility_fixture();
    let mut camera = Camera::new();
    camera.fit_structure(
        structure.center_of_mass(),
        structure.bounding_sphere_radius(),
    );
    let lighting = Lighting::default();

    let mut all_buf = Framebuffer::new(80, 48);
    all_buf.clear((0, 0, 0));
    render_structure(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        &camera,
        &mut all_buf,
        &lighting,
        Visibility::ALL,
        LodMode::Full,
    );

    let mut hidden_buf = Framebuffer::new(80, 48);
    hidden_buf.clear((0, 0, 0));
    render_structure(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        &camera,
        &mut hidden_buf,
        &lighting,
        Visibility::default(),
        LodMode::Full,
    );

    let all = count_drawn_pixels(&all_buf, (0, 0, 0));
    let hidden = count_drawn_pixels(&hidden_buf, (0, 0, 0));
    assert!(all > 0);
    assert!(
        hidden < all,
        "default visibility should hide waters: all={all} hidden={hidden}"
    );
}

#[test]
fn test_hiding_hydrogens_draws_fewer_vdw_pixels() {
    let structure = visibility_fixture();
    let mut camera = Camera::new();
    camera.fit_structure(
        structure.center_of_mass(),
        structure.bounding_sphere_radius(),
    );
    let lighting = Lighting::default();

    let mut with_h = Framebuffer::new(80, 48);
    with_h.clear((0, 0, 0));
    render_structure(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        &camera,
        &mut with_h,
        &lighting,
        Visibility {
            show_waters: false,
            show_hydrogens: true,
        },
        LodMode::Full,
    );

    let mut no_h = Framebuffer::new(80, 48);
    no_h.clear((0, 0, 0));
    render_structure(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        &camera,
        &mut no_h,
        &lighting,
        Visibility {
            show_waters: false,
            show_hydrogens: false,
        },
        LodMode::Full,
    );

    let a = count_drawn_pixels(&with_h, (0, 0, 0));
    let b = count_drawn_pixels(&no_h, (0, 0, 0));
    assert!(
        b < a,
        "hiding hydrogens should draw fewer pixels: with_h={a} no_h={b}"
    );
}

#[test]
fn test_lod_auto_thresholds() {
    assert_eq!(LodMode::Auto.resolve(0), LodLevel::Full);
    assert_eq!(
        LodMode::Auto.resolve(LOD_BACKBONE_ATOMS - 1),
        LodLevel::Full
    );
    assert_eq!(
        LodMode::Auto.resolve(LOD_BACKBONE_ATOMS),
        LodLevel::Backbone
    );
    assert_eq!(LodMode::Auto.resolve(LOD_CALPHA_ATOMS), LodLevel::CAlpha);
    assert_eq!(LodMode::Full.resolve(1_000_000), LodLevel::Full);
    assert_eq!(LodMode::Auto.next(), LodMode::Full);
    assert_eq!(LodMode::CAlpha.next(), LodMode::Auto);
}

#[test]
fn test_atom_passes_lod_keeps_backbone_and_ligands() {
    let c = element_by_symbol("C");
    let zn = element_by_symbol("Zn");
    let ca = Atom::new(0, 1, "CA", c, Vec3::ZERO, 0.0, "ALA", 1, "A", false);
    let cb = Atom::new(1, 2, "CB", c, Vec3::ZERO, 0.0, "ALA", 1, "A", false);
    let ligand = Atom::new(2, 3, "ZN", zn, Vec3::ZERO, 0.0, "ZN", 99, "A", true);
    let mut ala = Residue::new(1, "ALA", "A");
    ala.atom_indices = vec![0, 1];
    let mut zn_res = Residue::new(99, "ZN", "A");
    zn_res.atom_indices = vec![2];

    assert!(atom_passes_lod(&ca, Some(&ala), LodLevel::Full));
    assert!(atom_passes_lod(&cb, Some(&ala), LodLevel::Full));
    assert!(atom_passes_lod(&ca, Some(&ala), LodLevel::Backbone));
    assert!(!atom_passes_lod(&cb, Some(&ala), LodLevel::Backbone));
    assert!(atom_passes_lod(&ligand, Some(&zn_res), LodLevel::Backbone));
    assert!(atom_passes_lod(&ca, Some(&ala), LodLevel::CAlpha));
    assert!(!atom_passes_lod(&cb, Some(&ala), LodLevel::CAlpha));
    assert!(atom_passes_lod(&ligand, Some(&zn_res), LodLevel::CAlpha));
}

#[test]
fn test_lod_backbone_draws_fewer_vdw_pixels_than_full() {
    let mut structure = Structure::new("lod");
    let mut chain = Chain::new("A");
    let mut res = Residue::new(1, "ALA", "A");
    let c = element_by_symbol("C");
    let ca = structure.add_atom(Atom::new(
        0,
        1,
        "CA",
        c,
        Vec3::new(0.0, 0.0, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    ));
    let cb = structure.add_atom(Atom::new(
        0,
        2,
        "CB",
        c,
        Vec3::new(8.0, 0.0, 0.0),
        10.0,
        "ALA",
        1,
        "A",
        false,
    ));
    res.atom_indices.extend([ca, cb]);
    chain.residues.push(res);
    structure.add_chain(chain);

    let mut camera = Camera::new();
    camera.fit_structure(
        structure.center_of_mass(),
        structure.bounding_sphere_radius(),
    );
    let lighting = Lighting::default();

    let mut full = Framebuffer::new(80, 48);
    full.clear((0, 0, 0));
    render_structure(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        &camera,
        &mut full,
        &lighting,
        Visibility::ALL,
        LodMode::Full,
    );
    let mut bb = Framebuffer::new(80, 48);
    bb.clear((0, 0, 0));
    render_structure(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        &camera,
        &mut bb,
        &lighting,
        Visibility::ALL,
        LodMode::Backbone,
    );
    let a = count_drawn_pixels(&full, (0, 0, 0));
    let b = count_drawn_pixels(&bb, (0, 0, 0));
    assert!(b < a, "backbone LOD should skip CB: full={a} bb={b}");
}

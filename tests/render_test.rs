use termpdb::math::Vec3;
use termpdb::model::{Atom, Chain, Element, Residue, SecondaryStructure, Structure};
use termpdb::render::{
    Camera, ColorScheme, Framebuffer, Lighting, PixelColor, color_for_atom, draw_cylinder,
    draw_line_3d, draw_sphere, draw_triangle_3d,
};

#[test]
fn test_framebuffer_creation_and_clear() {
    let mut fb = Framebuffer::new(20, 20);
    assert_eq!(fb.width, 20);
    assert_eq!(fb.height, 20);

    let bg: PixelColor = (10, 20, 30);
    fb.clear(bg);

    for y in 0..20 {
        for x in 0..20 {
            assert_eq!(fb.get_pixel(x, y), Some(bg));
            assert_eq!(fb.get_depth(x, y), Some(f32::INFINITY));
        }
    }
}

#[test]
fn test_framebuffer_z_buffer_occlusion() {
    let mut fb = Framebuffer::new(10, 10);
    fb.clear((0, 0, 0));

    let red: PixelColor = (255, 0, 0);
    let blue: PixelColor = (0, 0, 255);
    let green: PixelColor = (0, 255, 0);

    // Far pixel at z = 10.0
    assert!(fb.set_pixel(5, 5, 10.0, blue));
    assert_eq!(fb.get_pixel(5, 5), Some(blue));
    assert_eq!(fb.get_depth(5, 5), Some(10.0));

    // Near pixel at z = 5.0 should overwrite far pixel
    assert!(fb.set_pixel(5, 5, 5.0, red));
    assert_eq!(fb.get_pixel(5, 5), Some(red));
    assert_eq!(fb.get_depth(5, 5), Some(5.0));

    // Even farther pixel at z = 15.0 should be rejected and NOT overwrite
    assert!(!fb.set_pixel(5, 5, 15.0, green));
    assert_eq!(fb.get_pixel(5, 5), Some(red));
    assert_eq!(fb.get_depth(5, 5), Some(5.0));
}

#[test]
fn test_framebuffer_out_of_bounds() {
    let mut fb = Framebuffer::new(10, 10);
    fb.clear((0, 0, 0));

    assert!(!fb.set_pixel(-1, 5, 1.0, (255, 255, 255)));
    assert!(!fb.set_pixel(10, 5, 1.0, (255, 255, 255)));
    assert!(!fb.set_pixel(5, -1, 1.0, (255, 255, 255)));
    assert!(!fb.set_pixel(5, 10, 1.0, (255, 255, 255)));

    assert_eq!(fb.get_pixel(15, 15), None);
    assert_eq!(fb.get_depth(15, 15), None);
}

#[test]
fn test_framebuffer_half_block_cells() {
    let mut fb = Framebuffer::new(4, 4);
    fb.clear((0, 0, 0));

    let top_color: PixelColor = (255, 0, 0);
    let bottom_color: PixelColor = (0, 255, 0);

    // Terminal row 0 corresponds to pixel y=0 (top) and y=1 (bottom)
    fb.set_pixel(1, 0, 1.0, top_color);
    fb.set_pixel(1, 1, 1.0, bottom_color);

    let cell = fb.cell_at(1, 0);
    assert_eq!(cell.0, '▀');
    assert_eq!(cell.1, top_color);
    assert_eq!(cell.2, Some(bottom_color));

    let blocks = fb.get_half_blocks();
    assert_eq!(blocks.len(), 2); // 4 pixels height / 2 = 2 terminal rows
    assert_eq!(blocks[0].len(), 4); // 4 terminal cols
    assert_eq!(blocks[0][1], ('▀', top_color, Some(bottom_color)));
}

#[test]
fn test_framebuffer_resize() {
    let mut fb = Framebuffer::new(10, 10);
    fb.clear((10, 10, 10));
    assert_eq!(fb.width, 10);
    assert_eq!(fb.height, 10);

    fb.resize(20, 30);
    assert_eq!(fb.width, 20);
    assert_eq!(fb.height, 30);
    assert_eq!(fb.get_pixel(19, 29), Some((0, 0, 0)));
}

#[test]
fn test_camera_initialization_and_matrices() {
    let camera = Camera::new();
    assert!(camera.distance > 0.0);
    assert!(camera.near > 0.0);
    assert!(camera.far > camera.near);

    let view = camera.view_matrix();
    let proj = camera.proj_matrix();

    // Invertible matrices
    assert!(view.inverse().is_some());
    assert!(proj.inverse().is_some());
}

#[test]
fn test_camera_fit_structure() {
    let mut camera = Camera::new();
    let center = Vec3::new(10.0, 20.0, 30.0);
    let radius = 15.0;

    camera.fit_structure(center, radius);

    assert_eq!(camera.target, center);
    assert!(camera.distance > radius);
    assert!(camera.near < camera.distance);
    assert!(camera.far > camera.distance);
}

#[test]
fn test_camera_orbit_pan_zoom() {
    let mut camera = Camera::new();
    let initial_dist = camera.distance;
    let initial_target = camera.target;

    camera.orbit(10.0, 5.0);
    camera.pan(2.0, -3.0);
    assert_ne!(camera.target, initial_target);

    camera.zoom(0.5);
    assert_ne!(camera.distance, initial_dist);

    camera.reset();
    assert_eq!(camera.target, Vec3::ZERO);
}

#[test]
fn test_camera_world_to_screen() {
    let mut camera = Camera::new();
    camera.aspect = 1.0;
    camera.fit_structure(Vec3::ZERO, 10.0);

    let width = 80;
    let height = 48;

    // Origin (target) should project to approximately screen center
    let screen_pt = camera.world_to_screen(Vec3::ZERO, width, height);
    assert!(screen_pt.is_some());
    let (sx, sy, sz) = screen_pt.unwrap();
    assert!((sx - 40.0).abs() < 1.0);
    assert!((sy - 24.0).abs() < 1.0);
    assert!(sz > 0.0);

    // Point far behind camera should return None
    let behind = camera.eye_position() + (camera.eye_position() - camera.target);
    let behind_proj = camera.world_to_screen(behind, width, height);
    assert!(behind_proj.is_none());
}

#[test]
fn test_lighting_shade_and_fog() {
    let lighting = Lighting::default();
    let base_color: PixelColor = (200, 100, 50);

    // Normal directly facing the light vs normal facing away
    let n_front = lighting.light_dir;
    let n_back = -lighting.light_dir;

    let lit_front = lighting.shade(n_front, 10.0, base_color, 5.0, 50.0);
    let lit_back = lighting.shade(n_back, 10.0, base_color, 5.0, 50.0);

    // Front should be brighter than back
    assert!(lit_front.0 > lit_back.0);
    assert!(lit_front.1 > lit_back.1);
    assert!(lit_front.2 > lit_back.2);

    // Back should still have non-zero ambient light
    assert!(lit_back.0 > 0);

    // Depth fog: farther depth should be dimmer
    let lit_far = lighting.shade(n_front, 45.0, base_color, 5.0, 50.0);
    assert!(lit_front.0 >= lit_far.0);
}

#[test]
fn test_color_schemes_cpk() {
    let carbon = Atom::new(
        0,
        1,
        "CA",
        Element {
            symbol: "C",
            name: "Carbon",
            atomic_number: 6,
            covalent_radius: 0.76,
            vdw_radius: 1.7,
            cpk_color: (144, 144, 144),
        },
        Vec3::ZERO,
        20.0,
        "ALA",
        1,
        "A",
        false,
    );
    let oxygen = Atom::new(
        1,
        2,
        "O",
        Element {
            symbol: "O",
            name: "Oxygen",
            atomic_number: 8,
            covalent_radius: 0.66,
            vdw_radius: 1.52,
            cpk_color: (255, 13, 13),
        },
        Vec3::ZERO,
        20.0,
        "ALA",
        1,
        "A",
        false,
    );
    let nitrogen = Atom::new(
        2,
        3,
        "N",
        Element {
            symbol: "N",
            name: "Nitrogen",
            atomic_number: 7,
            covalent_radius: 0.71,
            vdw_radius: 1.55,
            cpk_color: (48, 80, 248),
        },
        Vec3::ZERO,
        20.0,
        "ALA",
        1,
        "A",
        false,
    );

    let struct_empty = Structure::new("test");

    assert_eq!(
        color_for_atom(&carbon, None, &struct_empty, ColorScheme::Cpk),
        (144, 144, 144)
    );
    assert_eq!(
        color_for_atom(&oxygen, None, &struct_empty, ColorScheme::Cpk),
        (255, 13, 13)
    );
    assert_eq!(
        color_for_atom(&nitrogen, None, &struct_empty, ColorScheme::Cpk),
        (48, 80, 248)
    );
}

#[test]
fn test_color_schemes_secondary_structure() {
    let atom = Atom::new(
        0,
        1,
        "CA",
        Element::unknown(),
        Vec3::ZERO,
        10.0,
        "ALA",
        1,
        "A",
        false,
    );
    let struct_empty = Structure::new("test");

    let mut res_helix = Residue::new(1, "ALA", "A");
    res_helix.secondary_structure = SecondaryStructure::Helix;

    let mut res_sheet = Residue::new(2, "VAL", "A");
    res_sheet.secondary_structure = SecondaryStructure::Sheet;

    let mut res_coil = Residue::new(3, "GLY", "A");
    res_coil.secondary_structure = SecondaryStructure::Coil;

    let c_helix = color_for_atom(
        &atom,
        Some(&res_helix),
        &struct_empty,
        ColorScheme::SecondaryStructure,
    );
    let c_sheet = color_for_atom(
        &atom,
        Some(&res_sheet),
        &struct_empty,
        ColorScheme::SecondaryStructure,
    );
    let c_coil = color_for_atom(
        &atom,
        Some(&res_coil),
        &struct_empty,
        ColorScheme::SecondaryStructure,
    );

    assert_ne!(c_helix, c_sheet);
    assert_ne!(c_helix, c_coil);
    assert_ne!(c_sheet, c_coil);
}

#[test]
fn test_color_schemes_chain_and_cycle() {
    let atom_a = Atom::new(
        0,
        1,
        "CA",
        Element::unknown(),
        Vec3::ZERO,
        10.0,
        "ALA",
        1,
        "A",
        false,
    );
    let atom_b = Atom::new(
        1,
        2,
        "CA",
        Element::unknown(),
        Vec3::ZERO,
        10.0,
        "ALA",
        1,
        "B",
        false,
    );
    let mut structure = Structure::new("test");
    structure.add_chain(Chain::new("A"));
    structure.add_chain(Chain::new("B"));

    let c_a = color_for_atom(&atom_a, None, &structure, ColorScheme::Chain);
    let c_b = color_for_atom(&atom_b, None, &structure, ColorScheme::Chain);
    assert_ne!(c_a, c_b);

    // Test cycle
    let scheme = ColorScheme::Cpk;
    let next_scheme = scheme.next();
    assert_eq!(next_scheme, ColorScheme::Rainbow);
    assert_eq!(scheme.prev().next(), scheme);
    assert!(!scheme.name().is_empty());
}

#[test]
fn test_color_schemes_hydrophobicity_and_bfactor() {
    let mut structure = Structure::new("test");
    let atom_cold = Atom::new(
        0,
        1,
        "CA",
        Element::unknown(),
        Vec3::ZERO,
        10.0,
        "ILE",
        1,
        "A",
        false,
    );
    let atom_hot = Atom::new(
        1,
        2,
        "CA",
        Element::unknown(),
        Vec3::ZERO,
        90.0,
        "ARG",
        2,
        "A",
        false,
    );
    structure.add_atom(atom_cold.clone());
    structure.add_atom(atom_hot.clone());

    let res_hydrophobic = Residue::new(1, "ILE", "A");
    let res_hydrophilic = Residue::new(2, "ARG", "A");

    let c_hp = color_for_atom(
        &atom_cold,
        Some(&res_hydrophobic),
        &structure,
        ColorScheme::Hydrophobicity,
    );
    let c_hl = color_for_atom(
        &atom_hot,
        Some(&res_hydrophilic),
        &structure,
        ColorScheme::Hydrophobicity,
    );
    assert_ne!(c_hp, c_hl);

    let c_cold = color_for_atom(&atom_cold, None, &structure, ColorScheme::BFactor);
    let c_hot = color_for_atom(&atom_hot, None, &structure, ColorScheme::BFactor);
    assert_ne!(c_cold, c_hot);
}

#[test]
fn test_rasterizer_draw_sphere() {
    let mut fb = Framebuffer::new(30, 30);
    fb.clear((0, 0, 0));
    let lighting = Lighting::default();

    draw_sphere(&mut fb, (15.0, 15.0, 10.0), 8.0, (200, 50, 50), &lighting);

    // Center pixel should be drawn and have depth approx center_z - radius
    let center_color = fb.get_pixel(15, 15).unwrap();
    assert_ne!(center_color, (0, 0, 0));
    let center_depth = fb.get_depth(15, 15).unwrap();
    assert!((center_depth - 2.0).abs() < 0.5); // 10.0 - 8.0 = 2.0

    // Pixel far outside radius should be background and infinite depth
    assert_eq!(fb.get_pixel(0, 0).unwrap(), (0, 0, 0));
    assert_eq!(fb.get_depth(0, 0).unwrap(), f32::INFINITY);
}

#[test]
fn test_rasterizer_sphere_occlusion() {
    let mut fb = Framebuffer::new(30, 30);
    fb.clear((0, 0, 0));
    let lighting = Lighting::default();

    // Far blue sphere at z = 20.0
    draw_sphere(&mut fb, (15.0, 15.0, 20.0), 6.0, (0, 0, 255), &lighting);

    // Near red sphere at z = 10.0 overlapping center
    draw_sphere(&mut fb, (15.0, 15.0, 10.0), 6.0, (255, 0, 0), &lighting);

    // Center should be red (near sphere)
    let center_color = fb.get_pixel(15, 15).unwrap();
    assert!(center_color.0 > center_color.2); // red > blue
}

#[test]
fn test_rasterizer_draw_line_3d() {
    let mut fb = Framebuffer::new(20, 20);
    fb.clear((0, 0, 0));

    let p1 = (2.0, 2.0, 5.0);
    let p2 = (18.0, 2.0, 15.0);
    let line_color: PixelColor = (0, 255, 0);

    draw_line_3d(&mut fb, p1, p2, line_color);

    assert_eq!(fb.get_pixel(2, 2), Some(line_color));
    assert_eq!(fb.get_pixel(18, 2), Some(line_color));
    assert_eq!(fb.get_pixel(10, 2), Some(line_color));

    // Midpoint depth should be approx 10.0
    let mid_depth = fb.get_depth(10, 2).unwrap();
    assert!((mid_depth - 10.0).abs() < 1.0);
}

#[test]
fn test_rasterizer_draw_cylinder() {
    let mut fb = Framebuffer::new(30, 30);
    fb.clear((0, 0, 0));
    let lighting = Lighting::default();

    let p1 = (5.0, 15.0, 10.0);
    let p2 = (25.0, 15.0, 10.0);
    let cyl_color: PixelColor = (255, 200, 0);

    draw_cylinder(&mut fb, p1, p2, 3.0, cyl_color, &lighting);

    // Check pixels along the cylinder axis
    assert_ne!(fb.get_pixel(15, 15).unwrap(), (0, 0, 0));
    assert_ne!(fb.get_pixel(15, 13).unwrap(), (0, 0, 0));
    assert_ne!(fb.get_pixel(15, 17).unwrap(), (0, 0, 0));

    // Outside the cylinder thickness should remain background
    assert_eq!(fb.get_pixel(15, 5).unwrap(), (0, 0, 0));
}

#[test]
fn test_rasterizer_draw_triangle_3d() {
    let mut fb = Framebuffer::new(30, 30);
    fb.clear((0, 0, 0));
    let lighting = Lighting::default();

    let v1 = (15.0, 5.0, 10.0);
    let v2 = (5.0, 25.0, 10.0);
    let v3 = (25.0, 25.0, 10.0);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let tri_color: PixelColor = (100, 150, 250);

    draw_triangle_3d(&mut fb, v1, v2, v3, normal, tri_color, &lighting);

    // Centroid of the triangle (15, 18) should be filled
    assert_ne!(fb.get_pixel(15, 18).unwrap(), (0, 0, 0));
    assert!((fb.get_depth(15, 18).unwrap() - 10.0).abs() < 0.5);

    // Outside triangle (e.g. top corner 0, 0) should remain empty
    assert_eq!(fb.get_pixel(0, 0).unwrap(), (0, 0, 0));
}

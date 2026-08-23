use termpdb::render::buffer::Framebuffer;
use termpdb::render::lighting::Lighting;
use termpdb::render::rasterizer::{draw_cylinder, draw_sphere};

#[test]
fn test_offscreen_sphere_culling() {
    let mut fb = Framebuffer::new(40, 40);
    fb.clear((0, 0, 0));

    let lighting = Lighting::default();

    // Far off-screen center at (-50, -50), radius 10
    draw_sphere(&mut fb, (-50.0, -50.0, 5.0), 10.0, (255, 0, 0), &lighting);

    // Far off-screen center at (100, 100), radius 10
    draw_sphere(&mut fb, (100.0, 100.0, 5.0), 10.0, (0, 255, 0), &lighting);

    // Framebuffer should remain completely untouched
    assert!(fb.pixels.iter().all(|&p| p == (0, 0, 0)));
    assert!(fb.depth.iter().all(|&d| d.is_infinite()));
}

#[test]
fn test_onscreen_sphere_rasterization() {
    let mut fb = Framebuffer::new(40, 40);
    fb.clear((0, 0, 0));

    let lighting = Lighting::default();

    // Draw sphere centered at (20, 20) with radius 8.0 and depth 10.0
    draw_sphere(&mut fb, (20.0, 20.0, 10.0), 8.0, (200, 100, 50), &lighting);

    // Center pixel (20, 20) with sample point (20.5, 20.5) depth is close to cz - r = 2.0
    let center_depth = fb.get_depth(20, 20).unwrap();
    assert!((2.0..2.05).contains(&center_depth));
    assert_ne!(fb.get_pixel(20, 20).unwrap(), (0, 0, 0));

    // Pixel far outside radius should remain empty
    assert_eq!(fb.get_pixel(0, 0).unwrap(), (0, 0, 0));
    assert!(fb.get_depth(0, 0).unwrap().is_infinite());
}

#[test]
fn test_cylinder_bounding_clipping() {
    let mut fb = Framebuffer::new(30, 30);
    fb.clear((0, 0, 0));

    let lighting = Lighting::default();

    // Draw 3D cylinder from (5, 5, 10) to (25, 25, 10) with radius 2.0
    draw_cylinder(
        &mut fb,
        (5.0, 5.0, 10.0),
        (25.0, 25.0, 10.0),
        2.0,
        (100, 200, 255),
        &lighting,
    );

    // Midpoint (15, 15) should be drawn
    assert_ne!(fb.get_pixel(15, 15).unwrap(), (0, 0, 0));
    assert!(fb.get_depth(15, 15).unwrap() < 12.0);
}

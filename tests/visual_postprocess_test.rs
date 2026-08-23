use termpdb::render::buffer::Framebuffer;
use termpdb::render::postprocess::{apply_postprocessing, PostProcessConfig};

#[test]
fn test_silhouette_outline_detects_edges() {
    let mut fb = Framebuffer::new(10, 10);
    fb.clear((0, 0, 0));

    // Place a 4x4 square at depth 5.0 with color (200, 200, 200)
    for y in 3i32..7 {
        for x in 3i32..7 {
            fb.set_pixel(x, y, 5.0, (200, 200, 200));
        }
    }

    let config = PostProcessConfig {
        outline: true,
        ssao: false,
        outline_threshold: 0.1,
        ssao_radius: 2,
    };

    apply_postprocessing(&mut fb, &config);

    // Border of square (3, 3) should be darkened because neighbors outside are INFINITY
    let border_color = fb.get_pixel(3, 3).unwrap();
    // Center of square (4, 4) or (5, 5) surrounded by depth 5.0 should remain original color
    let center_color = fb.get_pixel(4, 4).unwrap();

    assert!(border_color.0 < center_color.0, "Border pixel should be darkened as an outline");
    assert_eq!(center_color, (200, 200, 200), "Interior pixel should retain full brightness");
}

#[test]
fn test_ssao_darkens_cavities() {
    let mut fb = Framebuffer::new(10, 10);
    fb.clear((0, 0, 0));

    // Create a flat surface at depth 10.0
    for y in 0i32..10 {
        for x in 0i32..10 {
            fb.set_pixel(x, y, 10.0, (180, 180, 180));
        }
    }

    // Create a "wall" around pixel (5, 5) with closer depth 8.0
    for dy in -2isize..=2 {
        for dx in -2isize..=2 {
            if dx != 0 || dy != 0 {
                let x = (5 + dx) as i32;
                let y = (5 + dy) as i32;
                fb.set_pixel(x, y, 8.0, (180, 180, 180));
            }
        }
    }

    let config = PostProcessConfig {
        outline: false,
        ssao: true,
        outline_threshold: 0.1,
        ssao_radius: 2,
    };

    apply_postprocessing(&mut fb, &config);

    let cavity_pixel = fb.get_pixel(5, 5).unwrap();
    // Cavity pixel (5,5) at depth 10.0 surrounded by walls at depth 8.0 should be occluded & darkened
    assert!(cavity_pixel.0 < 180, "Cavity pixel should have lower brightness due to SSAO");
}

//! Guard-rail tests for hostile render inputs.
//!
//! These pin the two bounds that keep extreme `--width`/`--height`/`--ssaa`
//! flags and far-off-screen geometry from turning into OOM aborts or
//! multi-second rasterizer hangs: the framebuffer pixel budget and
//! viewport-clipped line stepping.

use std::time::Instant;

use termpdb::math::Vec3;
use termpdb::model::{Atom, Element, Structure};
use termpdb::render::buffer::MAX_FRAMEBUFFER_PIXELS;
use termpdb::render::{
    BrailleBuffer, ColorScheme, Framebuffer, LodMode, PixelColor, RenderMode, Visibility,
    clip_segment_to_screen, draw_dashed_line_3d, draw_line_3d, draw_overlay_line,
    export_ansi_with_visibility, fit_render_size, render_supersampled,
};

fn tiny_structure() -> Structure {
    let mut structure = Structure::new("");
    let carbon = Element {
        atomic_number: 6,
        symbol: "C",
        name: "Carbon",
        covalent_radius: 0.77,
        vdw_radius: 1.70,
        cpk_color: (144, 144, 144),
    };
    for i in 0..4 {
        structure.add_atom(Atom::new(
            i,
            i as i32 + 1,
            "CA",
            carbon,
            Vec3::new(i as f32, (i % 2) as f32, 0.0),
            0.0,
            "ALA",
            i as i32 + 1,
            "A",
            false,
        ));
    }
    structure
}

#[test]
fn test_fit_render_size_leaves_normal_requests_unchanged() {
    assert_eq!(fit_render_size(800, 600, 2), (800, 600, 2));
    assert_eq!(fit_render_size(1920, 1080, 2), (1920, 1080, 2));
}

#[test]
fn test_fit_render_size_reduces_ssaa_before_dimensions() {
    // 1920x1080 base fits; ssaa is reduced until the product fits the budget.
    let (w, h, s) = fit_render_size(1920, 1080, 255);
    assert_eq!((w, h), (1920, 1080), "dimensions must be preserved");
    assert!(s < 255, "ssaa must be reduced");
    assert!(w * h * s * s <= MAX_FRAMEBUFFER_PIXELS);
}

#[test]
fn test_fit_render_size_scales_extreme_dimensions_proportionally() {
    let (w, h, s) = fit_render_size(65535, 65535, 255);
    assert_eq!(s, 1, "ssaa alone cannot rescue a huge base size");
    assert!(w * h <= MAX_FRAMEBUFFER_PIXELS);
    assert!((w > 0 && h > 0));
    // Aspect ratio preserved (both axes scaled by the same factor).
    let ratio = w as f64 / h as f64;
    assert!((ratio - 1.0).abs() < 0.01, "square input must stay square");
}

#[test]
fn test_fit_render_size_handles_zero_dimensions() {
    assert_eq!(fit_render_size(0, 600, 3), (0, 600, 3));
    assert_eq!(fit_render_size(800, 0, 3), (800, 0, 3));
}

#[test]
fn test_supersampled_export_with_absurd_flags_is_bounded() {
    // Previously requested a ~terabyte framebuffer and aborted with OOM.
    let structure = tiny_structure();
    let rgba = render_supersampled(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        u16::MAX as usize,
        u16::MAX as usize,
        255,
        Visibility::default(),
        LodMode::Auto,
    );
    assert_eq!(rgba.len() % 4, 0);
    assert!(
        rgba.len() <= MAX_FRAMEBUFFER_PIXELS * 4,
        "output buffer must respect the pixel budget"
    );
}

#[test]
fn test_ansi_export_with_absurd_flags_is_bounded() {
    // 65535 columns x 65535 rows = 8.6 Gpx previously; must now clamp.
    let structure = tiny_structure();
    let ansi = export_ansi_with_visibility(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        u16::MAX,
        u16::MAX,
        Visibility::default(),
        LodMode::Auto,
    );
    assert!(!ansi.is_empty());
}

/// Draws a line whose endpoints are tens of millions of pixels off-screen but
/// which crosses the viewport. Returns lit-pixel count and elapsed time.
fn draw_far_offscreen_line(draw: impl FnOnce(&mut Framebuffer)) -> (usize, std::time::Duration) {
    let mut fb = Framebuffer::new(120, 80);
    let start = Instant::now();
    draw(&mut fb);
    let elapsed = start.elapsed();
    let lit = fb.pixels.iter().filter(|&&p| p != (0, 0, 0)).count();
    (lit, elapsed)
}

#[test]
fn test_line_steppers_survive_extreme_screen_coordinates() {
    const GRAY: PixelColor = (200, 200, 200);

    // Regime 1 -- realistic pass-through magnitude (~10^7 px): the clipped
    // line must preserve its geometry, lighting roughly one pixel per step
    // along the viewport diagonal.
    let p1 = (3.0e7, 2.0e7, 1.0);
    let p2 = (-3.0e7, -2.0e7, 2.0);

    let (lit, elapsed) = draw_far_offscreen_line(|fb| draw_line_3d(fb, p1, p2, GRAY));
    assert!(
        (10..=300).contains(&lit),
        "clipped diagonal should light a bounded number of pixels, got {lit}"
    );
    assert!(elapsed.as_secs() < 5, "draw_line_3d hung for {elapsed:?}");

    // Regime 2 -- absurd magnitude (~10^9 px): f32 precision legitimately
    // degrades the endpoints here, but stepping must terminate instantly
    // instead of iterating billions of times.
    let p1 = (3.0e9, 2.0e9, 1.0);
    let p2 = (-3.0e9, -2.0e9, 2.0);

    let (lit, elapsed) = draw_far_offscreen_line(|fb| draw_line_3d(fb, p1, p2, GRAY));
    assert!(lit <= 300, "unbounded stepping lit {lit} pixels");
    assert!(elapsed.as_secs() < 5, "draw_line_3d hung for {elapsed:?}");

    let (_, elapsed) = draw_far_offscreen_line(|fb| {
        draw_dashed_line_3d(fb, p1, p2, GRAY, 4.0, 2.0);
    });
    assert!(
        elapsed.as_secs() < 5,
        "draw_dashed_line_3d hung for {elapsed:?}"
    );

    let (_, elapsed) = draw_far_offscreen_line(|fb| draw_overlay_line(fb, p1, p2, GRAY));
    assert!(
        elapsed.as_secs() < 5,
        "draw_overlay_line hung for {elapsed:?}"
    );

    // Braille canvas: subpixel space, same coordinate magnitudes.
    let mut canvas = BrailleBuffer::new(240, 160);
    let start = Instant::now();
    canvas.draw_line_3d(p1, p2, GRAY);
    let braille_elapsed = start.elapsed();
    assert!(
        braille_elapsed.as_secs() < 5,
        "braille draw_line_3d hung for {braille_elapsed:?}"
    );

    // Fully off-screen segment clips to None without touching anything.
    let mut fb = Framebuffer::new(120, 80);
    draw_line_3d(&mut fb, (-1.0e7, -1.0e7, 0.0), (-2.0e7, -2.0e7, 1.0), GRAY);
    assert!(fb.pixels.iter().all(|&p| p == (0, 0, 0)));
}

#[test]
fn test_clip_segment_to_screen_basic_geometry() {
    // Fully inside: unchanged endpoints (within float epsilon).
    let clipped = clip_segment_to_screen((10.0, 10.0, 1.0), (50.0, 30.0, 2.0), 100, 100);
    let Some(((ax, ay, _), (bx, by, _))) = clipped else {
        panic!("visible segment must not be clipped away");
    };
    assert!((ax - 10.0).abs() < 1e-4 && (ay - 10.0).abs() < 1e-4);
    assert!((bx - 50.0).abs() < 1e-4 && (by - 30.0).abs() < 1e-4);

    // Fully outside one corner: rejected.
    assert!(clip_segment_to_screen((-5.0, -5.0, 0.0), (-9.0, -9.0, 0.0), 100, 100).is_none());

    // Crossing the left edge with the far endpoint inside: only the left end
    // is pulled onto x = 0; the right endpoint stays put.
    let clipped = clip_segment_to_screen((-40.0, 50.0, 1.0), (60.0, 50.0, 2.0), 100, 100);
    let Some(((ax, ay, az), (bx, by, bz))) = clipped else {
        panic!("crossing segment must survive clipping");
    };
    assert_eq!(ax, 0.0);
    assert_eq!(ay, 50.0);
    assert_eq!(bx, 60.0);
    assert_eq!(by, 50.0);
    // Depth interpolates linearly: t0 = 40/100 shifts z from 1.0 to 1.4.
    assert!((az - 1.4).abs() < 1e-4);
    assert_eq!(bz, 2.0);
}

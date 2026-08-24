//! Near-plane behavior tests for spheres and cylinders.
//!
//! Regression guards for the wholesale primitive culling bug: primitives whose
//! reference point is nearer than the camera near plane used to vanish
//! entirely, making space-filling molecules pop out of existence while
//! zooming. Cylinders must now clip in view space and spheres must survive
//! while their surface crosses the near plane.

use termpdb::math::Vec3;
use termpdb::model::{Atom, Element, Structure};
use termpdb::render::{
    Camera, ColorScheme, Framebuffer, Lighting, LodMode, RenderMode, Visibility, render_structure,
};

/// Camera aimed at the origin from +Z, `cam.distance` away, with the near
/// plane pushed past the target so geometry at the origin is "behind near".
fn camera_with_near(distance: f32, near: f32) -> Camera {
    let mut cam = Camera::new();
    cam.target = Vec3::ZERO;
    cam.orientation = termpdb::math::Quat::identity();
    cam.distance = distance;
    cam.aspect = 80.0 / 48.0;
    cam.near = near;
    cam
}

#[test]
fn test_project_segment_clips_against_near_plane() {
    let cam = camera_with_near(10.0, 12.0);
    let mats = cam.matrices();

    // Points ON the view axis at distances 7 and 14 (camera at +10 looking
    // toward -Z: depth = 10 - z).
    let a = Vec3::new(0.0, 0.0, 3.0); // depth 7  (behind near)
    let b = Vec3::new(0.0, 0.0, -4.0); // depth 14 (in front)

    let result = cam.project_segment(&mats, a, b, 80, 48);
    let Some((pa, pb)) = result else {
        panic!("segment straddling the near plane must survive clipping");
    };
    // The behind-near endpoint is clipped exactly onto the near plane...
    assert!(
        (pa.2 - 12.0).abs() < 1e-4,
        "clipped depth {}, want 12",
        pa.2
    );
    // ...and the front endpoint keeps its true depth.
    assert!((pb.2 - 14.0).abs() < 1e-4);
}

#[test]
fn test_project_segment_fully_behind_near_is_culled() {
    let cam = camera_with_near(10.0, 12.0);
    let mats = cam.matrices();
    // Depths 5 and 6: entirely between camera and near plane.
    let a = Vec3::new(0.0, 0.0, 5.0);
    let b = Vec3::new(0.0, 0.0, 4.0);
    assert!(cam.project_segment(&mats, a, b, 80, 48).is_none());
}

#[test]
fn test_project_segment_fully_in_front_is_unchanged() {
    let cam = camera_with_near(10.0, 12.0);
    let mats = cam.matrices();
    let a = Vec3::new(0.0, 0.0, -3.0); // depth 13
    let b = Vec3::new(0.0, 0.0, -4.0); // depth 14
    let Some((pa, pb)) = cam.project_segment(&mats, a, b, 80, 48) else {
        panic!("fully visible segment must not be clipped");
    };
    assert!((pa.2 - 13.0).abs() < 1e-4 && (pb.2 - 14.0).abs() < 1e-4);
}

#[test]
fn test_project_sphere_crossing_near_plane_still_draws() {
    let cam = camera_with_near(10.0, 12.0);
    let mats = cam.matrices();

    // Center at the target (depth 10 < near 12), radius 5: the surface
    // reaches depth 15, well past the near plane.
    let result = cam.project_sphere(&mats, Vec3::ZERO, 5.0, 80, 48);
    let Some(((sx, sy, depth), screen_r)) = result else {
        panic!("sphere crossing the near plane must not be culled");
    };
    assert_eq!(depth, 12.0, "z-buffer depth clamps up to the near plane");
    assert!(screen_r > 1.0, "silhouette must remain substantial");
    let _ = (sx, sy);
}

#[test]
fn test_project_sphere_entirely_before_near_is_culled() {
    let cam = camera_with_near(10.0, 12.0);
    let mats = cam.matrices();
    // Center depth 10 + radius 1 <= near 12: nothing visible remains.
    assert!(cam.project_sphere(&mats, Vec3::ZERO, 1.0, 80, 48).is_none());
}

#[test]
fn test_vdw_render_with_center_behind_near_draws_surface() {
    // Integration regression: a single big atom sitting at the camera target
    // with the near plane pushed past it. The pre-fix code culled the atom
    // outright and produced an empty frame.
    let mut structure = Structure::new("nearplane");
    structure.add_atom(Atom::new(
        0,
        1,
        "XX",
        Element {
            symbol: "Xx",
            name: "Testium",
            atomic_number: 0,
            covalent_radius: 1.0,
            vdw_radius: 5.0,
            cpk_color: (200, 60, 60),
        },
        Vec3::ZERO,
        0.0,
        "XX",
        1,
        "A",
        false,
    ));

    let cam = camera_with_near(10.0, 12.0);
    let fb_size = (80usize, 48usize);
    let mut framebuffer = Framebuffer::new(fb_size.0, fb_size.1 * 2);
    framebuffer.clear((0, 0, 0));
    render_structure(
        &structure,
        RenderMode::Vdw,
        ColorScheme::Cpk,
        &cam,
        &mut framebuffer,
        &Lighting::default(),
        Visibility::default(),
        LodMode::Full,
    );

    let lit = framebuffer
        .pixels
        .iter()
        .filter(|&&p| p != (0, 0, 0))
        .count();
    assert!(
        lit > 100,
        "atom surface crossing the near plane must rasterize; got {lit} lit pixels"
    );
}

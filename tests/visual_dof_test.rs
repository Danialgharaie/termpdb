use termpdb::math::Vec3;
use termpdb::render::Lighting;

#[test]
fn test_dof_focal_plane_attenuation() {
    let mut lighting = Lighting::default();
    lighting.dof_focus = Some(5.0); // Focal distance at depth 5.0
    lighting.dof_range = 2.0;       // In-focus depth radius +-2.0

    let base_color = (200, 200, 200);
    let normal = Vec3::new(0.0, 0.0, 1.0);

    // At focal plane depth (5.0), object is fully in focus
    let in_focus = lighting.shade(normal, 5.0, base_color, 0.0, 20.0);

    // Far from focal plane (depth 18.0), object is attenuated
    let out_of_focus = lighting.shade(normal, 18.0, base_color, 0.0, 20.0);

    assert!(in_focus.0 > out_of_focus.0, "In-focus object should be brighter/clearer than out-of-focus object");
}

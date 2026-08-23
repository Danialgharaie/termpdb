use termpdb::math::Vec3;
use termpdb::render::Lighting;

#[test]
fn test_specular_shading_and_highlights() {
    let lighting = Lighting::new(
        Vec3::new(0.0, 0.0, 1.0), // Light pointing directly from viewer
        0.2,                      // Ambient
        0.8,                      // Diffuse
        0.5,                      // Specular
        16.0,                     // Shininess
        0.0,                      // No fog
    );

    let base_color = (100, 100, 100);
    // Normal pointing directly to viewer (+Z) aligns with H
    let normal_facing = Vec3::new(0.0, 0.0, 1.0);
    let shaded_facing = lighting.shade(normal_facing, 0.0, base_color, 0.0, 10.0);

    // Normal pointing 90 deg away (+X)
    let normal_away = Vec3::new(1.0, 0.0, 0.0);
    let shaded_away = lighting.shade(normal_away, 0.0, base_color, 0.0, 10.0);

    // Shaded facing should have specular glint added to it
    assert!(shaded_facing.0 > shaded_away.0);
    assert!(shaded_facing.0 > 100); // 100 * 1.0 + 255 * 0.5 ~ 228
}

#[test]
fn test_specular_zero_intensity() {
    let lighting = Lighting::new(
        Vec3::new(0.0, 0.0, 1.0),
        0.2,
        0.8,
        0.0, // Zero specular
        16.0,
        0.0,
    );

    let base = (100, 100, 100);
    let facing = lighting.shade(Vec3::new(0.0, 0.0, 1.0), 0.0, base, 0.0, 10.0);
    // (0.2 + 0.8 * 1.0) * 100 = 100
    assert_eq!(facing, (100, 100, 100));
}

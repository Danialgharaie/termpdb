use termpdb::math::Vec3;
use termpdb::math::kabsch::{kabsch_align, svd_3x3};

#[test]
fn test_svd_3x3_orthogonal() {
    // Identity matrix
    let a = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let (u, s, v) = svd_3x3(a);

    // Singular values should be 1.0
    assert!((s[0] - 1.0).abs() < 1e-4);
    assert!((s[1] - 1.0).abs() < 1e-4);
    assert!((s[2] - 1.0).abs() < 1e-4);

    // U * S * V^T should equal A
    let mut recon = [[0.0f32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                recon[i][j] += u[i][k] * s[k] * v[j][k];
            }
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            assert!((recon[i][j] - a[i][j]).abs() < 1e-4);
        }
    }
}

#[test]
fn test_kabsch_align_identity_and_pure_translation() {
    let p = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];

    // Identity test
    let res_id = kabsch_align(&p, &p).expect("Kabsch should align identical sets");
    assert!(res_id.rmsd < 1e-5, "RMSD of identical sets should be ~0");

    // Pure translation by (5, -3, 2)
    let shift = Vec3::new(5.0, -3.0, 2.0);
    let q: Vec<Vec3> = p.iter().map(|&v| v + shift).collect();

    let res_shift = kabsch_align(&p, &q).expect("Kabsch should align translated sets");
    assert!(
        res_shift.rmsd < 1e-5,
        "RMSD of pure translation should be ~0"
    );
    assert!((res_shift.translation.x - 5.0).abs() < 1e-4);
    assert!((res_shift.translation.y - (-3.0)).abs() < 1e-4);
    assert!((res_shift.translation.z - 2.0).abs() < 1e-4);
}

#[test]
fn test_kabsch_align_90deg_rotation() {
    let p = vec![
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(2.0, 1.0, 0.0),
    ];

    // Rotate 90 degrees around Z axis: (x, y, z) -> (-y, x, z) + translation (1, 1, 1)
    let q: Vec<Vec3> = p
        .iter()
        .map(|v| Vec3::new(-v.y + 1.0, v.x + 1.0, v.z + 1.0))
        .collect();

    let res = kabsch_align(&p, &q).expect("Kabsch should align 90-deg rotated sets");
    assert!(res.rmsd < 1e-4, "RMSD should be ~0 for rigid rotation");

    // Check transformed points match target
    for i in 0..p.len() {
        let transformed = res.transform_point(p[i]);
        assert!(transformed.distance(&q[i]) < 1e-3);
    }
}

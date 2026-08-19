use std::f32::consts::PI;
use termpdb::error::{Result, TermPdbError};
use termpdb::math::{CatmullRomSpline, Mat4, Quat, Vec3};

const EPSILON: f32 = 1e-4;

fn approx_eq(a: f32, b: f32) -> bool {
    (a - b).abs() < EPSILON
}

fn vec3_approx_eq(a: Vec3, b: Vec3) -> bool {
    approx_eq(a.x, b.x) && approx_eq(a.y, b.y) && approx_eq(a.z, b.z)
}

#[test]
fn test_vec3_basic_construction_and_index() {
    let v = Vec3::new(1.0, 2.0, 3.0);
    assert_eq!(v.x, 1.0);
    assert_eq!(v.y, 2.0);
    assert_eq!(v.z, 3.0);

    assert_eq!(v[0], 1.0);
    assert_eq!(v[1], 2.0);
    assert_eq!(v[2], 3.0);

    let mut vm = v;
    vm[0] = 10.0;
    vm[1] = 20.0;
    vm[2] = 30.0;
    assert_eq!(vm, Vec3::new(10.0, 20.0, 30.0));

    let zero = Vec3::zero();
    assert_eq!(zero, Vec3::new(0.0, 0.0, 0.0));

    let splat = Vec3::splat(5.0);
    assert_eq!(splat, Vec3::new(5.0, 5.0, 5.0));
}

#[test]
fn test_vec3_arithmetic_operators() {
    let a = Vec3::new(1.0, 2.0, 3.0);
    let b = Vec3::new(4.0, 5.0, 6.0);

    assert_eq!(a + b, Vec3::new(5.0, 7.0, 9.0));
    assert_eq!(a - b, Vec3::new(-3.0, -3.0, -3.0));
    assert_eq!(a * 2.0, Vec3::new(2.0, 4.0, 6.0));
    assert_eq!(2.0 * a, Vec3::new(2.0, 4.0, 6.0));
    assert_eq!(b / 2.0, Vec3::new(2.0, 2.5, 3.0));
    assert_eq!(-a, Vec3::new(-1.0, -2.0, -3.0));
    assert_eq!(a * b, Vec3::new(4.0, 10.0, 18.0));

    let mut c = a;
    c += b;
    assert_eq!(c, Vec3::new(5.0, 7.0, 9.0));

    let mut d = b;
    d -= a;
    assert_eq!(d, Vec3::new(3.0, 3.0, 3.0));

    let mut e = a;
    e *= 3.0;
    assert_eq!(e, Vec3::new(3.0, 6.0, 9.0));

    let mut f = Vec3::new(6.0, 8.0, 10.0);
    f /= 2.0;
    assert_eq!(f, Vec3::new(3.0, 4.0, 5.0));
}

#[test]
fn test_vec3_dot_and_cross() {
    let x = Vec3::new(1.0, 0.0, 0.0);
    let y = Vec3::new(0.0, 1.0, 0.0);
    let z = Vec3::new(0.0, 0.0, 1.0);

    // Dot products
    assert_eq!(x.dot(y), 0.0);
    assert_eq!(x.dot(x), 1.0);
    assert_eq!(
        Vec3::new(1.0, 2.0, 3.0).dot(Vec3::new(4.0, -5.0, 6.0)),
        12.0
    );

    // Cross products (Right hand rule)
    assert_eq!(x.cross(y), z);
    assert_eq!(y.cross(z), x);
    assert_eq!(z.cross(x), y);
    assert_eq!(y.cross(x), -z);
    assert_eq!(x.cross(x), Vec3::zero());
}

#[test]
fn test_vec3_norm_and_normalize() {
    let v = Vec3::new(3.0, 4.0, 0.0);
    assert_eq!(v.norm_squared(), 25.0);
    assert_eq!(v.norm(), 5.0);

    let unit = v.normalize();
    assert!(vec3_approx_eq(unit, Vec3::new(0.6, 0.8, 0.0)));
    assert!(approx_eq(unit.norm(), 1.0));

    // Zero vector normalization should be safe
    let zero = Vec3::zero();
    assert_eq!(zero.normalize(), Vec3::zero());
}

#[test]
fn test_vec3_lerp_and_distance() {
    let a = Vec3::new(0.0, 0.0, 0.0);
    let b = Vec3::new(10.0, 20.0, 30.0);

    assert_eq!(a.lerp(b, 0.0), a);
    assert_eq!(a.lerp(b, 1.0), b);
    assert_eq!(a.lerp(b, 0.5), Vec3::new(5.0, 10.0, 15.0));

    assert_eq!(a.distance(&b), (100.0 + 400.0 + 900.0_f32).sqrt());
    assert_eq!(a.distance_squared(&b), 1400.0);

    let v1 = Vec3::new(1.0, 5.0, 2.0);
    let v2 = Vec3::new(3.0, 2.0, 4.0);
    assert_eq!(v1.min(v2), Vec3::new(1.0, 2.0, 2.0));
    assert_eq!(v1.max(v2), Vec3::new(3.0, 5.0, 4.0));
}

#[test]
fn test_mat4_identity_and_translation() {
    let id = Mat4::identity();
    let p = Vec3::new(3.0, 4.0, 5.0);
    let v = Vec3::new(1.0, 0.0, -1.0);

    assert_eq!(id.transform_point(p), p);
    assert_eq!(id.transform_vector(v), v);
    assert_eq!(id * id, id);

    let trans = Mat4::from_translation(Vec3::new(10.0, -5.0, 2.0));
    assert_eq!(trans.transform_point(p), Vec3::new(13.0, -1.0, 7.0));
    // Vectors are invariant under translation
    assert_eq!(trans.transform_vector(v), v);
}

#[test]
fn test_mat4_scale_and_rotation() {
    let scale = Mat4::from_scale(Vec3::new(2.0, 0.5, 3.0));
    let p = Vec3::new(1.0, 4.0, 2.0);
    assert_eq!(scale.transform_point(p), Vec3::new(2.0, 2.0, 6.0));
    assert_eq!(scale.transform_vector(p), Vec3::new(2.0, 2.0, 6.0));

    // Rotate 90 degrees around Z axis: (1, 0, 0) -> (0, 1, 0)
    let rot_z = Mat4::from_rotation_z(PI / 2.0);
    let rx = rot_z.transform_point(Vec3::new(1.0, 0.0, 0.0));
    assert!(vec3_approx_eq(rx, Vec3::new(0.0, 1.0, 0.0)));

    // Rotate 90 degrees around X axis: (0, 1, 0) -> (0, 0, 1)
    let rot_x = Mat4::from_rotation_x(PI / 2.0);
    let ry = rot_x.transform_point(Vec3::new(0.0, 1.0, 0.0));
    assert!(vec3_approx_eq(ry, Vec3::new(0.0, 0.0, 1.0)));

    // Rotate 90 degrees around Y axis: (0, 0, 1) -> (1, 0, 0)
    let rot_y = Mat4::from_rotation_y(PI / 2.0);
    let rz = rot_y.transform_point(Vec3::new(0.0, 0.0, 1.0));
    assert!(vec3_approx_eq(rz, Vec3::new(1.0, 0.0, 0.0)));
}

#[test]
fn test_mat4_multiplication_associativity() {
    let t = Mat4::from_translation(Vec3::new(1.0, 2.0, 3.0));
    let r = Mat4::from_rotation_z(PI / 4.0);
    let s = Mat4::from_scale(Vec3::splat(2.0));

    let combined = t * r * s;
    let p = Vec3::new(1.0, 0.0, 0.0);

    let res1 = combined.transform_point(p);
    let res2 = t.transform_point(r.transform_point(s.transform_point(p)));
    assert!(vec3_approx_eq(res1, res2));
}

#[test]
fn test_mat4_look_at() {
    let eye = Vec3::new(0.0, 0.0, 10.0);
    let target = Vec3::new(0.0, 0.0, 0.0);
    let up = Vec3::new(0.0, 1.0, 0.0);

    let view = Mat4::look_at(eye, target, up);

    // Eye transforms to origin in view space
    let eye_view = view.transform_point(eye);
    assert!(vec3_approx_eq(eye_view, Vec3::zero()));

    // Target transforms to (0, 0, -10)
    let target_view = view.transform_point(target);
    assert!(vec3_approx_eq(target_view, Vec3::new(0.0, 0.0, -10.0)));
}

#[test]
fn test_mat4_perspective_and_orthographic() {
    let proj_persp = Mat4::perspective(PI / 2.0, 1.0, 1.0, 100.0);
    // Point on near plane along -Z: (0, 0, -1) -> NDC z = -1
    let near_pt = proj_persp.transform_point(Vec3::new(0.0, 0.0, -1.0));
    assert!(approx_eq(near_pt.z, -1.0));

    // Point on far plane along -Z: (0, 0, -100) -> NDC z = 1
    let far_pt = proj_persp.transform_point(Vec3::new(0.0, 0.0, -100.0));
    assert!(approx_eq(far_pt.z, 1.0));

    let proj_ortho = Mat4::orthographic(-10.0, 10.0, -5.0, 5.0, 1.0, 100.0);
    let left_bottom_near = proj_ortho.transform_point(Vec3::new(-10.0, -5.0, -1.0));
    assert!(approx_eq(left_bottom_near.x, -1.0));
    assert!(approx_eq(left_bottom_near.y, -1.0));
    assert!(approx_eq(left_bottom_near.z, -1.0));

    let right_top_far = proj_ortho.transform_point(Vec3::new(10.0, 5.0, -100.0));
    assert!(approx_eq(right_top_far.x, 1.0));
    assert!(approx_eq(right_top_far.y, 1.0));
    assert!(approx_eq(right_top_far.z, 1.0));
}

#[test]
fn test_mat4_transpose_and_inverse() {
    let m = Mat4::from_translation(Vec3::new(3.0, -2.0, 5.0)) * Mat4::from_rotation_y(0.7);
    let inv = m.inverse().expect("matrix should be invertible");
    let identity = m * inv;

    let p = Vec3::new(7.0, 8.0, 9.0);
    let p_transformed = identity.transform_point(p);
    assert!(vec3_approx_eq(p, p_transformed));
}

#[test]
fn test_quat_identity_and_axis_angle() {
    let q_id = Quat::identity();
    let v = Vec3::new(1.0, 2.0, 3.0);
    assert!(vec3_approx_eq(q_id.rotate_vec3(v), v));

    // 90 deg rotation around Z
    let q_z90 = Quat::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), PI / 2.0);
    let v_rot = q_z90.rotate_vec3(Vec3::new(1.0, 0.0, 0.0));
    assert!(vec3_approx_eq(v_rot, Vec3::new(0.0, 1.0, 0.0)));

    // 180 deg rotation around Y
    let q_y180 = Quat::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), PI);
    let v_rot_y = q_y180.rotate_vec3(Vec3::new(1.0, 0.0, 0.0));
    assert!(vec3_approx_eq(v_rot_y, Vec3::new(-1.0, 0.0, 0.0)));
}

#[test]
fn test_quat_to_mat4_equivalence() {
    let q = Quat::from_axis_angle(Vec3::new(1.0, 2.0, 3.0).normalize(), 1.25);
    let m = q.to_mat4();

    let p = Vec3::new(4.0, -2.0, 5.0);
    let from_quat = q.rotate_vec3(p);
    let from_mat = m.transform_vector(p);

    assert!(vec3_approx_eq(from_quat, from_mat));
}

#[test]
fn test_quat_multiplication_and_inverse() {
    let q1 = Quat::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), PI / 2.0);
    let q2 = Quat::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), PI / 2.0);

    let q_comb = q1 * q2;
    let v = Vec3::new(0.0, 0.0, 1.0);

    // (q1 * q2) * v == q1 * (q2 * v)
    let res1 = q_comb.rotate_vec3(v);
    let res2 = q1.rotate_vec3(q2.rotate_vec3(v));
    assert!(vec3_approx_eq(res1, res2));

    let q_inv = q_comb.inverse();
    let restored = q_inv.rotate_vec3(res1);
    assert!(vec3_approx_eq(restored, v));
}

#[test]
fn test_quat_slerp_and_drag() {
    let q0 = Quat::identity();
    let q1 = Quat::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), PI / 2.0);

    let mid = q0.slerp(&q1, 0.5);
    let expected = Quat::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), PI / 4.0);

    let v = Vec3::new(1.0, 0.0, 0.0);
    assert!(vec3_approx_eq(mid.rotate_vec3(v), expected.rotate_vec3(v)));

    // Drag produces valid rotation quaternion
    let q_drag = Quat::from_drag(10.0, 5.0, 0.01);
    assert!(approx_eq(q_drag.norm(), 1.0));
}

#[test]
fn test_spline_interpolation_and_continuity() {
    let points = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 2.0, 0.0),
        Vec3::new(2.0, 0.0, 0.0),
        Vec3::new(3.0, 2.0, 0.0),
    ];

    let spline = CatmullRomSpline::new(points.clone());

    // Interpolate at endpoints of interior segments
    // For 4 points (indices 0, 1, 2, 3), interior points are at t values corresponding to points[1] and points[2]
    let start = spline.interpolate(0.0);
    assert!(vec3_approx_eq(start, points[0]));

    let end = spline.interpolate(1.0);
    assert!(vec3_approx_eq(end, points[points.len() - 1]));

    // Tangent should be continuous and non-zero
    let t_mid = spline.tangent(0.5);
    assert!(t_mid.norm() > 0.0);
    assert!(approx_eq(t_mid.norm(), 1.0));

    // Smooth curve generation
    let samples = spline.generate_smooth_curve(10);
    assert!(!samples.is_empty());
    assert!(vec3_approx_eq(samples[0], points[0]));
    assert!(vec3_approx_eq(
        *samples.last().unwrap(),
        *points.last().unwrap()
    ));
}

#[test]
fn test_spline_edge_cases() {
    // Empty points
    let empty_spline = CatmullRomSpline::new(vec![]);
    assert_eq!(empty_spline.interpolate(0.5), Vec3::zero());
    assert!(empty_spline.generate_smooth_curve(10).is_empty());

    // Single point
    let single_point = vec![Vec3::new(5.0, 5.0, 5.0)];
    let single_spline = CatmullRomSpline::new(single_point.clone());
    assert_eq!(single_spline.interpolate(0.0), single_point[0]);
    assert_eq!(single_spline.interpolate(1.0), single_point[0]);

    // Two points
    let two_points = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0)];
    let two_spline = CatmullRomSpline::new(two_points.clone());
    let mid = two_spline.interpolate(0.5);
    assert!(vec3_approx_eq(mid, Vec3::new(5.0, 0.0, 0.0)));
}

#[test]
fn test_error_types_and_result() {
    fn make_parse_error() -> Result<()> {
        Err(TermPdbError::ParseError("invalid ATOM record".into()))
    }

    fn make_io_error() -> Result<()> {
        Err(TermPdbError::Io(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "file not found",
        )))
    }

    assert!(make_parse_error().is_err());
    let err_str = format!("{}", make_parse_error().unwrap_err());
    assert!(err_str.contains("invalid ATOM record"));

    assert!(make_io_error().is_err());
    let io_err_str = format!("{}", make_io_error().unwrap_err());
    assert!(io_err_str.contains("file not found"));
}

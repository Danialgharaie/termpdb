use crate::math::vec3::Vec3;
use std::ops::Mul;

/// A 4x4 matrix represented in column-major order.
///
/// Layout in memory:
/// m[0]  m[4]  m[8]   m[12]
/// m[1]  m[5]  m[9]   m[13]
/// m[2]  m[6]  m[10]  m[14]
/// m[3]  m[7]  m[11]  m[15]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat4 {
    pub m: [f32; 16],
}

impl Default for Mat4 {
    fn default() -> Self {
        Self::identity()
    }
}

impl Mat4 {
    pub const fn zeros() -> Self {
        Self { m: [0.0; 16] }
    }

    pub const fn identity() -> Self {
        let mut m = [0.0; 16];
        m[0] = 1.0;
        m[5] = 1.0;
        m[10] = 1.0;
        m[15] = 1.0;
        Self { m }
    }

    pub fn get(&self, row: usize, col: usize) -> f32 {
        self.m[col * 4 + row]
    }

    pub fn set(&mut self, row: usize, col: usize, val: f32) {
        self.m[col * 4 + row] = val;
    }

    pub fn from_translation(v: Vec3) -> Self {
        let mut m = Self::identity().m;
        m[12] = v.x;
        m[13] = v.y;
        m[14] = v.z;
        Self { m }
    }

    pub fn from_scale(s: Vec3) -> Self {
        let mut m = [0.0; 16];
        m[0] = s.x;
        m[5] = s.y;
        m[10] = s.z;
        m[15] = 1.0;
        Self { m }
    }

    pub fn from_rotation_x(rad: f32) -> Self {
        let (s, c) = rad.sin_cos();
        let mut m = Self::identity().m;
        m[5] = c;
        m[6] = s;
        m[9] = -s;
        m[10] = c;
        Self { m }
    }

    pub fn from_rotation_y(rad: f32) -> Self {
        let (s, c) = rad.sin_cos();
        let mut m = Self::identity().m;
        m[0] = c;
        m[2] = -s;
        m[8] = s;
        m[10] = c;
        Self { m }
    }

    pub fn from_rotation_z(rad: f32) -> Self {
        let (s, c) = rad.sin_cos();
        let mut m = Self::identity().m;
        m[0] = c;
        m[1] = s;
        m[4] = -s;
        m[5] = c;
        Self { m }
    }

    pub fn mul(&self, other: &Self) -> Self {
        let mut out = [0.0; 16];
        for col in 0..4 {
            for row in 0..4 {
                let mut sum = 0.0;
                for k in 0..4 {
                    sum += self.get(row, k) * other.get(k, col);
                }
                out[col * 4 + row] = sum;
            }
        }
        Self { m: out }
    }

    pub fn transform_point(&self, p: Vec3) -> Vec3 {
        let x = self.m[0] * p.x + self.m[4] * p.y + self.m[8] * p.z + self.m[12];
        let y = self.m[1] * p.x + self.m[5] * p.y + self.m[9] * p.z + self.m[13];
        let z = self.m[2] * p.x + self.m[6] * p.y + self.m[10] * p.z + self.m[14];
        let w = self.m[3] * p.x + self.m[7] * p.y + self.m[11] * p.z + self.m[15];

        if w != 0.0 && (w - 1.0).abs() > 1e-7 {
            let inv_w = 1.0 / w;
            Vec3::new(x * inv_w, y * inv_w, z * inv_w)
        } else {
            Vec3::new(x, y, z)
        }
    }

    pub fn transform_vector(&self, v: Vec3) -> Vec3 {
        let x = self.m[0] * v.x + self.m[4] * v.y + self.m[8] * v.z;
        let y = self.m[1] * v.x + self.m[5] * v.y + self.m[9] * v.z;
        let z = self.m[2] * v.x + self.m[6] * v.y + self.m[10] * v.z;
        Vec3::new(x, y, z)
    }

    pub fn look_at(eye: Vec3, target: Vec3, up: Vec3) -> Self {
        let mut f = (target - eye).normalize();
        if f.norm_squared() < 1e-8 {
            f = Vec3::new(0.0, 0.0, -1.0);
        }

        let mut s = f.cross(up.normalize()).normalize();
        if s.norm_squared() < 1e-8 {
            // Pick an arbitrary perpendicular vector if up is collinear with f
            let alt_up = if f.y.abs() < 0.9 {
                Vec3::new(0.0, 1.0, 0.0)
            } else {
                Vec3::new(1.0, 0.0, 0.0)
            };
            s = f.cross(alt_up).normalize();
        }
        let u = s.cross(f);

        let mut m = [0.0; 16];
        // Column 0
        m[0] = s.x;
        m[1] = u.x;
        m[2] = -f.x;
        m[3] = 0.0;

        // Column 1
        m[4] = s.y;
        m[5] = u.y;
        m[6] = -f.y;
        m[7] = 0.0;

        // Column 2
        m[8] = s.z;
        m[9] = u.z;
        m[10] = -f.z;
        m[11] = 0.0;

        // Column 3
        m[12] = -s.dot(eye);
        m[13] = -u.dot(eye);
        m[14] = f.dot(eye);
        m[15] = 1.0;

        Self { m }
    }

    pub fn perspective(fov_y_rad: f32, aspect: f32, z_near: f32, z_far: f32) -> Self {
        let tan_half = (fov_y_rad * 0.5).tan();
        let f = 1.0 / tan_half;

        let mut m = [0.0; 16];
        m[0] = f / aspect;
        m[5] = f;
        m[10] = (z_far + z_near) / (z_near - z_far);
        m[11] = -1.0;
        m[14] = (2.0 * z_far * z_near) / (z_near - z_far);
        m[15] = 0.0;

        Self { m }
    }

    pub fn orthographic(left: f32, right: f32, bottom: f32, top: f32, near: f32, far: f32) -> Self {
        let mut m = [0.0; 16];
        m[0] = 2.0 / (right - left);
        m[5] = 2.0 / (top - bottom);
        m[10] = -2.0 / (far - near);
        m[12] = -(right + left) / (right - left);
        m[13] = -(top + bottom) / (top - bottom);
        m[14] = -(far + near) / (far - near);
        m[15] = 1.0;

        Self { m }
    }

    pub fn transpose(&self) -> Self {
        let mut out = [0.0; 16];
        for col in 0..4 {
            for row in 0..4 {
                out[row * 4 + col] = self.m[col * 4 + row];
            }
        }
        Self { m: out }
    }

    pub fn inverse(&self) -> Option<Self> {
        let m = &self.m;
        let mut inv = [0.0f32; 16];

        inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15]
            + m[9] * m[7] * m[14]
            + m[13] * m[6] * m[11]
            - m[13] * m[7] * m[10];

        inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15]
            - m[8] * m[7] * m[14]
            - m[12] * m[6] * m[11]
            + m[12] * m[7] * m[10];

        inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15]
            + m[8] * m[7] * m[13]
            + m[12] * m[5] * m[11]
            - m[12] * m[7] * m[9];

        inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14]
            - m[8] * m[6] * m[13]
            - m[12] * m[5] * m[10]
            + m[12] * m[6] * m[9];

        inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15]
            - m[9] * m[3] * m[14]
            - m[13] * m[2] * m[11]
            + m[13] * m[3] * m[10];

        inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15]
            + m[8] * m[3] * m[14]
            + m[12] * m[2] * m[11]
            - m[12] * m[3] * m[10];

        inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15]
            - m[8] * m[3] * m[13]
            - m[12] * m[1] * m[11]
            + m[12] * m[3] * m[9];

        inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14]
            + m[8] * m[2] * m[13]
            + m[12] * m[1] * m[10]
            - m[12] * m[2] * m[9];

        inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15]
            + m[5] * m[3] * m[14]
            + m[13] * m[2] * m[7]
            - m[13] * m[3] * m[6];

        inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15]
            - m[4] * m[3] * m[14]
            - m[12] * m[2] * m[7]
            + m[12] * m[3] * m[6];

        inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15]
            + m[4] * m[3] * m[13]
            + m[12] * m[1] * m[7]
            - m[12] * m[3] * m[5];

        inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14]
            - m[4] * m[2] * m[13]
            - m[12] * m[1] * m[6]
            + m[12] * m[2] * m[5];

        inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11]
            - m[5] * m[3] * m[10]
            - m[9] * m[2] * m[7]
            + m[9] * m[3] * m[6];

        inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11]
            + m[4] * m[3] * m[10]
            + m[8] * m[2] * m[7]
            - m[8] * m[3] * m[6];

        inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11]
            - m[4] * m[3] * m[9]
            - m[8] * m[1] * m[7]
            + m[8] * m[3] * m[5];

        inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10]
            + m[4] * m[2] * m[9]
            + m[8] * m[1] * m[6]
            - m[8] * m[2] * m[5];

        let det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
        if det.abs() < 1e-8 {
            return None;
        }

        let inv_det = 1.0 / det;
        for x in &mut inv {
            *x *= inv_det;
        }

        Some(Mat4 { m: inv })
    }
}

impl Mul for Mat4 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Mat4::mul(&self, &rhs)
    }
}

impl Mul<&Mat4> for Mat4 {
    type Output = Self;
    fn mul(self, rhs: &Self) -> Self::Output {
        Mat4::mul(&self, rhs)
    }
}

impl Mul<Mat4> for &Mat4 {
    type Output = Mat4;
    fn mul(self, rhs: Mat4) -> Self::Output {
        Mat4::mul(self, &rhs)
    }
}

impl Mul<&Mat4> for &Mat4 {
    type Output = Mat4;
    fn mul(self, rhs: &Mat4) -> Self::Output {
        Mat4::mul(self, rhs)
    }
}

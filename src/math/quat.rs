use crate::math::mat4::Mat4;
use crate::math::vec3::Vec3;
use std::ops::{Mul, MulAssign};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quat {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32,
}

impl Default for Quat {
    fn default() -> Self {
        Self::identity()
    }
}

impl Quat {
    pub const fn new(x: f32, y: f32, z: f32, w: f32) -> Self {
        Self { x, y, z, w }
    }

    pub const fn identity() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 1.0,
        }
    }

    pub fn from_axis_angle(axis: Vec3, angle_rad: f32) -> Self {
        let half = angle_rad * 0.5;
        let (sin, cos) = half.sin_cos();
        let n = axis.normalize();
        Self {
            x: n.x * sin,
            y: n.y * sin,
            z: n.z * sin,
            w: cos,
        }
    }

    pub fn from_euler(pitch: f32, yaw: f32, roll: f32) -> Self {
        let qx = Self::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), pitch);
        let qy = Self::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), yaw);
        let qz = Self::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), roll);
        qy * qx * qz
    }

    pub fn from_drag(dx: f32, dy: f32, sensitivity: f32) -> Self {
        let q_yaw = Self::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), dx * sensitivity);
        let q_pitch = Self::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), dy * sensitivity);
        q_pitch * q_yaw
    }

    pub fn norm_squared(&self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w
    }

    pub fn norm(&self) -> f32 {
        self.norm_squared().sqrt()
    }

    pub fn normalize(&self) -> Self {
        let n = self.norm();
        if n > 1e-8 {
            Self {
                x: self.x / n,
                y: self.y / n,
                z: self.z / n,
                w: self.w / n,
            }
        } else {
            Self::identity()
        }
    }

    pub fn conjugate(&self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: self.w,
        }
    }

    pub fn inverse(&self) -> Self {
        let n_sq = self.norm_squared();
        if n_sq > 1e-8 {
            let inv_n_sq = 1.0 / n_sq;
            Self {
                x: -self.x * inv_n_sq,
                y: -self.y * inv_n_sq,
                z: -self.z * inv_n_sq,
                w: self.w * inv_n_sq,
            }
        } else {
            Self::identity()
        }
    }

    pub fn mul(&self, other: &Self) -> Self {
        Self {
            x: self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            y: self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
            z: self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w,
            w: self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
        }
    }

    pub fn rotate_vec3(&self, v: Vec3) -> Vec3 {
        let u = Vec3::new(self.x, self.y, self.z);
        let s = self.w;
        let t = u.cross(v) * 2.0;
        v + (t * s) + u.cross(t)
    }

    pub fn to_mat4(&self) -> Mat4 {
        let q = self.normalize();
        let x2 = q.x * 2.0;
        let y2 = q.y * 2.0;
        let z2 = q.z * 2.0;
        let xx = q.x * x2;
        let xy = q.x * y2;
        let xz = q.x * z2;
        let yy = q.y * y2;
        let yz = q.y * z2;
        let zz = q.z * z2;
        let wx = q.w * x2;
        let wy = q.w * y2;
        let wz = q.w * z2;

        let mut m = [0.0; 16];
        // Column 0
        m[0] = 1.0 - (yy + zz);
        m[1] = xy + wz;
        m[2] = xz - wy;
        m[3] = 0.0;

        // Column 1
        m[4] = xy - wz;
        m[5] = 1.0 - (xx + zz);
        m[6] = yz + wx;
        m[7] = 0.0;

        // Column 2
        m[8] = xz + wy;
        m[9] = yz - wx;
        m[10] = 1.0 - (xx + yy);
        m[11] = 0.0;

        // Column 3
        m[12] = 0.0;
        m[13] = 0.0;
        m[14] = 0.0;
        m[15] = 1.0;

        Mat4 { m }
    }

    pub fn slerp(&self, other: &Self, t: f32) -> Self {
        let mut dot = self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w;
        let mut other_w = other.w;
        let mut other_x = other.x;
        let mut other_y = other.y;
        let mut other_z = other.z;

        if dot < 0.0 {
            dot = -dot;
            other_w = -other_w;
            other_x = -other_x;
            other_y = -other_y;
            other_z = -other_z;
        }

        if dot > 0.9995 {
            let q = Self {
                x: self.x + t * (other_x - self.x),
                y: self.y + t * (other_y - self.y),
                z: self.z + t * (other_z - self.z),
                w: self.w + t * (other_w - self.w),
            };
            return q.normalize();
        }

        let theta_0 = dot.clamp(-1.0, 1.0).acos();
        let sin_theta_0 = theta_0.sin();
        let theta = theta_0 * t;
        let sin_theta = theta.sin();

        let s0 = ((1.0 - t) * theta_0).sin() / sin_theta_0;
        let s1 = sin_theta / sin_theta_0;

        Self {
            x: s0 * self.x + s1 * other_x,
            y: s0 * self.y + s1 * other_y,
            z: s0 * self.z + s1 * other_z,
            w: s0 * self.w + s1 * other_w,
        }
    }
}

impl Mul for Quat {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Quat::mul(&self, &rhs)
    }
}

impl Mul<&Quat> for Quat {
    type Output = Self;
    fn mul(self, rhs: &Self) -> Self::Output {
        Quat::mul(&self, rhs)
    }
}

impl Mul<Quat> for &Quat {
    type Output = Quat;
    fn mul(self, rhs: Quat) -> Self::Output {
        Quat::mul(self, &rhs)
    }
}

impl Mul<&Quat> for &Quat {
    type Output = Quat;
    fn mul(self, rhs: &Quat) -> Self::Output {
        Quat::mul(self, rhs)
    }
}

impl MulAssign for Quat {
    fn mul_assign(&mut self, rhs: Self) {
        *self = Quat::mul(self, &rhs);
    }
}


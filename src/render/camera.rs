//! Orbit / arcball camera and 3D projection transforms.
//!
//! Manages camera target, distance, orientation quaternion, perspective projection,
//! and world-to-screen coordinate mapping.

use crate::math::{Mat4, Quat, Vec3};

/// Interactive 3D camera with arcball/turntable orbit, pan, zoom, and projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// Look-at target center in world space
    pub target: Vec3,
    /// Distance from target along view axis
    pub distance: f32,
    /// Rotation orientation quaternion
    pub orientation: Quat,
    /// Vertical field of view in radians
    pub fov: f32,
    /// Aspect ratio (width / height)
    pub aspect: f32,
    /// Near clipping plane distance
    pub near: f32,
    /// Far clipping plane distance
    pub far: f32,
}

impl Default for Camera {
    fn default() -> Self {
        Self::new()
    }
}

impl Camera {
    /// Creates a camera with standard default parameters.
    pub fn new() -> Self {
        Self {
            target: Vec3::ZERO,
            distance: 50.0,
            orientation: Quat::identity(),
            fov: 45.0_f32.to_radians(),
            aspect: 1.0,
            near: 0.1,
            far: 1000.0,
        }
    }

    /// Fits the camera to frame a structure with the given center and bounding sphere radius.
    pub fn fit_structure(&mut self, center: Vec3, radius: f32) {
        self.target = center;
        self.orientation = Quat::identity();

        let fov_half = self.fov * 0.5;
        let sin_half = fov_half.sin();
        let dist = if sin_half > 1e-4 {
            (radius / sin_half).max(radius * 1.5)
        } else {
            radius * 2.5
        };

        self.distance = dist;
        self.near = (dist - radius * 2.5).max(0.1);
        self.far = dist + radius * 2.5;
    }

    /// Rotates the camera around the target using mouse/cursor drag offsets `(dx, dy)`.
    pub fn orbit(&mut self, dx: f32, dy: f32) {
        let delta_q = Quat::from_drag(dx, dy, 0.01);
        self.orientation = (delta_q * self.orientation).normalize();
    }

    /// Pans the camera target along local camera right and up axes.
    pub fn pan(&mut self, dx: f32, dy: f32) {
        let right = self.orientation.rotate_vec3(Vec3::X);
        let up = self.orientation.rotate_vec3(Vec3::Y);
        self.target += right * dx + up * dy;
    }

    /// Zooms the camera by adjusting distance to target.
    pub fn zoom(&mut self, delta: f32) {
        self.distance = (self.distance * (1.0 - delta * 0.1)).clamp(0.5, 10000.0);
    }

    /// Resets the camera orientation, target, and distance to default.
    pub fn reset(&mut self) {
        self.target = Vec3::ZERO;
        self.distance = 50.0;
        self.orientation = Quat::identity();
    }

    /// Returns the camera eye position in world space.
    pub fn eye_position(&self) -> Vec3 {
        self.target
            + self
                .orientation
                .rotate_vec3(Vec3::new(0.0, 0.0, self.distance))
    }

    /// Computes the 4x4 view matrix.
    pub fn view_matrix(&self) -> Mat4 {
        let eye = self.eye_position();
        let up = self.orientation.rotate_vec3(Vec3::Y);
        Mat4::look_at(eye, self.target, up)
    }

    /// Computes the 4x4 perspective projection matrix.
    pub fn proj_matrix(&self) -> Mat4 {
        Mat4::perspective(self.fov, self.aspect, self.near, self.far)
    }

    /// Transforms a 3D world position to screen coordinates `(screen_x, screen_y, view_depth)`.
    ///
    /// Returns `None` if the point is behind the camera or outside the near/far planes.
    pub fn world_to_screen(
        &self,
        world_pos: Vec3,
        width: usize,
        height: usize,
    ) -> Option<(f32, f32, f32)> {
        let view = self.view_matrix();
        let view_pos = view.transform_point(world_pos);
        let view_depth = -view_pos.z;

        if view_depth <= 0.0 || view_depth < self.near || view_depth > self.far {
            return None;
        }

        let proj = self.proj_matrix();
        let vp = proj.mul(&view);

        let m = &vp.m;
        let x = m[0] * world_pos.x + m[4] * world_pos.y + m[8] * world_pos.z + m[12];
        let y = m[1] * world_pos.x + m[5] * world_pos.y + m[9] * world_pos.z + m[13];
        let w = m[3] * world_pos.x + m[7] * world_pos.y + m[11] * world_pos.z + m[15];

        if w <= 0.0 || w.is_nan() {
            return None;
        }

        let inv_w = 1.0 / w;
        let ndc_x = x * inv_w;
        let ndc_y = y * inv_w;

        let screen_x = (ndc_x + 1.0) * 0.5 * (width as f32);
        let screen_y = (1.0 - ndc_y) * 0.5 * (height as f32);

        Some((screen_x, screen_y, view_depth))
    }
}

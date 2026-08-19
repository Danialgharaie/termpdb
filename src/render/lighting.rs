//! Directional lighting, Lambertian diffuse shading, and depth-cue fog.
//!
//! Simulates single-source directional lighting with ambient component and depth cueing
//! to enhance depth perception in 3D terminal rendering.

use crate::math::Vec3;
use crate::render::buffer::PixelColor;

/// Lighting model configuring light direction, ambient/diffuse terms, and depth fog.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Lighting {
    /// Normalized direction pointing towards the light source
    pub light_dir: Vec3,
    /// Ambient light intensity in 0.0..=1.0
    pub ambient_intensity: f32,
    /// Diffuse light intensity in 0.0..=1.0
    pub diffuse_intensity: f32,
    /// Depth cue fog factor in 0.0..=1.0 (0.0 = no fog, 1.0 = full attenuation at max depth)
    pub depth_cue_factor: f32,
}

impl Default for Lighting {
    fn default() -> Self {
        Self {
            light_dir: Vec3::new(0.5, 0.7, 1.0).normalize(),
            ambient_intensity: 0.35,
            diffuse_intensity: 0.65,
            depth_cue_factor: 0.5,
        }
    }
}

impl Lighting {
    /// Creates a new `Lighting` configuration.
    pub fn new(light_dir: Vec3, ambient: f32, diffuse: f32, depth_cue: f32) -> Self {
        Self {
            light_dir: light_dir.normalize(),
            ambient_intensity: ambient.clamp(0.0, 1.0),
            diffuse_intensity: diffuse.clamp(0.0, 1.0),
            depth_cue_factor: depth_cue.clamp(0.0, 1.0),
        }
    }

    /// Computes shaded color using surface normal, depth, and bounding depth range.
    pub fn shade(
        &self,
        normal: Vec3,
        depth: f32,
        base_color: PixelColor,
        min_depth: f32,
        max_depth: f32,
    ) -> PixelColor {
        let n = normal.normalize();
        let l = self.light_dir;
        let n_dot_l = n.dot(l).max(0.0);

        let intensity = (self.ambient_intensity + self.diffuse_intensity * n_dot_l).clamp(0.0, 1.0);

        let fog = if max_depth > min_depth + 1e-4 {
            let t = ((depth - min_depth) / (max_depth - min_depth)).clamp(0.0, 1.0);
            1.0 - self.depth_cue_factor * t
        } else {
            1.0
        };

        let factor = (intensity * fog).clamp(0.0, 1.0);

        let r = (base_color.0 as f32 * factor).round().clamp(0.0, 255.0) as u8;
        let g = (base_color.1 as f32 * factor).round().clamp(0.0, 255.0) as u8;
        let b = (base_color.2 as f32 * factor).round().clamp(0.0, 255.0) as u8;

        (r, g, b)
    }

    /// Computes shaded color without depth range cueing.
    pub fn compute_shade(&self, normal: Vec3, depth: f32, base_color: PixelColor) -> PixelColor {
        self.shade(normal, depth, base_color, 0.0, 0.0)
    }
}

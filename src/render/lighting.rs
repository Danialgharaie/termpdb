//! Directional lighting, Lambertian diffuse shading, Blinn-Phong specular,
//! depth-cue fog, and Depth-of-Field (DoF) focal plane cueing.

use crate::math::Vec3;
use crate::render::buffer::PixelColor;

/// Lighting model configuring light direction, ambient/diffuse/specular terms, depth fog, and DoF.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Lighting {
    /// Normalized direction pointing towards the light source.
    pub light_dir: Vec3,
    /// Ambient light intensity in 0.0..=1.0.
    pub ambient_intensity: f32,
    /// Diffuse (Lambertian) light intensity in 0.0..=1.0.
    pub diffuse_intensity: f32,
    /// Specular (Blinn-Phong) highlight intensity in 0.0..=1.0.
    pub specular_intensity: f32,
    /// Specular shininess exponent (higher = tighter highlight).
    pub shininess: f32,
    /// Depth cue fog factor in 0.0..=1.0 (0.0 = no fog, 1.0 = full attenuation at max depth).
    pub depth_cue_factor: f32,
    /// Optional Depth-of-Field focal plane distance (in view Z coordinates).
    pub dof_focus: Option<f32>,
    /// In-focus depth radius (distance from focal plane before cueing begins).
    pub dof_range: f32,
    /// Precomputed half-vector H = normalize(L + V), where V is the view direction (+Z).
    half_dir: Vec3,
}

impl Default for Lighting {
    fn default() -> Self {
        Self::new(Vec3::new(0.5, 0.7, 1.0), 0.35, 0.65, 0.25, 24.0, 0.5)
    }
}

impl Lighting {
    /// Creates a new Lighting configuration, recomputing the specular half-vector.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        light_dir: Vec3,
        ambient: f32,
        diffuse: f32,
        specular: f32,
        shininess: f32,
        depth_cue: f32,
    ) -> Self {
        let light_dir = light_dir.normalize();
        let half_dir = (light_dir + Vec3::new(0.0, 0.0, 1.0)).normalize();
        Self {
            light_dir,
            ambient_intensity: ambient.clamp(0.0, 1.0),
            diffuse_intensity: diffuse.clamp(0.0, 1.0),
            specular_intensity: specular.clamp(0.0, 1.0),
            shininess: shininess.max(1.0),
            depth_cue_factor: depth_cue.clamp(0.0, 1.0),
            dof_focus: None,
            dof_range: 5.0,
            half_dir,
        }
    }

    /// Computes shaded color using surface normal, depth, bounding depth range, and DoF.
    pub fn shade(
        &self,
        normal: Vec3,
        depth: f32,
        base_color: PixelColor,
        min_depth: f32,
        max_depth: f32,
    ) -> PixelColor {
        let n = normal;
        let l = self.light_dir;
        let n_dot_l = n.dot(l).max(0.0);

        let spec = if self.specular_intensity > 0.0 {
            let n_dot_h = n.dot(self.half_dir).max(0.0);
            self.specular_intensity * n_dot_h.powi(self.shininess as i32)
        } else {
            0.0
        };

        let intensity = (self.ambient_intensity + self.diffuse_intensity * n_dot_l).clamp(0.0, 1.0);

        let fog = if max_depth > min_depth + 1e-4 {
            let t = ((depth - min_depth) / (max_depth - min_depth)).clamp(0.0, 1.0);
            1.0 - self.depth_cue_factor * t
        } else {
            1.0
        };

        let dof_cue = if let Some(focus) = self.dof_focus {
            let dist = (depth - focus).abs();
            if dist <= self.dof_range {
                1.0
            } else {
                let excess = dist - self.dof_range;
                (1.0 - (excess / (self.dof_range * 3.0 + 1e-4)).clamp(0.0, 0.65)).clamp(0.35, 1.0)
            }
        } else {
            1.0
        };

        let factor = (intensity * fog * dof_cue).clamp(0.0, 1.0);
        let spec_factor = (spec * fog * dof_cue).clamp(0.0, 1.0);

        let r = (base_color.0 as f32 * factor + 255.0 * spec_factor)
            .round()
            .clamp(0.0, 255.0) as u8;
        let g = (base_color.1 as f32 * factor + 255.0 * spec_factor)
            .round()
            .clamp(0.0, 255.0) as u8;
        let b = (base_color.2 as f32 * factor + 255.0 * spec_factor)
            .round()
            .clamp(0.0, 255.0) as u8;

        (r, g, b)
    }

    /// Computes shaded color without depth range cueing.
    pub fn compute_shade(&self, normal: Vec3, depth: f32, base_color: PixelColor) -> PixelColor {
        self.shade(normal, depth, base_color, 0.0, 0.0)
    }
}

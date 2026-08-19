//! 3D Molecular structural representations (Trace, Ball & Stick, Ribbon, VDW).
//!
//! Generates 3D geometric visual primitives from molecular data models:
//! - **Trace**: C-alpha / backbone line/cylinder trace with ligand spheres.
//! - **Ball & Stick**: All-atom analytical spheres and covalent bond cylinders.
//! - **Ribbon**: Secondary structure cartoon (helices as thick tubes, sheets as arrow ribbons, coils as smooth tubes).
//! - **VDW**: Van der Waals space-filling analytical spheres.

pub mod ball_stick;
pub mod ribbon;
pub mod trace;
pub mod vdw;

use std::collections::HashMap;

pub use ball_stick::render_ball_stick;
pub use ribbon::render_ribbon;
pub use trace::render_trace;
pub use vdw::render_vdw;

use crate::model::{Residue, Structure};
use crate::render::buffer::{Framebuffer, PixelColor};
use crate::render::camera::Camera;
use crate::render::color::ColorScheme;
use crate::render::lighting::Lighting;

/// Available molecular representation rendering modes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, clap::ValueEnum)]
pub enum RenderMode {
    /// Backbone CA / Nucleic phosphorus trace
    #[value(name = "trace")]
    Trace,
    /// All-atom ball and stick
    #[value(name = "ball-and-stick", alias = "ball_and_stick")]
    BallAndStick,
    /// Secondary structure cartoon ribbon
    #[default]
    #[value(name = "ribbon")]
    Ribbon,
    /// Space-filling Van der Waals spheres
    #[value(name = "vdw")]
    Vdw,
}

impl RenderMode {
    /// Returns the human-readable display name of the render mode.
    pub fn name(&self) -> &'static str {
        match self {
            RenderMode::Trace => "Trace",
            RenderMode::BallAndStick => "Ball & Stick",
            RenderMode::Ribbon => "Ribbon",
            RenderMode::Vdw => "VDW",
        }
    }

    /// Returns an array of all available render modes in cycle order.
    pub fn all() -> &'static [RenderMode] {
        &[
            RenderMode::Trace,
            RenderMode::BallAndStick,
            RenderMode::Ribbon,
            RenderMode::Vdw,
        ]
    }

    /// Cycles to the next render mode.
    pub fn next(&self) -> Self {
        let all = Self::all();
        let idx = all.iter().position(|&m| m == *self).unwrap_or(0);
        all[(idx + 1) % all.len()]
    }

    /// Cycles to the previous render mode.
    pub fn prev(&self) -> Self {
        let all = Self::all();
        let idx = all.iter().position(|&m| m == *self).unwrap_or(0);
        all[(idx + all.len() - 1) % all.len()]
    }
}

/// Computes the projected screen pixel radius from a 3D world radius and view depth.
pub fn project_radius(world_radius: f32, view_depth: f32, fov: f32, height: usize) -> f32 {
    let tan_half = (fov * 0.5).tan();
    if view_depth <= 1e-4 || tan_half <= 1e-4 {
        0.0
    } else {
        (world_radius / (view_depth * tan_half)) * (height as f32 * 0.5)
    }
}

/// Linearly interpolates between two RGB pixel colors.
pub fn lerp_rgb(c1: PixelColor, c2: PixelColor, t: f32) -> PixelColor {
    let t = t.clamp(0.0, 1.0);
    let r = (c1.0 as f32 + t * (c2.0 as f32 - c1.0 as f32)).round() as u8;
    let g = (c1.1 as f32 + t * (c2.1 as f32 - c1.1 as f32)).round() as u8;
    let b = (c1.2 as f32 + t * (c2.2 as f32 - c1.2 as f32)).round() as u8;
    (r, g, b)
}

/// Builds a fast lookup map from atom index to its parent residue.
pub fn build_atom_residue_map(structure: &Structure) -> HashMap<usize, &Residue> {
    let mut map = HashMap::with_capacity(structure.atoms.len());
    for chain in &structure.chains {
        for res in &chain.residues {
            for &atom_idx in &res.atom_indices {
                map.insert(atom_idx, res);
            }
        }
    }
    map
}

/// Orchestrates 3D rendering of a molecular structure into the given framebuffer.
pub fn render_structure(
    structure: &Structure,
    mode: RenderMode,
    color_scheme: ColorScheme,
    camera: &Camera,
    buffer: &mut Framebuffer,
    lighting: &Lighting,
) {
    if structure.atoms.is_empty() {
        return;
    }

    match mode {
        RenderMode::Trace => trace::render_trace(structure, color_scheme, camera, buffer, lighting),
        RenderMode::BallAndStick => {
            ball_stick::render_ball_stick(structure, color_scheme, camera, buffer, lighting)
        }
        RenderMode::Ribbon => {
            ribbon::render_ribbon(structure, color_scheme, camera, buffer, lighting)
        }
        RenderMode::Vdw => vdw::render_vdw(structure, color_scheme, camera, buffer, lighting),
    }
}

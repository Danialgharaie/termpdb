//! Color schemes and per-atom color assignment.
//!
//! Supports standard CPK element colors, rainbow sequence spectrum, per-chain palettes,
//! secondary structure assignment, B-factor / pLDDT confidence heatmaps, and hydrophobicity scales.

use crate::model::{Atom, Residue, SecondaryStructure, Structure};
use crate::render::buffer::PixelColor;

/// Available color schemes for molecular visualization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ColorScheme {
    /// Element-based CPK coloring (Carbon grey, Oxygen red, Nitrogen blue, etc.)
    #[default]
    Cpk,
    /// N-to-C terminal rainbow gradient (Blue -> Green -> Yellow -> Red)
    Rainbow,
    /// Unique color per polymer chain
    Chain,
    /// Color by secondary structure (Helix magenta, Sheet yellow, Coil cyan/grey)
    SecondaryStructure,
    /// Temperature factor (B-factor / pLDDT) heatmap
    BFactor,
    /// Kyte-Doolittle hydrophobicity scale (Hydrophobic orange, Hydrophilic blue)
    Hydrophobicity,
}

impl ColorScheme {
    /// Returns the human-readable display name of the color scheme.
    pub fn name(&self) -> &'static str {
        match self {
            ColorScheme::Cpk => "CPK",
            ColorScheme::Rainbow => "Rainbow",
            ColorScheme::Chain => "Chain",
            ColorScheme::SecondaryStructure => "Secondary Structure",
            ColorScheme::BFactor => "B-Factor",
            ColorScheme::Hydrophobicity => "Hydrophobicity",
        }
    }

    /// Returns an array of all available color schemes in cycle order.
    pub fn all() -> &'static [ColorScheme] {
        &[
            ColorScheme::Cpk,
            ColorScheme::Rainbow,
            ColorScheme::Chain,
            ColorScheme::SecondaryStructure,
            ColorScheme::BFactor,
            ColorScheme::Hydrophobicity,
        ]
    }

    /// Cycles to the next color scheme.
    pub fn next(&self) -> Self {
        let all = Self::all();
        let idx = all.iter().position(|&s| s == *self).unwrap_or(0);
        all[(idx + 1) % all.len()]
    }

    /// Cycles to the previous color scheme.
    pub fn prev(&self) -> Self {
        let all = Self::all();
        let idx = all.iter().position(|&s| s == *self).unwrap_or(0);
        all[(idx + all.len() - 1) % all.len()]
    }
}

/// Palette of distinct colors for polymer chains.
const CHAIN_COLORS: &[PixelColor] = &[
    (0, 168, 204),   // Teal / Cyan
    (243, 129, 129), // Coral / Salmon
    (255, 184, 76),  // Amber / Gold
    (155, 89, 182),  // Purple
    (46, 204, 113),  // Emerald Green
    (231, 76, 60),   // Crimson Red
    (52, 152, 219),  // Sky Blue
    (241, 196, 15),  // Bright Yellow
    (230, 126, 34),  // Orange
    (26, 188, 156),  // Turquoise
    (214, 162, 232), // Lavender
    (149, 175, 192), // Slate Blue
];

/// Converts HSV color to RGB. `h_deg` in [0, 360), `s` in [0, 1], `v` in [0, 1].
fn hsv_to_rgb(h_deg: f32, s: f32, v: f32) -> PixelColor {
    let h = (h_deg.rem_euclid(360.0)) / 60.0;
    let c = v * s;
    let x = c * (1.0 - ((h % 2.0) - 1.0).abs());
    let m = v - c;

    let (r1, g1, b1) = match h as u32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };

    (
        ((r1 + m) * 255.0).round().clamp(0.0, 255.0) as u8,
        ((g1 + m) * 255.0).round().clamp(0.0, 255.0) as u8,
        ((b1 + m) * 255.0).round().clamp(0.0, 255.0) as u8,
    )
}

/// Linearly interpolates between two RGB colors.
fn lerp_rgb(c1: PixelColor, c2: PixelColor, t: f32) -> PixelColor {
    let t = t.clamp(0.0, 1.0);
    let r = (c1.0 as f32 + t * (c2.0 as f32 - c1.0 as f32)).round() as u8;
    let g = (c1.1 as f32 + t * (c2.1 as f32 - c1.1 as f32)).round() as u8;
    let b = (c1.2 as f32 + t * (c2.2 as f32 - c1.2 as f32)).round() as u8;
    (r, g, b)
}

/// Computes the RGB color for an atom given the active color scheme.
pub fn color_for_atom(
    atom: &Atom,
    residue: Option<&Residue>,
    structure: &Structure,
    scheme: ColorScheme,
) -> PixelColor {
    match scheme {
        ColorScheme::Cpk => atom.element.cpk_color,

        ColorScheme::Rainbow => {
            let total = structure.atom_count();
            let t = if total > 1 {
                (atom.index as f32) / ((total - 1) as f32)
            } else {
                0.0
            };
            // Blue (240 deg) to Red (0 deg)
            let hue = 240.0 * (1.0 - t.clamp(0.0, 1.0));
            hsv_to_rgb(hue, 1.0, 1.0)
        }

        ColorScheme::Chain => {
            let chain_idx = structure
                .chains
                .iter()
                .position(|c| c.id == atom.chain_id)
                .unwrap_or_else(|| {
                    // Fallback to simple hash of chain_id
                    let mut h: usize = 0;
                    for b in atom.chain_id.bytes() {
                        h = h.wrapping_mul(31).wrapping_add(b as usize);
                    }
                    h
                });
            CHAIN_COLORS[chain_idx % CHAIN_COLORS.len()]
        }

        ColorScheme::SecondaryStructure => {
            let ss = residue
                .map(|r| r.secondary_structure)
                .unwrap_or(SecondaryStructure::Coil);

            match ss {
                SecondaryStructure::Helix => (235, 60, 150), // Magenta / Purple
                SecondaryStructure::Sheet => (255, 200, 0),  // Yellow / Gold
                SecondaryStructure::Coil => (120, 170, 220), // Cyan / Soft Blue
            }
        }

        ColorScheme::BFactor => {
            let (min_b, max_b) = structure.b_factor_range();
            let t = if max_b > min_b + 1e-4 {
                ((atom.b_factor - min_b) / (max_b - min_b)).clamp(0.0, 1.0)
            } else {
                0.5
            };
            // Blue (240 deg) -> Cyan -> Green -> Yellow -> Red (0 deg)
            let hue = 240.0 * (1.0 - t);
            hsv_to_rgb(hue, 1.0, 1.0)
        }

        ColorScheme::Hydrophobicity => {
            let score = residue.map(|r| r.hydrophobicity_score()).unwrap_or(0.0);
            // Range from -4.5 (hydrophilic) to +4.5 (hydrophobic)
            let norm = ((score + 4.5) / 9.0).clamp(0.0, 1.0);

            let blue = (40, 100, 240);
            let neutral = (220, 220, 220);
            let orange_red = (230, 60, 20);

            if norm < 0.5 {
                lerp_rgb(blue, neutral, norm * 2.0)
            } else {
                lerp_rgb(neutral, orange_red, (norm - 0.5) * 2.0)
            }
        }
    }
}

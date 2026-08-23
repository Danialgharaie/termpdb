//! Color schemes and per-atom color assignment.
//!
//! Supports standard CPK element colors, rainbow sequence spectrum, per-chain palettes,
//! secondary structure assignment, B-factor / pLDDT confidence heatmaps, hydrophobicity scales,
//! electrostatic potential, and curated terminal themes (Catppuccin, Nord, Tokyo Night, Gruvbox).

use std::collections::HashMap;

use crate::model::{Atom, Residue, SecondaryStructure, Structure};
use crate::render::buffer::PixelColor;

/// Available color schemes for molecular visualization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, clap::ValueEnum)]
pub enum ColorScheme {
    /// Element-based CPK coloring (Carbon grey, Oxygen red, Nitrogen blue, etc.)
    #[default]
    #[value(name = "cpk")]
    Cpk,
    /// N-to-C terminal rainbow gradient (Blue -> Green -> Yellow -> Red)
    #[value(name = "rainbow")]
    Rainbow,
    /// Unique color per polymer chain
    #[value(name = "chain")]
    Chain,
    /// Color by secondary structure (Helix magenta, Sheet yellow, Coil cyan/grey)
    #[value(name = "secondary-structure", alias = "ss")]
    SecondaryStructure,
    /// Temperature factor (B-factor) heatmap
    #[value(name = "bfactor", alias = "b-factor")]
    BFactor,
    /// AlphaFold pLDDT confidence scale (Very high blue, confident cyan, low yellow, very low orange)
    #[value(name = "plddt")]
    Plddt,
    /// Kyte-Doolittle hydrophobicity scale (Hydrophobic orange, Hydrophilic blue)
    #[value(name = "hydrophobicity")]
    Hydrophobicity,
    /// Electrostatic potential (Acidic/Negative red, Basic/Positive blue, Neutral white)
    #[value(name = "electrostatic", alias = "charge")]
    Electrostatic,
    /// Catppuccin Mocha pastel palette
    #[value(name = "catppuccin")]
    Catppuccin,
    /// Arctic Nord palette
    #[value(name = "nord")]
    Nord,
    /// Tokyo Night neon palette
    #[value(name = "tokyo-night", alias = "tokyonight")]
    TokyoNight,
    /// Warm retro Gruvbox palette
    #[value(name = "gruvbox")]
    Gruvbox,
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
            ColorScheme::Plddt => "AlphaFold pLDDT",
            ColorScheme::Hydrophobicity => "Hydrophobicity",
            ColorScheme::Electrostatic => "Electrostatic",
            ColorScheme::Catppuccin => "Catppuccin",
            ColorScheme::Nord => "Nord",
            ColorScheme::TokyoNight => "Tokyo Night",
            ColorScheme::Gruvbox => "Gruvbox",
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
            ColorScheme::Plddt,
            ColorScheme::Hydrophobicity,
            ColorScheme::Electrostatic,
            ColorScheme::Catppuccin,
            ColorScheme::Nord,
            ColorScheme::TokyoNight,
            ColorScheme::Gruvbox,
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

const CATPPUCCIN_COLORS: &[PixelColor] = &[
    (203, 166, 247), // Mauve
    (245, 194, 231), // Pink
    (242, 205, 205), // Flamingo
    (243, 139, 168), // Red
    (250, 179, 135), // Peach
    (249, 226, 175), // Yellow
    (166, 227, 161), // Green
    (148, 226, 213), // Teal
    (137, 220, 235), // Sky
    (116, 199, 236), // Sapphire
    (137, 180, 250), // Blue
    (180, 190, 254), // Lavender
];

const NORD_COLORS: &[PixelColor] = &[
    (136, 192, 208), // Frost Blue
    (143, 188, 187), // Frost Teal
    (129, 161, 193), // Frost Slate
    (94, 129, 172),  // Frost Deep
    (191, 97, 106),  // Aurora Red
    (208, 135, 112), // Aurora Orange
    (235, 203, 139), // Aurora Yellow
    (163, 190, 140), // Aurora Green
    (180, 142, 173), // Aurora Purple
];

const TOKYO_NIGHT_COLORS: &[PixelColor] = &[
    (125, 207, 255), // Neon Cyan
    (122, 162, 247), // Blue
    (187, 154, 247), // Purple
    (247, 118, 142), // Red/Pink
    (255, 158, 100), // Orange
    (224, 175, 104), // Yellow
    (158, 206, 106), // Green
    (65, 166, 181),  // Teal
];

const GRUVBOX_COLORS: &[PixelColor] = &[
    (251, 73, 52),   // Red
    (184, 187, 38),  // Green
    (250, 189, 47),  // Yellow
    (131, 165, 152), // Blue
    (211, 134, 155), // Purple
    (142, 192, 124), // Aqua
    (254, 128, 25),  // Orange
    (235, 219, 178), // Light Gray
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

/// Structure-wide statistics needed by some color schemes.
#[derive(Debug, Clone)]
pub struct ColorStats {
    /// Total atom count (Rainbow scheme).
    pub atom_count: usize,
    /// Min/max B-factor across the active view (BFactor scheme).
    pub b_factor_range: (f32, f32),
    /// chain id -> chain position, for the Chain scheme.
    pub chain_index: HashMap<String, usize>,
}

impl ColorStats {
    /// Computes the stats for the active view in a single O(n) pass.
    pub fn for_structure(structure: &Structure) -> Self {
        Self {
            atom_count: structure.atom_count(),
            b_factor_range: structure.b_factor_range(),
            chain_index: structure
                .chains()
                .iter()
                .enumerate()
                .map(|(i, c)| (c.id.clone(), i))
                .collect(),
        }
    }
}

/// Computes the RGB color for an atom using precomputed ColorStats.
pub fn color_for_atom_with_stats(
    atom: &Atom,
    residue: Option<&Residue>,
    scheme: ColorScheme,
    stats: &ColorStats,
) -> PixelColor {
    let get_chain_idx = || {
        stats
            .chain_index
            .get(&atom.chain_id)
            .copied()
            .unwrap_or_else(|| {
                let mut h: usize = 0;
                for b in atom.chain_id.bytes() {
                    h = h.wrapping_mul(31).wrapping_add(b as usize);
                }
                h
            })
    };

    match scheme {
        ColorScheme::Cpk => atom.element.cpk_color,

        ColorScheme::Rainbow => {
            let total = stats.atom_count;
            let t = if total > 1 {
                (atom.index as f32) / ((total - 1) as f32)
            } else {
                0.0
            };
            let hue = 240.0 * (1.0 - t.clamp(0.0, 1.0));
            hsv_to_rgb(hue, 1.0, 1.0)
        }

        ColorScheme::Chain => {
            let idx = get_chain_idx();
            CHAIN_COLORS[idx % CHAIN_COLORS.len()]
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
            let (min_b, max_b) = stats.b_factor_range;
            let t = if max_b > min_b + 1e-4 {
                ((atom.b_factor - min_b) / (max_b - min_b)).clamp(0.0, 1.0)
            } else {
                0.5
            };
            let hue = 240.0 * (1.0 - t);
            hsv_to_rgb(hue, 1.0, 1.0)
        }

        ColorScheme::Plddt => {
            // AlphaFold official pLDDT thresholds:
            // > 90: Very high (Dark Blue)
            // 70 - 90: Confident (Cyan)
            // 50 - 70: Low (Yellow)
            // < 50: Very low (Orange)
            let plddt = atom.b_factor;
            if plddt >= 90.0 {
                (0, 83, 214)
            } else if plddt >= 70.0 {
                (101, 203, 243)
            } else if plddt >= 50.0 {
                (255, 219, 19)
            } else {
                (255, 125, 69)
            }
        }

        ColorScheme::Hydrophobicity => {
            let score = residue.map(|r| r.hydrophobicity_score()).unwrap_or(0.0);
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

        ColorScheme::Electrostatic => {
            let res_name = atom.res_name.trim().to_ascii_uppercase();
            match res_name.as_str() {
                // Negative (Acidic): Red
                "ASP" | "GLU" => (235, 45, 45),
                // Positive (Basic): Blue
                "LYS" | "ARG" | "HIS" => (45, 95, 235),
                // Polar: Cyan / Green
                "SER" | "THR" | "ASN" | "GLN" | "CYS" | "TYR" => (50, 190, 160),
                // Hydrophobic / Nonpolar: White / Neutral
                _ => (220, 220, 220),
            }
        }

        ColorScheme::Catppuccin => {
            let idx = get_chain_idx();
            CATPPUCCIN_COLORS[idx % CATPPUCCIN_COLORS.len()]
        }

        ColorScheme::Nord => {
            let idx = get_chain_idx();
            NORD_COLORS[idx % NORD_COLORS.len()]
        }

        ColorScheme::TokyoNight => {
            let idx = get_chain_idx();
            TOKYO_NIGHT_COLORS[idx % TOKYO_NIGHT_COLORS.len()]
        }

        ColorScheme::Gruvbox => {
            let idx = get_chain_idx();
            GRUVBOX_COLORS[idx % GRUVBOX_COLORS.len()]
        }
    }
}

/// Computes the RGB color for an atom given the active color scheme.
pub fn color_for_atom(
    atom: &Atom,
    residue: Option<&Residue>,
    structure: &Structure,
    scheme: ColorScheme,
) -> PixelColor {
    let stats = ColorStats::for_structure(structure);
    color_for_atom_with_stats(atom, residue, scheme, &stats)
}

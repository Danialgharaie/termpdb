//! Nucleic acid (DNA/RNA) double-helix cartoon ladder representation.
//!
//! Generates smooth sugar-phosphate backbone ribbons and rectangular planar base-pair slabs:
//! - Adenine (A/DA): Forest Green `(50, 180, 70)`
//! - Thymine (T/DT): Crimson Red `(220, 50, 50)`
//! - Guanine (G/DG): Golden Yellow `(230, 190, 30)`
//! - Cytosine (C/DC): Bright Cyan `(40, 190, 220)`
//! - Uracil (U): Coral / Rose `(220, 80, 120)`

use crate::math::Vec3;
use crate::model::atom::Atom;
use crate::model::residue::Residue;
use crate::render::buffer::PixelColor;
use crate::render::representations::ribbon::RibbonPrimitive;

/// Checks if a residue name is a standard or modified DNA/RNA nucleotide.
pub fn is_nucleic_residue(res_name: &str) -> bool {
    let s = res_name.trim().to_ascii_uppercase();
    matches!(
        s.as_str(),
        "A" | "C"
            | "G"
            | "T"
            | "U"
            | "DA"
            | "DC"
            | "DG"
            | "DT"
            | "DI"
            | "I"
            | "5MC"
            | "OMC"
            | "1MA"
            | "2MG"
            | "7MG"
            | "H2U"
            | "PSU"
    )
}

/// Returns the publication standard color for a nucleic acid base.
pub fn base_color(res_name: &str) -> PixelColor {
    let s = res_name.trim().to_ascii_uppercase();
    match s.as_str() {
        "A" | "DA" | "1MA" => (50, 190, 80),          // Green
        "T" | "DT" => (225, 50, 50),                  // Red
        "G" | "DG" | "2MG" | "7MG" => (235, 195, 35), // Yellow
        "C" | "DC" | "5MC" | "OMC" => (40, 190, 230), // Cyan
        "U" | "H2U" | "PSU" => (230, 85, 130),        // Rose
        _ => (160, 160, 160),                         // Gray
    }
}

/// Generates a rectangular base slab primitive for a nucleic acid residue.
pub fn build_base_slab(
    residue: &Residue,
    atoms: &[Atom],
    sugar_pos: Vec3,
) -> Option<RibbonPrimitive> {
    // Find base nitrogen / ring atoms (N1, N9, N3, C2, C4, C6, etc.)
    let mut base_centroid = Vec3::ZERO;
    let mut count = 0;

    for &idx in &residue.atom_indices {
        let Some(atom) = atoms.get(idx) else { continue };
        let name = atom.name.trim();
        // Ignore phosphate and ribose atoms (P, OP1, OP2, O5', C5', C4', O4', C3', O3', C2', O2', C1')
        if name.starts_with('P')
            || name.starts_with("OP")
            || name.ends_with('\'')
            || name.ends_with('*')
        {
            continue;
        }
        base_centroid += atom.pos;
        count += 1;
    }

    if count == 0 {
        return None;
    }

    base_centroid /= count as f32;
    let color = base_color(&residue.name);

    // Create a thick cylinder or box from the sugar position to the base centroid
    Some(RibbonPrimitive::Cylinder {
        a: sugar_pos,
        b: base_centroid,
        r_world: 0.35,
        min_r: 1.0,
        color,
    })
}

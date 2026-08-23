//! Non-covalent interaction and disulfide bridge detection.

use crate::model::Structure;
use crate::render::buffer::PixelColor;

/// Types of molecular interactions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InteractionKind {
    /// Covalent disulfide bridge between two Cysteine sulfur (SG) atoms.
    Disulfide,
    /// Polar donor-acceptor hydrogen bond (N/O ... N/O).
    HydrogenBond,
    /// Ionic salt bridge between basic (Lys/Arg/His) and acidic (Asp/Glu) residues.
    SaltBridge,
}

impl InteractionKind {
    /// Returns default rendering color for this interaction kind.
    pub fn default_color(&self) -> PixelColor {
        match self {
            InteractionKind::Disulfide => (255, 215, 0),      // Gold
            InteractionKind::HydrogenBond => (100, 240, 255), // Cyan
            InteractionKind::SaltBridge => (255, 80, 200),    // Magenta
        }
    }
}

/// A detected molecular interaction connecting two atoms.
#[derive(Debug, Clone, PartialEq)]
pub struct Interaction {
    pub atom1_idx: usize,
    pub atom2_idx: usize,
    pub kind: InteractionKind,
    pub distance: f32,
}

/// Detects disulfide bridges, hydrogen bonds, and salt bridges in a structure.
#[allow(clippy::needless_range_loop)]
pub fn detect_interactions(structure: &Structure) -> Vec<Interaction> {
    let atoms = structure.atoms();
    let mut interactions = Vec::new();

    for i in 0..atoms.len() {
        let a1 = &atoms[i];
        let sym1 = a1.element.symbol;
        let is_sg1 = a1.name.trim() == "SG" && a1.res_name.trim().eq_ignore_ascii_case("CYS");
        let is_polar1 = sym1 == "N" || sym1 == "O";

        for j in (i + 1)..atoms.len() {
            let a2 = &atoms[j];
            let sym2 = a2.element.symbol;
            let dist = a1.pos.distance(&a2.pos);

            // 1. Disulfide bonds: CYS SG to CYS SG within 1.8..=2.4 A
            if is_sg1
                && a2.name.trim() == "SG"
                && a2.res_name.trim().eq_ignore_ascii_case("CYS")
                && (1.8..=2.4).contains(&dist)
            {
                interactions.push(Interaction {
                    atom1_idx: i,
                    atom2_idx: j,
                    kind: InteractionKind::Disulfide,
                    distance: dist,
                });
                continue;
            }

            // 2. Hydrogen bond / Polar interaction
            if is_polar1 && (sym2 == "N" || sym2 == "O") {
                if a1.chain_id == a2.chain_id && a1.res_seq == a2.res_seq {
                    continue;
                }
                if a1.chain_id == a2.chain_id && (a1.res_seq - a2.res_seq).abs() == 1 {
                    let n1 = a1.name.trim();
                    let n2 = a2.name.trim();
                    if (n1 == "C" && n2 == "N")
                        || (n1 == "N" && n2 == "C")
                        || (n1 == "O" && n2 == "N")
                        || (n1 == "N" && n2 == "O")
                    {
                        continue;
                    }
                }

                if (2.4..=3.5).contains(&dist) {
                    interactions.push(Interaction {
                        atom1_idx: i,
                        atom2_idx: j,
                        kind: InteractionKind::HydrogenBond,
                        distance: dist,
                    });
                }
            }
        }
    }

    interactions
}

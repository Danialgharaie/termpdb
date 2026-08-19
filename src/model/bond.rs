//! Covalent bond graph and spatial hash grid bond detection.

use std::collections::HashMap;

use crate::model::atom::Atom;

/// Bond order representation (single, double, triple, aromatic).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    Other(f32),
}

impl BondOrder {
    /// Returns the numerical bond order value.
    pub fn as_f32(&self) -> f32 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
            BondOrder::Other(o) => *o,
        }
    }

    /// Converts from numerical bond order value.
    pub fn from_f32(order: f32) -> Self {
        if (order - 1.0).abs() < 1e-3 {
            BondOrder::Single
        } else if (order - 2.0).abs() < 1e-3 {
            BondOrder::Double
        } else if (order - 3.0).abs() < 1e-3 {
            BondOrder::Triple
        } else if (order - 1.5).abs() < 1e-3 {
            BondOrder::Aromatic
        } else {
            BondOrder::Other(order)
        }
    }
}

/// Represents a chemical bond between two atoms in the structure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bond {
    /// Index of the first atom in the `Structure.atoms` array
    pub atom1_idx: usize,
    /// Index of the second atom in the `Structure.atoms` array
    pub atom2_idx: usize,
    /// Numerical bond order (1.0 = single, 2.0 = double, 1.5 = aromatic)
    pub order: f32,
}

impl Bond {
    /// Creates a new bond between two atom indices with a given bond order.
    pub fn new(atom1_idx: usize, atom2_idx: usize, order: f32) -> Self {
        Self {
            atom1_idx,
            atom2_idx,
            order,
        }
    }

    /// Creates a standard single bond (order 1.0).
    pub fn single(atom1_idx: usize, atom2_idx: usize) -> Self {
        Self::new(atom1_idx, atom2_idx, 1.0)
    }

    /// Returns the other atom index given one end of the bond.
    pub fn other(&self, atom_idx: usize) -> Option<usize> {
        if self.atom1_idx == atom_idx {
            Some(self.atom2_idx)
        } else if self.atom2_idx == atom_idx {
            Some(self.atom1_idx)
        } else {
            None
        }
    }
}

/// Fast O(N) spatial hash grid for detecting covalent bonds between atoms.
pub struct BondDetector;

impl BondDetector {
    /// Default covalent bond distance tolerance in Å (added to sum of covalent radii).
    pub const DEFAULT_TOLERANCE: f32 = 0.45;
    /// Default grid cell size in Å.
    pub const DEFAULT_CELL_SIZE: f32 = 3.5;
    /// Minimum valid bond distance in Å to filter out coincident/overlapping atoms.
    pub const MIN_BOND_DISTANCE: f32 = 0.40;

    /// Automatically detects all covalent bonds for the given slice of atoms
    /// using standard covalent radii and spatial hashing.
    pub fn detect_bonds(atoms: &[Atom]) -> Vec<Bond> {
        Self::detect_bonds_with_cutoff(atoms, Self::DEFAULT_TOLERANCE, Self::DEFAULT_CELL_SIZE)
    }

    /// Detects covalent bonds with custom tolerance and spatial cell size.
    pub fn detect_bonds_with_cutoff(atoms: &[Atom], tolerance: f32, cell_size: f32) -> Vec<Bond> {
        if atoms.len() < 2 {
            return Vec::new();
        }

        let inv_cell = 1.0 / cell_size;
        let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();

        // 1. Bin all atoms into 3D grid cells
        for (idx, atom) in atoms.iter().enumerate() {
            let gx = (atom.pos.x * inv_cell).floor() as i32;
            let gy = (atom.pos.y * inv_cell).floor() as i32;
            let gz = (atom.pos.z * inv_cell).floor() as i32;
            grid.entry((gx, gy, gz)).or_default().push(idx);
        }

        let mut bonds = Vec::new();

        // 2. Query each atom and its 27 neighbor cells
        for (idx1, atom1) in atoms.iter().enumerate() {
            let gx = (atom1.pos.x * inv_cell).floor() as i32;
            let gy = (atom1.pos.y * inv_cell).floor() as i32;
            let gz = (atom1.pos.z * inv_cell).floor() as i32;

            let r1 = atom1.covalent_radius();

            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let neighbor_cell = (gx + dx, gy + dy, gz + dz);
                        if let Some(cell_atoms) = grid.get(&neighbor_cell) {
                            for &idx2 in cell_atoms {
                                // Only check pairs once and avoid self-pairing
                                if idx2 <= idx1 {
                                    continue;
                                }

                                let atom2 = &atoms[idx2];
                                let dist = atom1.pos.distance(&atom2.pos);

                                if dist < Self::MIN_BOND_DISTANCE {
                                    continue;
                                }

                                let max_allowed = r1 + atom2.covalent_radius() + tolerance;
                                if dist <= max_allowed {
                                    bonds.push(Bond::single(idx1, idx2));
                                }
                            }
                        }
                    }
                }
            }
        }

        bonds
    }
}

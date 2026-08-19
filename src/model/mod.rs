//! Chemical and molecular data models.
//!
//! Provides the core domain representation for atoms, residues, chains, bonds,
//! elements, and complete macromolecular structures.

pub mod atom;
pub mod bond;
pub mod chain;
pub mod elements;
pub mod residue;

use std::collections::HashMap;

pub use atom::Atom;
pub use bond::{Bond, BondDetector, BondOrder};
pub use chain::Chain;
pub use elements::{ELEMENTS, Element, element_by_atomic_number, element_by_symbol};
pub use residue::{Residue, SecondaryStructure};

use crate::math::Vec3;

/// Complete molecular structure containing chains, atoms, bonds, and metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct Structure {
    /// Structure title or description (e.g. from PDB TITLE record)
    pub title: String,
    /// 4-letter PDB ID or identifier if known (e.g. "1CRN")
    pub id_code: Option<String>,
    /// Polymer chains in this structure
    pub chains: Vec<Chain>,
    /// Flat array of all atoms across all chains and heteroatoms
    pub atoms: Vec<Atom>,
    /// Covalent bonds connecting atoms
    pub bonds: Vec<Bond>,
    /// Arbitrary key-value metadata (e.g. resolution, experimental method, deposition date)
    pub metadata: HashMap<String, String>,
}

impl Structure {
    /// Creates a new empty `Structure` with the given title.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            id_code: None,
            chains: Vec::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            metadata: HashMap::new(),
        }
    }

    /// Creates a new empty `Structure` with ID code and title.
    pub fn with_id(id_code: impl Into<String>, title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            id_code: Some(id_code.into()),
            chains: Vec::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            metadata: HashMap::new(),
        }
    }

    /// Adds an atom to the structure, assigning its sequential index.
    /// Returns the index of the newly added atom.
    pub fn add_atom(&mut self, mut atom: Atom) -> usize {
        let idx = self.atoms.len();
        atom.index = idx;
        self.atoms.push(atom);
        idx
    }

    /// Adds a chain to the structure.
    pub fn add_chain(&mut self, chain: Chain) {
        self.chains.push(chain);
    }

    /// Adds a bond to the structure.
    pub fn add_bond(&mut self, bond: Bond) {
        self.bonds.push(bond);
    }

    /// Total number of atoms in the structure.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Total number of chains in the structure.
    pub fn chain_count(&self) -> usize {
        self.chains.len()
    }

    /// Total number of residues across all chains.
    pub fn residue_count(&self) -> usize {
        self.chains.iter().map(|c| c.residue_count()).sum()
    }

    /// Finds a chain by its ID string.
    pub fn get_chain(&self, id: &str) -> Option<&Chain> {
        self.chains.iter().find(|c| c.id == id)
    }

    /// Finds a mutable chain by its ID string.
    pub fn get_chain_mut(&mut self, id: &str) -> Option<&mut Chain> {
        self.chains.iter_mut().find(|c| c.id == id)
    }

    /// Returns references to all C-alpha (CA) atoms across all chains in sequence order.
    pub fn ca_atoms(&self) -> Vec<&Atom> {
        self.chains
            .iter()
            .flat_map(|c| c.ca_atoms(&self.atoms))
            .collect()
    }

    /// Computes the geometric center of mass (centroid) of all atoms in the structure.
    pub fn center_of_mass(&self) -> Vec3 {
        if self.atoms.is_empty() {
            return Vec3::ZERO;
        }

        let mut sum = Vec3::ZERO;
        for atom in &self.atoms {
            sum += atom.pos;
        }

        sum / (self.atoms.len() as f32)
    }

    /// Computes the radius of the minimum bounding sphere centered at the center of mass.
    pub fn bounding_sphere_radius(&self) -> f32 {
        if self.atoms.is_empty() {
            return 1.0;
        }

        let com = self.center_of_mass();
        let mut max_dist_sq: f32 = 0.0;

        for atom in &self.atoms {
            let dist_sq = atom.pos.distance_squared(&com);
            if dist_sq > max_dist_sq {
                max_dist_sq = dist_sq;
            }
        }

        let r = max_dist_sq.sqrt();
        if r < 1e-4 { 1.0 } else { r }
    }

    /// Translates all atom positions so that the center of mass is centered at `(0, 0, 0)`.
    pub fn center_and_normalize(&mut self) {
        if self.atoms.is_empty() {
            return;
        }

        let com = self.center_of_mass();
        for atom in &mut self.atoms {
            atom.pos -= com;
        }
    }

    /// Automatically detects and builds covalent bonds for all atoms using the spatial hash grid.
    pub fn build_bonds(&mut self) {
        self.bonds = BondDetector::detect_bonds(&self.atoms);
    }

    /// Returns the minimum and maximum B-factor in the structure (or `(0.0, 100.0)` if empty).
    pub fn b_factor_range(&self) -> (f32, f32) {
        if self.atoms.is_empty() {
            return (0.0, 100.0);
        }

        let mut min_b = f32::INFINITY;
        let mut max_b = f32::NEG_INFINITY;

        for atom in &self.atoms {
            if atom.b_factor < min_b {
                min_b = atom.b_factor;
            }
            if atom.b_factor > max_b {
                max_b = atom.b_factor;
            }
        }

        if min_b.is_infinite() || max_b.is_infinite() {
            (0.0, 100.0)
        } else {
            (min_b, max_b)
        }
    }
}

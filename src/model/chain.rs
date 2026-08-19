//! Molecular chain representation.

use crate::model::atom::Atom;
use crate::model::residue::Residue;

/// Represents a single macromolecular chain (e.g. Chain "A") containing residues.
#[derive(Debug, Clone, PartialEq)]
pub struct Chain {
    /// Chain identifier (e.g. "A", "B", "1")
    pub id: String,
    /// Ordered list of residues in this chain
    pub residues: Vec<Residue>,
}

impl Chain {
    /// Creates a new empty chain with the given ID.
    pub fn new(id: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            residues: Vec::new(),
        }
    }

    /// Returns references to all C-alpha (CA) atoms in sequential order along the chain.
    pub fn ca_atoms<'a>(&self, atoms: &'a [Atom]) -> Vec<&'a Atom> {
        self.residues
            .iter()
            .filter_map(|r| r.ca_atom(atoms))
            .collect()
    }

    /// Returns indices of all C-alpha (CA) atoms in sequential order along the chain.
    pub fn ca_atom_indices(&self, atoms: &[Atom]) -> Vec<usize> {
        self.residues
            .iter()
            .filter_map(|r| r.ca_atom_index(atoms))
            .collect()
    }

    /// Returns total number of atoms in this chain across all residues.
    pub fn atom_count(&self) -> usize {
        self.residues.iter().map(|r| r.atom_indices.len()).sum()
    }

    /// Returns the number of residues in this chain.
    pub fn residue_count(&self) -> usize {
        self.residues.len()
    }

    /// Finds a residue by its sequence number.
    pub fn get_residue(&self, seq: i32) -> Option<&Residue> {
        self.residues.iter().find(|r| r.seq == seq)
    }

    /// Finds a mutable residue by its sequence number.
    pub fn get_residue_mut(&mut self, seq: i32) -> Option<&mut Residue> {
        self.residues.iter_mut().find(|r| r.seq == seq)
    }
}

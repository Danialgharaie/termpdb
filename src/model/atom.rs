//! Atom data structure and helper predicates.

use crate::math::Vec3;
use crate::model::elements::Element;

/// Represents a single atom in a molecular structure.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    /// 0-based sequential index within the structure's atom array
    pub index: usize,
    /// Atom serial number from PDB/CIF file
    pub serial: i32,
    /// Atom name (e.g. "CA", "N", "C", "O", "CB", "P")
    pub name: String,
    /// Chemical element metadata
    pub element: Element,
    /// 3D position in Angstroms (Å)
    pub pos: Vec3,
    /// Temperature factor (B-factor) or pLDDT confidence score
    pub b_factor: f32,
    /// Crystallographic occupancy (typically 0.0 .. 1.0, default 1.0)
    pub occupancy: f32,
    /// Parent residue 3-letter / 1-letter name (e.g. "ALA", "GLY", "DA")
    pub res_name: String,
    /// Parent residue sequence number
    pub res_seq: i32,
    /// Chain identifier (e.g. "A", "B")
    pub chain_id: String,
    /// True for heteroatoms (ligands, ions, water), false for standard polymer residues
    pub is_hetatm: bool,
    /// Alternate location indicator if present (e.g. 'A', 'B')
    pub alt_loc: Option<char>,
    /// Formal electrostatic charge if specified
    pub charge: Option<i8>,
}

impl Atom {
    /// Creates a new `Atom` with standard defaults for optional fields.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        index: usize,
        serial: i32,
        name: impl Into<String>,
        element: Element,
        pos: Vec3,
        b_factor: f32,
        res_name: impl Into<String>,
        res_seq: i32,
        chain_id: impl Into<String>,
        is_hetatm: bool,
    ) -> Self {
        Self {
            index,
            serial,
            name: name.into(),
            element,
            pos,
            b_factor,
            occupancy: 1.0,
            res_name: res_name.into(),
            res_seq,
            chain_id: chain_id.into(),
            is_hetatm,
            alt_loc: None,
            charge: None,
        }
    }

    /// Returns `true` if this atom is a protein C-alpha (CA) backbone carbon.
    pub fn is_c_alpha(&self) -> bool {
        if self.is_hetatm {
            return false;
        }
        let trimmed = self.name.trim();
        trimmed.eq_ignore_ascii_case("CA")
            && (self.element.symbol.eq_ignore_ascii_case("C") || self.element.atomic_number == 6)
    }

    /// Returns `true` if this atom is part of a standard protein or nucleic acid backbone.
    pub fn is_backbone(&self) -> bool {
        if self.is_hetatm {
            return false;
        }
        let trimmed = self.name.trim();
        // Protein backbone
        if trimmed.eq_ignore_ascii_case("N")
            || trimmed.eq_ignore_ascii_case("CA")
            || trimmed.eq_ignore_ascii_case("C")
            || trimmed.eq_ignore_ascii_case("O")
            || trimmed.eq_ignore_ascii_case("OXT")
        {
            return true;
        }

        // Nucleic acid backbone
        matches!(
            trimmed,
            "P" | "OP1"
                | "OP2"
                | "O5'"
                | "O5*"
                | "C5'"
                | "C5*"
                | "C4'"
                | "C4*"
                | "O4'"
                | "O4*"
                | "C3'"
                | "C3*"
                | "O3'"
                | "O3*"
        )
    }

    /// Returns `true` if this atom is a hydrogen atom.
    pub fn is_hydrogen(&self) -> bool {
        self.element.atomic_number == 1
            || self.element.symbol.eq_ignore_ascii_case("H")
            || self.name.trim().starts_with('H')
    }

    /// Returns the covalent radius in Å.
    pub fn covalent_radius(&self) -> f32 {
        self.element.covalent_radius
    }

    /// Returns the Van der Waals radius in Å.
    pub fn vdw_radius(&self) -> f32 {
        self.element.vdw_radius
    }
}

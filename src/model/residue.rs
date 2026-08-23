//! Residue representations, amino acid classification, and secondary structure.

use crate::model::atom::Atom;

/// Secondary structure conformation for a residue or segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum SecondaryStructure {
    /// Alpha-helix or other helical conformation (3-10 helix, pi-helix)
    Helix,
    /// Beta-strand / Beta-sheet
    Sheet,
    /// Random coil, turn, loop, or unspecified conformation
    #[default]
    Coil,
}

impl SecondaryStructure {
    /// Returns `true` if this secondary structure is a helix.
    pub fn is_helix(&self) -> bool {
        matches!(self, SecondaryStructure::Helix)
    }

    /// Returns `true` if this secondary structure is a beta-sheet/strand.
    pub fn is_sheet(&self) -> bool {
        matches!(self, SecondaryStructure::Sheet)
    }

    /// Returns `true` if this secondary structure is a coil or loop.
    pub fn is_coil(&self) -> bool {
        matches!(self, SecondaryStructure::Coil)
    }
}

/// Represents a residue (amino acid, nucleotide, ligand, or water molecule) in a chain.
#[derive(Debug, Clone, PartialEq)]
pub struct Residue {
    /// Residue sequence number (res_seq)
    pub seq: i32,
    /// Insertion code if present (e.g. 'A')
    pub ins_code: Option<char>,
    /// Residue name (e.g. "ALA", "GLY", "DA", "HOH")
    pub name: String,
    /// Chain identifier
    pub chain_id: String,
    /// Indices of atoms belonging to this residue in the parent `Structure.atoms` array
    pub atom_indices: Vec<usize>,
    /// Assigned secondary structure
    pub secondary_structure: SecondaryStructure,
}

impl Residue {
    /// Creates a new residue with default coil secondary structure and empty atom list.
    pub fn new(seq: i32, name: impl Into<String>, chain_id: impl Into<String>) -> Self {
        Self {
            seq,
            ins_code: None,
            name: name.into(),
            chain_id: chain_id.into(),
            atom_indices: Vec::new(),
            secondary_structure: SecondaryStructure::Coil,
        }
    }

    /// Returns `true` if this residue is a standard or common modified amino acid.
    pub fn is_amino_acid(&self) -> bool {
        let trimmed = self.name.trim().to_ascii_uppercase();
        matches!(
            trimmed.as_str(),
            "ALA"
                | "ARG"
                | "ASN"
                | "ASP"
                | "CYS"
                | "GLN"
                | "GLU"
                | "GLY"
                | "HIS"
                | "ILE"
                | "LEU"
                | "LYS"
                | "MET"
                | "PHE"
                | "PRO"
                | "SER"
                | "THR"
                | "TRP"
                | "TYR"
                | "VAL"
                | "MSE"
                | "SEC"
                | "PYL"
                | "ASX"
                | "GLX"
                | "XLE"
                | "HYP"
                | "PCA"
        )
    }

    /// Returns `true` if this residue is a nucleic acid nucleotide (RNA / DNA).
    pub fn is_nucleic(&self) -> bool {
        let trimmed = self.name.trim().to_ascii_uppercase();
        matches!(
            trimmed.as_str(),
            "A" | "C"
                | "G"
                | "T"
                | "U"
                | "DA"
                | "DC"
                | "DG"
                | "DT"
                | "DI"
                | "DU"
                | "ADE"
                | "CYT"
                | "GUA"
                | "THY"
                | "URA"
                | "URI"
                | "+A"
                | "+C"
                | "+G"
                | "+T"
                | "+U"
        )
    }

    /// Returns `true` if `name` is a solvent / water residue name.
    pub fn name_is_water(name: &str) -> bool {
        matches!(
            name.trim().to_ascii_uppercase().as_str(),
            "HOH" | "WAT" | "H2O" | "DOD" | "TIP3" | "SOL"
        )
    }

    /// Returns `true` if this residue is a solvent / water molecule.
    pub fn is_water(&self) -> bool {
        Self::name_is_water(&self.name)
    }

    /// Returns the C-alpha (CA) atom reference if present in `atoms`.
    pub fn ca_atom<'a>(&self, atoms: &'a [Atom]) -> Option<&'a Atom> {
        self.atom_indices
            .iter()
            .filter_map(|&idx| atoms.get(idx))
            .find(|atom| atom.is_c_alpha())
    }

    /// Returns the structure atom index of the C-alpha (CA) atom if present.
    pub fn ca_atom_index(&self, atoms: &[Atom]) -> Option<usize> {
        self.ca_atom(atoms).map(|a| a.index)
    }

    /// Returns the single-letter code for amino acids / nucleotides, or 'X' for unknown.
    pub fn one_letter_code(&self) -> char {
        let trimmed = self.name.trim().to_ascii_uppercase();
        match trimmed.as_str() {
            "ALA" => 'A',
            "ARG" => 'R',
            "ASN" => 'N',
            "ASP" => 'D',
            "CYS" => 'C',
            "GLN" => 'Q',
            "GLU" => 'E',
            "GLY" => 'G',
            "HIS" => 'H',
            "ILE" => 'I',
            "LEU" => 'L',
            "LYS" => 'K',
            "MET" | "MSE" => 'M',
            "PHE" => 'F',
            "PRO" => 'P',
            "SER" => 'S',
            "THR" => 'T',
            "TRP" => 'W',
            "TYR" => 'Y',
            "VAL" => 'V',
            "SEC" => 'U',
            "PYL" => 'O',
            "A" | "DA" | "ADE" => 'A',
            "C" | "DC" | "CYT" => 'C',
            "G" | "DG" | "GUA" => 'G',
            "T" | "DT" | "THY" => 'T',
            "U" | "DU" | "URA" | "URI" => 'U',
            _ => 'X',
        }
    }

    /// Returns the Kyte-Doolittle hydrophobicity score (-4.5 to +4.5).
    /// Positive values are hydrophobic, negative values are hydrophilic.
    pub fn hydrophobicity_score(&self) -> f32 {
        let trimmed = self.name.trim().to_ascii_uppercase();
        match trimmed.as_str() {
            "ILE" => 4.5,
            "VAL" => 4.2,
            "LEU" => 3.8,
            "PHE" => 2.8,
            "CYS" => 2.5,
            "MET" | "MSE" => 1.9,
            "ALA" => 1.8,
            "GLY" => -0.4,
            "THR" => -0.7,
            "SER" => -0.8,
            "TRP" => -0.9,
            "TYR" => -1.3,
            "PRO" => -1.6,
            "HIS" => -3.2,
            "GLU" | "GLX" => -3.5,
            "GLN" => -3.5,
            "ASP" | "ASX" => -3.5,
            "ASN" => -3.5,
            "LYS" => -3.9,
            "ARG" => -4.5,
            _ => 0.0,
        }
    }
}

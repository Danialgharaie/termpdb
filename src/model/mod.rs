//! Chemical and molecular data models.
//!
//! Provides the core domain representation for atoms, residues, chains, bonds,
//! elements, and complete macromolecular structures.

pub mod align;
pub mod assembly;
pub mod atom;
pub mod bond;
pub mod chain;
pub mod dssp;
pub mod elements;
pub mod geometry;
pub mod interactions;
pub mod residue;
pub mod spatial;

use std::collections::HashMap;

pub use align::{SuperpositionResult, needleman_wunsch, pair_ca_coordinates, superimpose_structures};
pub use assembly::{Assembly, AssemblyGen, expand_model};
pub use atom::Atom;
pub use bond::{Bond, BondDetector, BondOrder};
pub use chain::Chain;
pub use dssp::{assign_dssp, calculate_dssp_hbond_energy};
pub use elements::{ELEMENTS, Element, element_by_atomic_number, element_by_symbol};
pub use geometry::{RamachandranRegion, calculate_bond_angle, calculate_dihedral_angle, classify_ramachandran};
pub use interactions::{Interaction, InteractionKind, detect_interactions};
pub use residue::{Residue, SecondaryStructure};
pub use spatial::SpatialGrid;

use crate::error::{Result, TermPdbError};
use crate::math::Vec3;

/// One NMR / docking / deposition model: a full coordinate set with its own
/// atoms, chains, and bonds. Bonds never cross models.
#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    /// File serial (`MODEL` / `pdbx_PDB_model_num`), not a compacted 0-based index.
    pub serial: i32,
    /// Polymer chains in this model
    pub chains: Vec<Chain>,
    /// Atoms in this model (indices are local to this model)
    pub atoms: Vec<Atom>,
    /// Covalent bonds connecting atoms in this model
    pub bonds: Vec<Bond>,
}

impl Model {
    /// Creates an empty model with the given file serial.
    pub fn new(serial: i32) -> Self {
        Self {
            serial,
            chains: Vec::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }
}

/// Complete molecular structure containing one or more models plus metadata.
///
/// Accessors such as [`Structure::atoms`] refer to the **active** model.
/// Files without `MODEL` records become a single model with serial `1`.
#[derive(Debug, Clone, PartialEq)]
pub struct Structure {
    /// Structure title or description (e.g. from PDB TITLE record)
    pub title: String,
    /// 4-letter PDB ID or identifier if known (e.g. "1CRN")
    pub id_code: Option<String>,
    models: Vec<Model>,
    active_index: usize,
    assemblies: Vec<Assembly>,
    /// `None` = deposited asymmetric unit. `Some(i)` indexes `assemblies`.
    active_assembly: Option<usize>,
    expanded: Option<Model>,
    /// Arbitrary key-value metadata (e.g. resolution, experimental method, deposition date)
    pub metadata: HashMap<String, String>,
}

impl Structure {
    /// Creates a new empty `Structure` with the given title and a single model (serial 1).
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            id_code: None,
            models: vec![Model::new(1)],
            active_index: 0,
            assemblies: Vec::new(),
            active_assembly: None,
            expanded: None,
            metadata: HashMap::new(),
        }
    }

    /// Creates a new empty `Structure` with ID code and title.
    pub fn with_id(id_code: impl Into<String>, title: impl Into<String>) -> Self {
        Self {
            title: title.into(),
            id_code: Some(id_code.into()),
            models: vec![Model::new(1)],
            active_index: 0,
            assemblies: Vec::new(),
            active_assembly: None,
            expanded: None,
            metadata: HashMap::new(),
        }
    }

    /// Replaces all models. Active model becomes the lowest serial (index 0 after sort).
    /// An empty list becomes a single empty model with serial 1.
    pub fn set_models(&mut self, mut models: Vec<Model>) {
        models.sort_by_key(|m| m.serial);
        if models.is_empty() {
            self.models = vec![Model::new(1)];
            self.active_index = 0;
        } else {
            self.models = models;
            self.active_index = 0;
        }
        self.refresh_view();
    }

    /// Replaces biological assemblies. Active view becomes the asymmetric unit.
    pub fn set_assemblies(&mut self, assemblies: Vec<Assembly>) {
        self.assemblies = assemblies;
        self.active_assembly = None;
        self.expanded = None;
    }

    fn refresh_view(&mut self) {
        self.expanded = None;
        let Some(i) = self.active_assembly else {
            return;
        };
        if i >= self.assemblies.len() || self.models.is_empty() {
            self.active_assembly = None;
            return;
        }
        let expanded = expand_model(&self.models[self.active_index], &self.assemblies[i]);
        match expanded {
            Ok(model) => self.expanded = Some(model),
            Err(_) => {
                self.active_assembly = None;
                self.expanded = None;
            }
        }
    }

    fn view_model(&self) -> &Model {
        self.expanded
            .as_ref()
            .unwrap_or(&self.models[self.active_index])
    }

    /// All models, sorted by file serial.
    pub fn models(&self) -> &[Model] {
        &self.models
    }

    /// File serials in ascending order.
    pub fn model_serials(&self) -> Vec<i32> {
        self.models.iter().map(|m| m.serial).collect()
    }

    /// Number of models in the file.
    pub fn model_count(&self) -> usize {
        self.models.len()
    }

    /// True when the file contains more than one model.
    pub fn has_multiple_models(&self) -> bool {
        self.models.len() > 1
    }

    /// Highest file serial among models (denominator for `Model 5/20` HUD).
    pub fn max_model_serial(&self) -> i32 {
        self.models.last().map(|m| m.serial).unwrap_or(1)
    }

    /// File serial of the active model.
    pub fn active_model_serial(&self) -> i32 {
        self.active_model().serial
    }

    /// The currently active model.
    pub fn active_model(&self) -> &Model {
        self.view_model()
    }

    /// The currently active model, mutably.
    pub fn active_model_mut(&mut self) -> &mut Model {
        &mut self.models[self.active_index]
    }

    /// Makes `serial` the active model. Errors if that serial is not in the file.
    pub fn set_active_model(&mut self, serial: i32) -> Result<()> {
        match self.models.iter().position(|m| m.serial == serial) {
            Some(idx) => {
                self.active_index = idx;
                self.refresh_view();
                Ok(())
            }
            None => {
                let available = self
                    .model_serials()
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>()
                    .join(", ");
                Err(TermPdbError::InvalidStructure(format!(
                    "Model {} not found (available: {})",
                    serial, available
                )))
            }
        }
    }

    /// Steps to the next model by serial order, wrapping to the first.
    pub fn next_model(&mut self) {
        let n = self.models.len();
        if n == 0 {
            return;
        }
        self.active_index = (self.active_index + 1) % n;
        self.refresh_view();
    }

    /// Steps to the previous model by serial order, wrapping to the last.
    pub fn prev_model(&mut self) {
        let n = self.models.len();
        if n == 0 {
            return;
        }
        self.active_index = (self.active_index + n - 1) % n;
        self.refresh_view();
    }

    /// Parsed biological assemblies (empty if the file had none).
    pub fn assemblies(&self) -> &[Assembly] {
        &self.assemblies
    }

    pub fn has_assemblies(&self) -> bool {
        !self.assemblies.is_empty()
    }

    /// File id of the active assembly, or `None` for the asymmetric unit.
    pub fn active_assembly_id(&self) -> Option<&str> {
        self.active_assembly
            .and_then(|i| self.assemblies.get(i).map(|a| a.id.as_str()))
    }

    pub fn assembly_ids(&self) -> Vec<&str> {
        self.assemblies.iter().map(|a| a.id.as_str()).collect()
    }

    /// `None`/`"asu"`/`"au"`/`"0"` selects the asymmetric unit. Other values match file ids.
    pub fn set_assembly(&mut self, id: Option<&str>) -> Result<()> {
        let Some(id) = id else {
            self.active_assembly = None;
            self.expanded = None;
            return Ok(());
        };
        let trimmed = id.trim();
        if trimmed.is_empty()
            || trimmed.eq_ignore_ascii_case("asu")
            || trimmed.eq_ignore_ascii_case("au")
            || trimmed == "0"
        {
            self.active_assembly = None;
            self.expanded = None;
            return Ok(());
        }
        match self
            .assemblies
            .iter()
            .position(|a| a.id.eq_ignore_ascii_case(trimmed))
        {
            Some(idx) => {
                self.active_assembly = Some(idx);
                self.refresh_view();
                Ok(())
            }
            None => {
                let available = self.assembly_ids().join(", ");
                Err(TermPdbError::InvalidStructure(format!(
                    "Assembly '{}' not found (available: {available})",
                    trimmed
                )))
            }
        }
    }

    /// ASU → first assembly → … → wrap to ASU.
    pub fn next_assembly(&mut self) {
        if self.assemblies.is_empty() {
            return;
        }
        self.active_assembly = match self.active_assembly {
            None => Some(0),
            Some(i) if i + 1 < self.assemblies.len() => Some(i + 1),
            Some(_) => None,
        };
        self.refresh_view();
    }

    pub fn prev_assembly(&mut self) {
        if self.assemblies.is_empty() {
            return;
        }
        self.active_assembly = match self.active_assembly {
            None => Some(self.assemblies.len() - 1),
            Some(0) => None,
            Some(i) => Some(i - 1),
        };
        self.refresh_view();
    }

    /// Atoms of the active view (assembly expansion or deposited model).
    pub fn atoms(&self) -> &[Atom] {
        &self.view_model().atoms
    }

    /// Mutable atoms of the deposited active model (invalidates a live assembly view).
    pub fn atoms_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.models[self.active_index].atoms
    }

    /// Chains of the active view.
    pub fn chains(&self) -> &[Chain] {
        &self.view_model().chains
    }

    /// Mutable chains of the deposited active model.
    pub fn chains_mut(&mut self) -> &mut Vec<Chain> {
        &mut self.models[self.active_index].chains
    }

    /// Bonds of the active view.
    pub fn bonds(&self) -> &[Bond] {
        &self.view_model().bonds
    }

    /// Mutable bonds of the deposited active model.
    pub fn bonds_mut(&mut self) -> &mut Vec<Bond> {
        &mut self.models[self.active_index].bonds
    }

    /// Ensures the active view has covalent bonds, detecting them once with an
    /// O(N) spatial hash if none were provided (e.g. a PDB/mmCIF file with no
    /// CONECT records viewed as the asymmetric unit). Assembly-expanded views
    /// already have bonds detected during expansion, so this is a no-op there.
    ///
    /// Bond detection is relatively expensive, so callers should invoke this
    /// once when the view is established (load, model switch, assembly switch)
    /// rather than letting the renderer re-detect on every frame.
    pub fn ensure_bonds(&mut self) {
        if self.bonds().is_empty() && !self.atoms().is_empty() {
            let detected = BondDetector::detect_bonds(self.atoms());
            *self.bonds_mut() = detected;
        }
    }

    /// Adds an atom to the active model, assigning its sequential index.
    /// Returns the index of the newly added atom.
    pub fn add_atom(&mut self, mut atom: Atom) -> usize {
        let idx = self.atoms().len();
        atom.index = idx;
        self.atoms_mut().push(atom);
        self.refresh_view();
        idx
    }

    /// Adds a chain to the active model.
    pub fn add_chain(&mut self, chain: Chain) {
        self.chains_mut().push(chain);
        self.refresh_view();
    }

    /// Adds a bond to the active model.
    pub fn add_bond(&mut self, bond: Bond) {
        self.bonds_mut().push(bond);
        self.refresh_view();
    }

    /// Total number of atoms in the active model.
    pub fn atom_count(&self) -> usize {
        self.atoms().len()
    }

    /// Total number of chains in the active model.
    pub fn chain_count(&self) -> usize {
        self.chains().len()
    }

    /// Total number of residues across all chains in the active model.
    pub fn residue_count(&self) -> usize {
        self.chains().iter().map(|c| c.residue_count()).sum()
    }

    /// Finds a chain by its ID string in the active model.
    pub fn get_chain(&self, id: &str) -> Option<&Chain> {
        self.chains().iter().find(|c| c.id == id)
    }

    /// Finds a mutable chain by its ID string in the active model.
    pub fn get_chain_mut(&mut self, id: &str) -> Option<&mut Chain> {
        self.chains_mut().iter_mut().find(|c| c.id == id)
    }

    /// Returns references to all C-alpha (CA) atoms in the active model in sequence order.
    pub fn ca_atoms(&self) -> Vec<&Atom> {
        let atoms = self.atoms();
        self.chains()
            .iter()
            .flat_map(|c| c.ca_atoms(atoms))
            .collect()
    }

    /// Computes the geometric center of mass (centroid) of the active model.
    pub fn center_of_mass(&self) -> Vec3 {
        let atoms = self.atoms();
        if atoms.is_empty() {
            return Vec3::ZERO;
        }

        let mut sum = Vec3::ZERO;
        for atom in atoms {
            sum += atom.pos;
        }

        sum / (atoms.len() as f32)
    }

    /// Computes the radius of the minimum bounding sphere of the active model.
    pub fn bounding_sphere_radius(&self) -> f32 {
        let atoms = self.atoms();
        if atoms.is_empty() {
            return 1.0;
        }

        let com = self.center_of_mass();
        let mut max_dist_sq: f32 = 0.0;

        for atom in atoms {
            let dist_sq = atom.pos.distance_squared(&com);
            if dist_sq > max_dist_sq {
                max_dist_sq = dist_sq;
            }
        }

        let r = max_dist_sq.sqrt();
        if r < 1e-4 { 1.0 } else { r }
    }

    /// Translates active-model atom positions so the center of mass is at `(0, 0, 0)`.
    pub fn center_and_normalize(&mut self) {
        if self.atoms().is_empty() {
            return;
        }

        let com = {
            let atoms = &self.models[self.active_index].atoms;
            if atoms.is_empty() {
                return;
            }
            let mut sum = Vec3::ZERO;
            for atom in atoms {
                sum += atom.pos;
            }
            sum / (atoms.len() as f32)
        };
        for atom in &mut self.models[self.active_index].atoms {
            atom.pos -= com;
        }
        self.refresh_view();
    }

    /// Automatically detects and builds covalent bonds for the deposited active model.
    pub fn build_bonds(&mut self) {
        let bonds = BondDetector::detect_bonds(&self.models[self.active_index].atoms);
        self.models[self.active_index].bonds = bonds;
        self.refresh_view();
    }

    /// Returns the minimum and maximum B-factor in the active model (or `(0.0, 100.0)` if empty).
    pub fn b_factor_range(&self) -> (f32, f32) {
        let atoms = self.atoms();
        if atoms.is_empty() {
            return (0.0, 100.0);
        }

        let mut min_b = f32::INFINITY;
        let mut max_b = f32::NEG_INFINITY;

        for atom in atoms {
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

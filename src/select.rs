//! Atom specifiers, selection (at most two atoms), and distance.

use crate::error::{Result, TermPdbError};
use crate::model::{Atom, Residue, Structure};
use crate::render::representations::build_atom_residue_map;
use crate::render::{Camera, Visibility};

/// Parsed atom identifier: optional chain, residue sequence, optional atom name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomSpec {
    pub chain_id: Option<String>,
    pub res_seq: i32,
    pub atom_name: Option<String>,
}

impl AtomSpec {
    pub fn new(chain_id: Option<&str>, res_seq: i32, atom_name: Option<&str>) -> Self {
        Self {
            chain_id: chain_id.map(|s| s.to_string()),
            res_seq,
            atom_name: atom_name.map(|s| s.to_string()),
        }
    }
}

/// At most two selected atom indices in the active model.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Selection {
    indices: Vec<usize>,
}

impl Selection {
    pub fn atoms(&self) -> &[usize] {
        &self.indices
    }

    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    pub fn clear(&mut self) {
        self.indices.clear();
    }

    /// Adds `idx`. Picking an already-selected atom deselects it. A 5th pick
    /// drops the oldest so a walking selection stays at most 4 atoms.
    pub fn pick(&mut self, idx: usize) {
        if let Some(pos) = self.indices.iter().position(|&i| i == idx) {
            self.indices.remove(pos);
            return;
        }
        if self.indices.len() >= 4 {
            self.indices.remove(0);
        }
        self.indices.push(idx);
    }

    /// Computes Euclidean distance between the first two selected atoms.
    pub fn distance(&self, structure: &Structure) -> Option<f32> {
        if self.indices.len() < 2 {
            return None;
        }
        atom_distance(structure, self.indices[0], self.indices[1]).ok()
    }

    /// Computes planar bond angle between the first three selected atoms.
    pub fn angle(&self, structure: &Structure) -> Option<f32> {
        if self.indices.len() < 3 {
            return None;
        }
        let atoms = structure.atoms();
        let p1 = atoms.get(self.indices[0])?.pos;
        let p2 = atoms.get(self.indices[1])?.pos;
        let p3 = atoms.get(self.indices[2])?.pos;
        Some(crate::model::geometry::calculate_bond_angle(p1, p2, p3))
    }

    /// Computes dihedral / torsion angle between 4 selected atoms.
    pub fn dihedral(&self, structure: &Structure) -> Option<f32> {
        if self.indices.len() < 4 {
            return None;
        }
        let atoms = structure.atoms();
        let p1 = atoms.get(self.indices[0])?.pos;
        let p2 = atoms.get(self.indices[1])?.pos;
        let p3 = atoms.get(self.indices[2])?.pos;
        let p4 = atoms.get(self.indices[3])?.pos;
        Some(crate::model::geometry::calculate_dihedral_angle(
            p1, p2, p3, p4,
        ))
    }

    /// Returns a full formatted geometry status description for the HUD status.
    pub fn status(&self, structure: &Structure) -> String {
        match self.indices.len() {
            0 => "No selection".to_string(),
            1 => format!("Selected: {}", atom_label(structure, self.indices[0])),
            2 => {
                let a = atom_label(structure, self.indices[0]);
                let b = atom_label(structure, self.indices[1]);
                if let Some(d) = self.distance(structure) {
                    format!("Distance: {a} · {b} = {d:.2} Å")
                } else {
                    format!("{a} · {b}")
                }
            }
            3 => {
                let a = atom_label(structure, self.indices[0]);
                let b = atom_label(structure, self.indices[1]);
                let c = atom_label(structure, self.indices[2]);
                if let Some(ang) = self.angle(structure) {
                    format!("Angle: {ang:.1}° ({a} · {b} · {c})")
                } else {
                    format!("{a} · {b} · {c}")
                }
            }
            _ => {
                let a = atom_label(structure, self.indices[0]);
                let b = atom_label(structure, self.indices[1]);
                let c = atom_label(structure, self.indices[2]);
                let d = atom_label(structure, self.indices[3]);
                if let Some(dih) = self.dihedral(structure) {
                    let rama = crate::model::geometry::classify_ramachandran(dih, 0.0);
                    format!(
                        "Dihedral: {dih:.1}° [{}] ({a} · {b} · {c} · {d})",
                        rama.name()
                    )
                } else {
                    format!("{a} · {b} · {c} · {d}")
                }
            }
        }
    }

    pub fn status_line(&self, structure: &Structure) -> Option<String> {
        if self.indices.is_empty() {
            return None;
        }
        let a = atom_label(structure, self.indices[0]);
        if self.indices.len() == 1 {
            return Some(a);
        }
        let b = atom_label(structure, self.indices[1]);
        if self.indices.len() == 2 {
            if let Some(d) = self.distance(structure) {
                return Some(format!("{a} · {b}  {d:.2} Å"));
            } else {
                return Some(format!("{a} · {b}"));
            }
        }
        if self.indices.len() == 3 {
            let c = atom_label(structure, self.indices[2]);
            if let Some(ang) = self.angle(structure) {
                return Some(format!("{a} · {b} · {c}  {ang:.1}°"));
            } else {
                return Some(format!("{a} · {b} · {c}"));
            }
        }
        let c = atom_label(structure, self.indices[2]);
        let d = atom_label(structure, self.indices[3]);
        if let Some(dih) = self.dihedral(structure) {
            let rama = crate::model::geometry::classify_ramachandran(dih, 0.0);
            Some(format!(
                "{a} · {b} · {c} · {d}  {dih:.1}° [{}]",
                rama.name()
            ))
        } else {
            Some(format!("{a} · {b} · {c} · {d}"))
        }
    }
}

/// Parses `A:12:CA`, `A:12`, `12`, `A/12/N`, or `A 12 CA`.
pub fn parse_atom_spec(input: &str) -> Result<AtomSpec> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Err(TermPdbError::ParseError("Empty atom specifier".to_string()));
    }

    let normalized = trimmed.replace('/', ":");
    let parts: Vec<&str> = if normalized.contains(':') {
        normalized.split(':').map(str::trim).collect()
    } else {
        normalized.split_whitespace().collect()
    };

    if parts.is_empty() || parts.iter().any(|p| p.is_empty()) {
        return Err(TermPdbError::ParseError(format!(
            "Invalid atom specifier '{input}'"
        )));
    }

    match parts.len() {
        1 => {
            let res_seq = parse_seq(parts[0], input)?;
            Ok(AtomSpec::new(None, res_seq, None))
        }
        2 => {
            if let Ok(res_seq) = parts[0].parse::<i32>() {
                Ok(AtomSpec::new(None, res_seq, Some(parts[1])))
            } else {
                let res_seq = parse_seq(parts[1], input)?;
                Ok(AtomSpec::new(Some(parts[0]), res_seq, None))
            }
        }
        3 => {
            let res_seq = parse_seq(parts[1], input)?;
            Ok(AtomSpec::new(Some(parts[0]), res_seq, Some(parts[2])))
        }
        _ => Err(TermPdbError::ParseError(format!(
            "Invalid atom specifier '{input}'"
        ))),
    }
}

fn parse_seq(s: &str, original: &str) -> Result<i32> {
    s.parse::<i32>()
        .map_err(|_| TermPdbError::ParseError(format!("Invalid residue number in '{original}'")))
}

/// Resolves a specifier to an atom index in the active model.
///
/// Explicit atom names are found even if currently hidden. Residue-only specs
/// skip hidden atoms and prefer CA, then P, then the first visible heavy atom.
pub fn resolve_atom(
    structure: &Structure,
    spec: &AtomSpec,
    visibility: Option<&Visibility>,
) -> Result<usize> {
    let mut found_residue: Option<&Residue> = None;

    for chain in structure.chains() {
        if let Some(id) = &spec.chain_id
            && !chain.id.eq_ignore_ascii_case(id)
        {
            continue;
        }
        if let Some(res) = chain.get_residue(spec.res_seq) {
            found_residue = Some(res);
            break;
        }
        if spec.chain_id.is_some() {
            break;
        }
    }

    let Some(res) = found_residue else {
        let where_ = match &spec.chain_id {
            Some(c) => format!("{c}:{}", spec.res_seq),
            None => spec.res_seq.to_string(),
        };
        return Err(TermPdbError::InvalidStructure(format!(
            "Residue {where_} not found"
        )));
    };

    let atoms = structure.atoms();
    if let Some(name) = &spec.atom_name {
        for &idx in &res.atom_indices {
            if let Some(atom) = atoms.get(idx)
                && atom.name.trim().eq_ignore_ascii_case(name.trim())
            {
                return Ok(idx);
            }
        }
        return Err(TermPdbError::InvalidStructure(format!(
            "Atom '{}' not found in {}{}",
            name,
            spec.chain_id
                .as_ref()
                .map(|c| format!("{c}:"))
                .unwrap_or_default(),
            spec.res_seq
        )));
    }

    if let Some(ca) = res.ca_atom(atoms)
        && atom_allowed(ca, Some(res), visibility)
    {
        return Ok(ca.index);
    }

    for &idx in &res.atom_indices {
        if let Some(atom) = atoms.get(idx)
            && atom.name.trim().eq_ignore_ascii_case("P")
            && atom_allowed(atom, Some(res), visibility)
        {
            return Ok(idx);
        }
    }

    for &idx in &res.atom_indices {
        if let Some(atom) = atoms.get(idx)
            && !atom.is_hydrogen()
            && atom_allowed(atom, Some(res), visibility)
        {
            return Ok(idx);
        }
    }

    for &idx in &res.atom_indices {
        if let Some(atom) = atoms.get(idx)
            && atom_allowed(atom, Some(res), visibility)
        {
            return Ok(idx);
        }
    }

    if visibility.is_none()
        && let Some(&idx) = res.atom_indices.first()
    {
        return Ok(idx);
    }

    Err(TermPdbError::InvalidStructure(format!(
        "Residue {} has no visible atoms",
        spec.res_seq
    )))
}

fn atom_allowed(atom: &Atom, residue: Option<&Residue>, visibility: Option<&Visibility>) -> bool {
    match visibility {
        Some(v) => v.atom_visible(atom, residue),
        None => true,
    }
}

/// Parses `SPEC,SPEC` and returns `A:1:CA  A:2:N  3.824`.
pub fn distance_report(structure: &Structure, pair: &str) -> Result<String> {
    let parts: Vec<&str> = pair
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect();
    if parts.len() != 2 {
        return Err(TermPdbError::ParseError(
            "Expected two atom specifiers separated by a comma (e.g. A:12:CA,A:40:N)".to_string(),
        ));
    }
    let a = parse_atom_spec(parts[0])?;
    let b = parse_atom_spec(parts[1])?;
    let i = resolve_atom(structure, &a, None)?;
    let j = resolve_atom(structure, &b, None)?;
    let d = atom_distance(structure, i, j)?;
    Ok(format!(
        "{}  {}  {:.3}",
        atom_label(structure, i),
        atom_label(structure, j),
        d
    ))
}

/// Parses `SPEC,SPEC,SPEC` and returns `A:1:CA  A:2:CA  A:3:CA  109.50°`.
pub fn angle_report(structure: &Structure, triplet: &str) -> Result<String> {
    let parts: Vec<&str> = triplet
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect();
    if parts.len() != 3 {
        return Err(TermPdbError::ParseError(
            "Expected three atom specifiers separated by commas (e.g. A:1:CA,A:2:CA,A:3:CA)"
                .to_string(),
        ));
    }
    let a = parse_atom_spec(parts[0])?;
    let b = parse_atom_spec(parts[1])?;
    let c = parse_atom_spec(parts[2])?;
    let i = resolve_atom(structure, &a, None)?;
    let j = resolve_atom(structure, &b, None)?;
    let k = resolve_atom(structure, &c, None)?;
    let atoms = structure.atoms();
    let ang =
        crate::model::geometry::calculate_bond_angle(atoms[i].pos, atoms[j].pos, atoms[k].pos);
    Ok(format!(
        "{}  {}  {}  {:.2}°",
        atom_label(structure, i),
        atom_label(structure, j),
        atom_label(structure, k),
        ang
    ))
}

/// Parses `SPEC,SPEC,SPEC,SPEC` and returns dihedral report with Ramachandran region.
pub fn dihedral_report(structure: &Structure, quartet: &str) -> Result<String> {
    let parts: Vec<&str> = quartet
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect();
    if parts.len() != 4 {
        return Err(TermPdbError::ParseError(
            "Expected four atom specifiers separated by commas (e.g. A:1:N,A:1:CA,A:1:C,A:2:N)"
                .to_string(),
        ));
    }
    let a = parse_atom_spec(parts[0])?;
    let b = parse_atom_spec(parts[1])?;
    let c = parse_atom_spec(parts[2])?;
    let d = parse_atom_spec(parts[3])?;
    let i = resolve_atom(structure, &a, None)?;
    let j = resolve_atom(structure, &b, None)?;
    let k = resolve_atom(structure, &c, None)?;
    let l = resolve_atom(structure, &d, None)?;
    let atoms = structure.atoms();
    let dih = crate::model::geometry::calculate_dihedral_angle(
        atoms[i].pos,
        atoms[j].pos,
        atoms[k].pos,
        atoms[l].pos,
    );
    let rama = crate::model::geometry::classify_ramachandran(dih, 0.0);
    Ok(format!(
        "{}  {}  {}  {}  {:.2}° [{}]",
        atom_label(structure, i),
        atom_label(structure, j),
        atom_label(structure, k),
        atom_label(structure, l),
        dih,
        rama.name()
    ))
}

/// Euclidean distance in Å between two atom indices in the active model.
pub fn atom_distance(structure: &Structure, i: usize, j: usize) -> Result<f32> {
    let atoms = structure.atoms();
    let a = atoms
        .get(i)
        .ok_or_else(|| TermPdbError::InvalidStructure(format!("Atom index {i} out of range")))?;
    let b = atoms
        .get(j)
        .ok_or_else(|| TermPdbError::InvalidStructure(format!("Atom index {j} out of range")))?;
    Ok(a.pos.distance(&b.pos))
}

/// `A:12:CA`-style label for HUD and CLI.
pub fn atom_label(structure: &Structure, idx: usize) -> String {
    match structure.atoms().get(idx) {
        Some(atom) => format!(
            "{}:{}:{}",
            if atom.chain_id.is_empty() {
                "A"
            } else {
                atom.chain_id.as_str()
            },
            atom.res_seq,
            atom.name.trim()
        ),
        None => format!("#{idx}"),
    }
}

/// Nearest visible atom to a framebuffer pixel, or `None` if nothing is close.
#[allow(clippy::too_many_arguments)]
pub fn pick_atom_at_screen(
    structure: &Structure,
    camera: &Camera,
    visibility: Visibility,
    width: usize,
    height: usize,
    sx: f32,
    sy: f32,
    max_radius: f32,
) -> Option<usize> {
    let residue_map = build_atom_residue_map(structure);
    let max_r_sq = max_radius * max_radius;
    let mut best: Option<(usize, f32, f32)> = None;

    for (idx, atom) in structure.atoms().iter().enumerate() {
        let residue = residue_map.get(&idx).copied();
        if !visibility.atom_visible(atom, residue) {
            continue;
        }
        let Some((px, py, depth)) = camera.world_to_screen(atom.pos, width, height) else {
            continue;
        };
        let dx = px - sx;
        let dy = py - sy;
        let dist_sq = dx * dx + dy * dy;
        if dist_sq > max_r_sq {
            continue;
        }
        let better = match best {
            None => true,
            Some((_, best_d, best_z)) => {
                dist_sq + 1e-3 < best_d || ((dist_sq - best_d).abs() < 1e-3 && depth < best_z)
            }
        };
        if better {
            best = Some((idx, dist_sq, depth));
        }
    }

    best.map(|(idx, _, _)| idx)
}

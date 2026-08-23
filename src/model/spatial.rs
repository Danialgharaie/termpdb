//! 3D Uniform Spatial Hash Grid for fast neighbor queries and occlusion culling.

use std::collections::HashMap;
use crate::math::Vec3;
use crate::model::atom::Atom;

/// A 3D spatial partitioning grid storing atom indices in cubic voxel bins.
#[derive(Debug, Clone)]
pub struct SpatialGrid<'a> {
    cell_size: f32,
    inv_cell_size: f32,
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
    atoms: &'a [Atom],
}

impl<'a> SpatialGrid<'a> {
    /// Builds a new `SpatialGrid` from a slice of atoms with the specified cubic voxel size.
    pub fn new(atoms: &'a [Atom], cell_size: f32) -> Self {
        let cell_size = if cell_size > 0.1 { cell_size } else { 3.5 };
        let inv_cell_size = 1.0 / cell_size;
        let mut cells: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::with_capacity(atoms.len() / 4);

        for (idx, atom) in atoms.iter().enumerate() {
            let key = (
                (atom.pos.x * inv_cell_size).floor() as i32,
                (atom.pos.y * inv_cell_size).floor() as i32,
                (atom.pos.z * inv_cell_size).floor() as i32,
            );
            cells.entry(key).or_default().push(idx);
        }

        Self {
            cell_size,
            inv_cell_size,
            cells,
            atoms,
        }
    }

    /// Returns the cubic cell size.
    pub fn cell_size(&self) -> f32 {
        self.cell_size
    }

    /// Returns the indices of all atoms located within Euclidean distance `radius` of `center`.
    pub fn neighbors_within(&self, center: Vec3, radius: f32) -> Vec<usize> {
        let mut result = Vec::new();
        let r_sq = radius * radius;

        let min_x = ((center.x - radius) * self.inv_cell_size).floor() as i32;
        let max_x = ((center.x + radius) * self.inv_cell_size).floor() as i32;
        let min_y = ((center.y - radius) * self.inv_cell_size).floor() as i32;
        let max_y = ((center.y + radius) * self.inv_cell_size).floor() as i32;
        let min_z = ((center.z - radius) * self.inv_cell_size).floor() as i32;
        let max_z = ((center.z + radius) * self.inv_cell_size).floor() as i32;

        for gx in min_x..=max_x {
            for gy in min_y..=max_y {
                for gz in min_z..=max_z {
                    if let Some(indices) = self.cells.get(&(gx, gy, gz)) {
                        for &idx in indices {
                            if idx < self.atoms.len() {
                                let dist_sq = self.atoms[idx].pos.distance_squared(&center);
                                if dist_sq <= r_sq {
                                    result.push(idx);
                                }
                            }
                        }
                    }
                }
            }
        }

        result
    }

    /// Computes a boolean mask of buried interior atoms (those with $\ge \text{min\_neighbors}$ within `cutoff_radius`).
    pub fn compute_buried_atoms(&self, cutoff_radius: f32, min_neighbors: usize) -> Vec<bool> {
        let mut buried = vec![false; self.atoms.len()];
        let r_sq = cutoff_radius * cutoff_radius;

        for (i, atom) in self.atoms.iter().enumerate() {
            let mut count = 0;

            let min_x = ((atom.pos.x - cutoff_radius) * self.inv_cell_size).floor() as i32;
            let max_x = ((atom.pos.x + cutoff_radius) * self.inv_cell_size).floor() as i32;
            let min_y = ((atom.pos.y - cutoff_radius) * self.inv_cell_size).floor() as i32;
            let max_y = ((atom.pos.y + cutoff_radius) * self.inv_cell_size).floor() as i32;
            let min_z = ((atom.pos.z - cutoff_radius) * self.inv_cell_size).floor() as i32;
            let max_z = ((atom.pos.z + cutoff_radius) * self.inv_cell_size).floor() as i32;

            'cell_loop: for gx in min_x..=max_x {
                for gy in min_y..=max_y {
                    for gz in min_z..=max_z {
                        if let Some(indices) = self.cells.get(&(gx, gy, gz)) {
                            for &other_idx in indices {
                                if other_idx != i && other_idx < self.atoms.len() {
                                    let dist_sq = self.atoms[other_idx].pos.distance_squared(&atom.pos);
                                    if dist_sq <= r_sq {
                                        count += 1;
                                        if count >= min_neighbors {
                                            buried[i] = true;
                                            break 'cell_loop;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        buried
    }
}

//! Contact Geometry: Bond Angles, Dihedral / Torsion Angles, and Ramachandran Conformations.

use crate::math::Vec3;

/// Ramachandran backbone conformation quadrant classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RamachandranRegion {
    /// Right-handed $\alpha$-helix region
    AlphaHelix,
    /// $\beta$-sheet / extended strand region
    BetaSheet,
    /// Left-handed $\alpha$-helix region
    LeftHandedAlpha,
    /// Disallowed or outlier region
    Outlier,
}

impl RamachandranRegion {
    /// Returns the human-readable display name of the Ramachandran region.
    pub fn name(&self) -> &'static str {
        match self {
            RamachandranRegion::AlphaHelix => "α-Helix (Favored)",
            RamachandranRegion::BetaSheet => "β-Sheet (Favored)",
            RamachandranRegion::LeftHandedAlpha => "Left-handed α",
            RamachandranRegion::Outlier => "Outlier / Disallowed",
        }
    }
}

/// Computes the planar bond angle at vertex $p_2$ formed by $p_1 - p_2 - p_3$ in degrees $\in [0^\circ, 180^\circ]$.
pub fn calculate_bond_angle(p1: Vec3, p2: Vec3, p3: Vec3) -> f32 {
    let v1 = p1 - p2;
    let v2 = p3 - p2;
    let len1 = v1.norm();
    let len2 = v2.norm();

    if len1 < 1e-6 || len2 < 1e-6 {
        return 0.0;
    }

    let dot = (v1.dot(v2) / (len1 * len2)).clamp(-1.0, 1.0);
    dot.acos().to_degrees()
}

/// Computes the dihedral (torsion) angle for sequence $p_1 - p_2 - p_3 - p_4$ in degrees $\in [-180^\circ, 180^\circ]$.
pub fn calculate_dihedral_angle(p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3) -> f32 {
    let b1 = p2 - p1;
    let b2 = p3 - p2;
    let b3 = p4 - p3;

    let b2_norm = b2.norm();
    if b2_norm < 1e-6 {
        return 0.0;
    }

    let n1 = b1.cross(b2);
    let n2 = b2.cross(b3);

    let m = n1.cross(b2 / b2_norm);

    let x = n1.dot(n2);
    let y = m.dot(n2);

    if x.abs() < 1e-6 && y.abs() < 1e-6 {
        return 0.0;
    }

    y.atan2(x).to_degrees()
}

/// Classifies a pair of backbone dihedral angles $(\phi, \psi)$ into a Ramachandran region.
pub fn classify_ramachandran(phi: f32, psi: f32) -> RamachandranRegion {
    // Normalization to [-180, 180]
    let phi = ((phi + 180.0).rem_euclid(360.0)) - 180.0;
    let psi = ((psi + 180.0).rem_euclid(360.0)) - 180.0;

    // Right-handed alpha helix region
    if (-160.0..=-20.0).contains(&phi) && (-100.0..=30.0).contains(&psi) {
        RamachandranRegion::AlphaHelix
    }
    // Beta sheet / extended region
    else if (-180.0..=-40.0).contains(&phi) && (40.0..=180.0).contains(&psi)
        || (-180.0..=-40.0).contains(&phi) && (-180.0..=-150.0).contains(&psi)
    {
        RamachandranRegion::BetaSheet
    }
    // Left-handed alpha helix
    else if (20.0..=100.0).contains(&phi) && (-20.0..=100.0).contains(&psi) {
        RamachandranRegion::LeftHandedAlpha
    } else {
        RamachandranRegion::Outlier
    }
}

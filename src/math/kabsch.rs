//! Kabsch 3D Structural Superposition & SVD Algorithm.
//!
//! Computes the optimal rigid-body rotation matrix $R$ and translation vector $t$
//! that minimizes the Root-Mean-Square Deviation (RMSD) between two paired 3D point sets.

use crate::math::Vec3;

/// Result of Kabsch structural superposition.
#[derive(Debug, Clone, PartialEq)]
pub struct KabschResult {
    /// Optimal $3 \times 3$ rotation matrix (row-major)
    pub rotation: [[f32; 3]; 3],
    /// Optimal translation vector $t = \bar{q} - R \bar{p}$
    pub translation: Vec3,
    /// Root-mean-square deviation across aligned point pairs after superposition
    pub rmsd: f32,
    /// Number of paired points used for alignment
    pub num_points: usize,
}

impl KabschResult {
    /// Applies the rotation and translation to a 3D point: $p' = R p + t$.
    #[inline]
    pub fn transform_point(&self, p: Vec3) -> Vec3 {
        let r = &self.rotation;
        let x = r[0][0] * p.x + r[0][1] * p.y + r[0][2] * p.z + self.translation.x;
        let y = r[1][0] * p.x + r[1][1] * p.y + r[1][2] * p.z + self.translation.y;
        let z = r[2][0] * p.x + r[2][1] * p.y + r[2][2] * p.z + self.translation.z;
        Vec3::new(x, y, z)
    }
}

/// Jacobi eigenvalue decomposition for a symmetric $3 \times 3$ matrix $A$.
/// Returns eigenvalues and eigenvector columns matrix $V$: $A V = V \Lambda$.
#[allow(clippy::needless_range_loop)]
pub fn jacobi_eigen_3x3(mut a: [[f32; 3]; 3]) -> ([f32; 3], [[f32; 3]; 3]) {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    for _ in 0..50 {
        // Find largest off-diagonal element
        let mut max_off = 0.0f32;
        let mut p = 0;
        let mut q = 1;

        if a[0][1].abs() > max_off {
            max_off = a[0][1].abs();
            p = 0;
            q = 1;
        }
        if a[0][2].abs() > max_off {
            max_off = a[0][2].abs();
            p = 0;
            q = 2;
        }
        if a[1][2].abs() > max_off {
            max_off = a[1][2].abs();
            p = 1;
            q = 2;
        }

        if max_off < 1e-7 {
            break;
        }

        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];

        let theta = (aqq - app) / (2.0 * apq);
        let t = if theta >= 0.0 {
            1.0 / (theta + (1.0 + theta * theta).sqrt())
        } else {
            -1.0 / (-theta + (1.0 + theta * theta).sqrt())
        };

        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        let tau = s / (1.0 + c);

        // Update A in-place
        a[p][p] -= t * apq;
        a[q][q] += t * apq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for r in 0..3 {
            if r != p && r != q {
                let arp = a[r][p];
                let arq = a[r][q];
                a[r][p] = arp - s * (arq + tau * arp);
                a[p][r] = a[r][p];
                a[r][q] = arq + s * (arp - tau * arq);
                a[q][r] = a[r][q];
            }
        }

        // Update V columns
        for r in 0..3 {
            let vrp = v[r][p];
            let vrq = v[r][q];
            v[r][p] = vrp - s * (vrq + tau * vrp);
            v[r][q] = vrq + s * (vrp - tau * vrq);
        }
    }

    let evals = [a[0][0], a[1][1], a[2][2]];
    (evals, v)
}

/// Singular Value Decomposition for a general $3 \times 3$ matrix $A$: $A = U \Sigma V^T$.
#[allow(clippy::needless_range_loop)]
pub fn svd_3x3(a: [[f32; 3]; 3]) -> ([[f32; 3]; 3], [f32; 3], [[f32; 3]; 3]) {
    // 1. ATA = A^T * A
    let mut ata = [[0.0f32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                ata[i][j] += a[k][i] * a[k][j];
            }
        }
    }

    // 2. Eigenvalues and eigenvectors of ATA: ATA * V = V * S^2
    let (evals, v_raw) = jacobi_eigen_3x3(ata);

    // Sort singular values descending
    let mut order = [0, 1, 2];
    order.sort_by(|&i, &j| {
        evals[j]
            .partial_cmp(&evals[i])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut s = [0.0f32; 3];
    let mut v = [[0.0f32; 3]; 3];
    for (new_idx, &old_idx) in order.iter().enumerate() {
        s[new_idx] = evals[old_idx].max(0.0).sqrt();
        for row in 0..3 {
            v[row][new_idx] = v_raw[row][old_idx];
        }
    }

    // 3. Compute columns of U: u_i = A * v_i / s_i
    let mut u = [[0.0f32; 3]; 3];
    for col in 0..3 {
        if s[col] > 1e-5 {
            let inv_s = 1.0 / s[col];
            for row in 0..3 {
                let mut sum = 0.0f32;
                for k in 0..3 {
                    sum += a[row][k] * v[k][col];
                }
                u[row][col] = sum * inv_s;
            }
        }
    }

    // Handle rank deficiency: complete orthogonal basis for U if needed
    if s[1] <= 1e-5 {
        // u0 is valid; compute u1 and u2 perpendicular to u0
        let u0 = Vec3::new(u[0][0], u[1][0], u[2][0]);
        let temp = if u0.x.abs() < 0.8 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let u1 = u0.cross(temp).normalize();
        let u2 = u0.cross(u1).normalize();
        u[0][1] = u1.x;
        u[1][1] = u1.y;
        u[2][1] = u1.z;
        u[0][2] = u2.x;
        u[1][2] = u2.y;
        u[2][2] = u2.z;
    } else if s[2] <= 1e-5 {
        // u0 and u1 are valid; u2 = u0 x u1
        let u0 = Vec3::new(u[0][0], u[1][0], u[2][0]);
        let u1 = Vec3::new(u[0][1], u[1][1], u[2][1]);
        let u2 = u0.cross(u1).normalize();
        u[0][2] = u2.x;
        u[1][2] = u2.y;
        u[2][2] = u2.z;
    }

    (u, s, v)
}

/// Determinant of a $3 \times 3$ matrix.
#[inline]
pub fn det_3x3(m: &[[f32; 3]; 3]) -> f32 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Computes the optimal Kabsch rigid-body alignment superimposing source points $P$ onto target points $Q$.
pub fn kabsch_align(p: &[Vec3], q: &[Vec3]) -> Option<KabschResult> {
    if p.len() < 3 || p.len() != q.len() {
        return None;
    }

    let n = p.len() as f32;

    // 1. Centroids
    let mut p_com = Vec3::ZERO;
    let mut q_com = Vec3::ZERO;
    for (&pi, &qi) in p.iter().zip(q.iter()) {
        p_com += pi;
        q_com += qi;
    }
    p_com /= n;
    q_com /= n;

    // 2. Cross-covariance matrix H = sum((p_i - p_com) * (q_i - q_com)^T)
    // H[i][j] = sum( (p_k - p_com)_i * (q_k - q_com)_j )
    let mut h = [[0.0f32; 3]; 3];
    for (&pi, &qi) in p.iter().zip(q.iter()) {
        let px = pi.x - p_com.x;
        let py = pi.y - p_com.y;
        let pz = pi.z - p_com.z;

        let qx = qi.x - q_com.x;
        let qy = qi.y - q_com.y;
        let qz = qi.z - q_com.z;

        h[0][0] += px * qx;
        h[0][1] += px * qy;
        h[0][2] += px * qz;

        h[1][0] += py * qx;
        h[1][1] += py * qy;
        h[1][2] += py * qz;

        h[2][0] += pz * qx;
        h[2][1] += pz * qy;
        h[2][2] += pz * qz;
    }

    // 3. SVD of H = U Sigma V^T
    let (u, _, v) = svd_3x3(h);

    // 4. Optimal rotation R = V * diag(1, 1, d) * U^T
    // Compute d = sign(det(V * U^T))
    let mut v_ut = [[0.0f32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                v_ut[i][j] += v[i][k] * u[j][k]; // U^T[k][j] = U[j][k]
            }
        }
    }

    let det = det_3x3(&v_ut);
    let d = if det < 0.0 { -1.0 } else { 1.0 };

    let mut r = [[0.0f32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = v[i][0] * u[j][0] + v[i][1] * u[j][1] + v[i][2] * d * u[j][2];
        }
    }

    // 5. Translation vector t = q_com - R * p_com
    let rx = r[0][0] * p_com.x + r[0][1] * p_com.y + r[0][2] * p_com.z;
    let ry = r[1][0] * p_com.x + r[1][1] * p_com.y + r[1][2] * p_com.z;
    let rz = r[2][0] * p_com.x + r[2][1] * p_com.y + r[2][2] * p_com.z;
    let translation = Vec3::new(q_com.x - rx, q_com.y - ry, q_com.z - rz);

    let res = KabschResult {
        rotation: r,
        translation,
        rmsd: 0.0,
        num_points: p.len(),
    };

    // 6. Compute RMSD
    let mut sum_sq = 0.0f32;
    for (&pi, &qi) in p.iter().zip(q.iter()) {
        let transformed = res.transform_point(pi);
        sum_sq += transformed.distance_squared(&qi);
    }
    let rmsd = (sum_sq / n).sqrt();

    Some(KabschResult { rmsd, ..res })
}

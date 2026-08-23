# Scientific Tooling & Analysis Architecture Design Spec

## 1. Overview
This specification details the technical design for **Scientific Tooling & Analysis** in `termpdb`:
1. **Kabsch 3D Structural Superposition & Alignment**: Optimal rigid-body rotation and translation via SVD/eigen-decomposition, with Needleman-Wunsch sequence alignment for automatic $C\alpha$ pairing and per-residue RMSD heatmaps.
2. **Contact Geometry & Ramachandran Analysis**: 1-to-4 atom selection queue computing Euclidean distances, 3-atom bond angles $\theta^\circ$, and 4-atom dihedral angles $\phi/\psi/\omega^\circ$ with Ramachandran quadrant classification.
3. **Pure Rust DSSP Secondary Structure Assignment**: Backbone electrostatic hydrogen-bonding energy calculation ($E = q_1 q_2 [1/r_{ON} + 1/r_{CH} - 1/r_{OH} - 1/r_{CN}] \cdot 332\text{ kcal}\cdot\text{\AA}/\text{mol}$) to assign $\alpha$-helices, $3_{10}$-helices, $\beta$-sheets, and coils when files lack annotations.

---

## 2. Algorithms & Technical Details

### 2.1 Kabsch Superposition & Sequence Alignment

#### Needleman-Wunsch Sequence Alignment
- Maps two amino acid sequences $A$ and $B$ to identify aligned residues:
  - BLOSUM62 substitution matrix + affine gap penalty ($g_{\text{open}} = -10, g_{\text{extend}} = -1$).
  - Produces aligned index pairs $(i_k, j_k)$ where $C\alpha(A_{i_k})$ corresponds to $C\alpha(B_{j_k})$.

#### Kabsch Algorithm (Minimizing RMSD)
Given paired coordinates $P = \{p_1, \dots, p_N\}$ and $Q = \{q_1, \dots, q_N\}$:
1. Compute centroids $\bar{p} = \frac{1}{N}\sum p_i$, $\bar{q} = \frac{1}{N}\sum q_i$.
2. Center coordinate sets: $x_i = p_i - \bar{p}$, $y_i = q_i - \bar{q}$.
3. Form covariance matrix $H = X^T Y = \sum_{i=1}^N x_i y_i^T$ ($3 \times 3$).
4. Compute SVD: $H = U \Sigma V^T$.
5. Optimal rotation matrix: $R = V \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & d \end{pmatrix} U^T$, where $d = \text{sign}(\det(V U^T))$ ensuring a proper rotation ($\det(R) = +1$).
6. Translation vector: $t = \bar{q} - R \bar{p}$.
7. Optimal RMSD: $\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N \| R x_i - y_i \|^2}$.

---

### 2.2 Contact Geometry & Dihedral Angles

#### Planar Bond Angle (3 Atoms: A, B, C)
$$\vec{v}_1 = A - B, \quad \vec{v}_2 = C - B$$
$$\theta = \arccos\left(\frac{\vec{v}_1 \cdot \vec{v}_2}{\|\vec{v}_1\| \|\vec{v}_2\|}\right) \times \frac{180^\circ}{\pi}$$

#### Dihedral / Torsion Angle (4 Atoms: A, B, C, D)
$$\vec{b}_1 = B - A, \quad \vec{b}_2 = C - B, \quad \vec{b}_3 = D - C$$
$$\vec{n}_1 = \vec{b}_1 \times \vec{b}_2, \quad \vec{n}_2 = \vec{b}_2 \times \vec{b}_3$$
$$\vec{m} = \vec{n}_1 \times \frac{\vec{b}_2}{\|\vec{b}_2\|}$$
$$x = \vec{n}_1 \cdot \vec{n}_2, \quad y = \vec{m} \cdot \vec{n}_2$$
$$\phi = \text{atan2}(y, x) \times \frac{180^\circ}{\pi} \in [-180^\circ, 180^\circ]$$

#### Ramachandran Region Classification
Given backbone $\phi$ and $\psi$:
- **Core $\beta$-sheet**: $\phi \in [-180^\circ, -45^\circ], \psi \in [45^\circ, 180^\circ]$
- **Core $\alpha$-helix**: $\phi \in [-120^\circ, -30^\circ], \psi \in [-60^\circ, 0^\circ]$
- **Left-handed $\alpha$-helix**: $\phi \in [30^\circ, 90^\circ], \psi \in [0^\circ, 90^\circ]$
- **Allowed / Outlier**: other regions.

---

### 2.3 Pure Rust DSSP Secondary Structure Assignment

#### Electrostatic Energy Formula
$$E = q_1 q_2 \left( \frac{1}{r_{ON}} + \frac{1}{r_{CH}} - \frac{1}{r_{OH}} - \frac{1}{r_{CN}} \right) \cdot 332.0\text{ kcal}\cdot\text{\AA}/\text{mol}$$
Where $q_1 = 0.42e, q_2 = 0.20e$, $H$ position is computed from $N, C_{\alpha}, C_{i-1}$ geometry if missing:
$$\vec{H} = \vec{N} + \frac{\vec{N} - \vec{C}_{\alpha}}{\|\vec{N} - \vec{C}_{\alpha}\|} + \frac{\vec{N} - \vec{C}_{i-1}}{\|\vec{N} - \vec{C}_{i-1}\|}$$
Normalized to $1.0\text{ \AA}$ $N-H$ bond length.

#### Hydrogen Bond Threshold
Residue $i$ and $j$ form an H-bond if $E(i, j) < -0.5\text{ kcal/mol}$.
- **$\alpha$-helix**: consecutive $i \to i+4$ H-bonds.
- **$3_{10}$-helix**: consecutive $i \to i+3$ H-bonds.
- **$\pi$-helix**: consecutive $i \to i+5$ H-bonds.
- **$\beta$-sheet**: parallel ($i \to j, j \to i+2$) or antiparallel ($i \to j, j \to i$) ladders.

---

## 3. Module Layout
- `src/math/kabsch.rs`: Covariance matrix, SVD for $3 \times 3$, optimal rotation and translation.
- `src/model/align.rs`: Needleman-Wunsch sequence alignment, coordinate pairing, multi-structure alignment.
- `src/model/geometry.rs`: Bond angle, dihedral angle, Ramachandran classification.
- `src/model/dssp.rs`: Electrostatic H-bond energy matrix, helical and sheet pattern detection.
- `src/select.rs`: Extended selection queue supporting up to 4 atoms.
- `src/tui/`: Updated HUD header and status overlays.
- `src/cli.rs`: CLI flags `--align`, `--dssp`, `--angle`, `--dihedral`.

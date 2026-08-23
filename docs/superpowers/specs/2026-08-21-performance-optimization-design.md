# Performance & Optimization Architecture Design Spec

## 1. Overview
This document specifies the technical architecture for the performance milestones in [`ROADMAP.md`](../../ROADMAP.md) Section 2:
1. **Multi-Threaded Tile/Chunk Rasterizer (`rayon`)**: Parallelizing software rendering across CPU cores for interactive 60+ FPS playback on large complexes (>100k atoms).
2. **Fast Vectorized Analytical Sphere & Cylinder Rasterization**: Optimized ray-sphere intersection with frustum bounding box culling, early-Z tests, and branchless depth buffer updates.
3. **Spatial Grid Indexing & Occlusion/Frustum Culling**: Uniform spatial hash grid for instantaneous neighbor queries, buried interior atom culling in VDW mode, and screen-space tile binning.

---

## 2. Architecture & Design

### 2.1 Multi-Threaded Rasterization Architecture
- In CPU software rasterization, writing to a single shared framebuffer concurrently causes data races on `pixels` and `depth`.
- **Approach**: Horizontal tile banding / chunk partitioning.
  - Divide the `Framebuffer` of height $H$ into $N$ horizontal chunks (e.g. bands of 16-32 rows).
  - Each chunk owns its slice of pixels and depth buffer (`&mut [PixelColor]` and `&mut [f32]`).
  - Using `rayon`, dispatch chunk rendering across all available hardware threads (`par_chunks_mut`).
  - Screen bounding box clipping ensures each atom/cylinder is only drawn into the chunks it intersects.

### 2.2 Vectorized Sphere & Ray Analytical Rasterizer
- Sphere rasterization is the #1 CPU consumer in VDW and Ball & Stick modes ($O(\sum \pi r_i^2)$ pixels).
- **Optimizations**:
  - **Screen-Space Frustum Culling**: Project sphere center $(x_c, y_c, z_c)$ and projected radius $R$. If $[x_c - R, x_c + R] \times [y_c - R, y_c + R]$ is disjoint from $[0, W) \times [0, H)$, discard immediately.
  - **Early Z Bounding**: If $z_{\text{front}} = z_c - r_{\text{world}} > \text{min\_depth\_in\_tile}$, skip.
  - **Precomputed Inverse & Squared Tables**: Precompute $R^2$ and $1/R$. The circle condition is $(dx^2 + dy^2) \le R^2$.
  - **Analytical Normal & Half-Vector Fast Path**: Compute local normal $n_z = \sqrt{1 - (dx^2+dy^2)/R^2}$, $n_x = dx/R$, $n_y = dy/R$ with single-precision vectorized reciprocal square root / sqrt.

### 2.3 Spatial Grid Indexing & Interior Atom Occlusion Culling
- In large proteins (e.g. ribosome >100k atoms), >70% of atoms are buried in the core and never visible from the outside in VDW mode.
- **Uniform Spatial Hash Grid (`SpatialGrid`)**:
  - Bin atoms into cubic cells of size $\approx 4.0\text{ \AA}$.
  - Compute burial score: count number of neighbors within $r_{\text{cutoff}} \approx 3.5\text{ \AA}$. Atoms with $\ge 12$ neighbors in 3D are completely enclosed (buried core).
  - In VDW mode, skip buried atoms unless viewing sliced/cavity sections or near clipping plane.

### 2.4 Parallel Postprocessing & Cache Building
- Parallelize `apply_postprocessing` (Sobel depth outlines and SSAO kernels) across row chunks via `rayon`.
- Parallelize `build_render_cache` and color mapping using `par_iter()`.

---

## 3. Interfaces & Data Structures

### `SpatialGrid`
```rust
pub struct SpatialGrid<'a> {
    cell_size: f32,
    inv_cell_size: f32,
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
    atoms: &'a [Atom],
}
```

### Tile Chunk Rasterizer
```rust
pub fn render_spheres_parallel(
    spheres: &[ProjectedSphere],
    buffer: &mut Framebuffer,
    lighting: &Lighting,
);
```

---

## 4. Verification & Testing
- Unit tests for `SpatialGrid` construction and neighbor queries.
- Bounding-box frustum culling accuracy tests.
- Multi-threaded tile chunk rasterizer rendering equivalence tests (parallel output == serial output).
- Performance benchmarks and regression tests verifying 100% deterministic visual parity.

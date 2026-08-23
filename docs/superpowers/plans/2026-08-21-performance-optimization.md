# Performance & Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement high-performance software rendering optimizations (Spatial Grid neighbor indexing & occlusion culling, vectorized sphere/cylinder rasterization with frustum clipping, multi-threaded Rayon tile-chunk rasterizer, and parallel post-processing).

**Architecture:** Add `rayon` dependency; introduce `termpdb::model::spatial` for uniform 3D hash grid indexing and buried-atom culling; optimize `src/render/rasterizer.rs` with analytical early-Z and screen-space frustum bounds; implement multi-threaded horizontal band rendering in `src/render/buffer.rs` and `src/render/postprocess.rs`.

**Tech Stack:** Rust 2024, `rayon` (1.10+), `ratatui` (0.29+), `crossterm` (0.28+).

## Global Constraints
- Pure Terminal ANSI / Half-Block & Unicode Braille (zero GPU dependencies, works over tmux/SSH)
- Deterministic 3D math and analytical rasterization (parallel rendering must produce identical pixels to serial rendering)
- 100% test pass rate across all unit and integration test suites with zero compiler warnings

---

### Task 1: Add Rayon Dependency & Parallel Framebuffer Chunking Scaffolding

**Files:**
- Modify: `Cargo.toml`
- Modify: `src/render/buffer.rs`
- Test: `tests/perf_buffer_test.rs`

- [x] **Step 1: Write unit tests in `tests/perf_buffer_test.rs` testing parallel horizontal row-band slicing and deterministic pixel writes**
- [x] **Step 2: Run `cargo test --test perf_buffer_test` to verify failure**
- [x] **Step 3: Add `rayon = "1.10"` to `Cargo.toml` and implement `Framebuffer::par_bands_mut`**
- [x] **Step 4: Run `cargo test --test perf_buffer_test` to verify pass**

---

### Task 2: Spatial Grid Indexing & Occlusion Culling

**Files:**
- Create: `src/model/spatial.rs`
- Modify: `src/model/mod.rs`
- Modify: `src/model/structure.rs`
- Modify: `src/render/representations/vdw.rs`
- Test: `tests/perf_spatial_test.rs`

- [x] **Step 1: Write unit tests in `tests/perf_spatial_test.rs` for `SpatialGrid` neighbor queries and buried-atom detection**
- [x] **Step 2: Run `cargo test --test perf_spatial_test` to verify failure**
- [x] **Step 3: Implement `SpatialGrid` uniform 3D hash grid in `src/model/spatial.rs`**
- [x] **Step 4: Add core/buried atom filter to VDW render cache to skip non-visible interior atoms**
- [x] **Step 5: Run `cargo test --test perf_spatial_test` to verify pass**

---

### Task 3: Vectorized Sphere & Cylinder Rasterizer with Frustum Bounds

**Files:**
- Modify: `src/render/rasterizer.rs`
- Modify: `src/render/representations/vdw.rs`
- Modify: `src/render/representations/ball_stick.rs`
- Test: `tests/perf_rasterizer_test.rs`

- [x] **Step 1: Write unit tests in `tests/perf_rasterizer_test.rs` for fast sphere rasterization with frustum clipping and early-Z**
- [x] **Step 2: Run `cargo test --test perf_rasterizer_test` to verify failure**
- [x] **Step 3: Optimize `draw_sphere_lit` and `draw_cylinder` in `src/render/rasterizer.rs` with precomputed inverse-radii and branchless depth updates**
- [x] **Step 4: Run `cargo test --test perf_rasterizer_test` to verify pass**

---

### Task 4: Parallel Sphere & Post-Processing Tile Pipeline

**Files:**
- Modify: `src/render/postprocess.rs`
- Modify: `src/render/representations/vdw.rs`
- Modify: `src/render/representations/mod.rs`
- Modify: `src/render/buffer.rs`
- Test: `tests/perf_parallel_render_test.rs`

- [x] **Step 1: Write unit tests in `tests/perf_parallel_render_test.rs` verifying that parallel rendering matches serial rendering bit-for-bit**
- [x] **Step 2: Run `cargo test --test perf_parallel_render_test` to verify failure**
- [x] **Step 3: Implement parallel row-band postprocessing for SSAO and outline edge detection**
- [x] **Step 4: Implement parallel sphere rasterization across horizontal framebuffer bands**
- [x] **Step 5: Run `cargo test --test perf_parallel_render_test` to verify pass**

---

### Task 5: Full Workspace Integration, Benchmarking & Verification

**Files:**
- Modify: `src/tui/app.rs`
- Modify: `ROADMAP.md`
- Test: All tests in workspace

- [x] **Step 1: Run full test suite with `cargo test` across all 18 test suites**
- [x] **Step 2: Run `cargo clippy --all-targets -- -D warnings` to verify clean build**
- [x] **Step 3: Update `ROADMAP.md` to reflect completed performance tasks**
- [x] **Step 4: Provide final verification report and mark goal complete**


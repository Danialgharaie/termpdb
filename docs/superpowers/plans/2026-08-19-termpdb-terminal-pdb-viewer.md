# termpdb Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `termpdb`, a high-performance terminal 3D macromolecular structure viewer in Rust supporting PDB/mmCIF files and RCSB fetching, with interactive arcball controls, 4 structural representations (Trace, Ball & Stick, Ribbon, VDW), 6 color schemes, and truecolor half-block/Braille Z-buffer rendering.

**Architecture:** A modular architecture separating structure parsing (`parser`), chemical data modeling & bond detection (`model`), 3D math & splines (`math`), analytical software rasterization with depth testing and directional lighting (`render`), Ratatui TUI event loop & custom viewport widget (`tui`), and CLI entry point (`cli`).

**Tech Stack:** Rust 2024, `ratatui` (0.29+), `crossterm` (0.28+), `clap` (4.5+), `ureq` (2.12+), `flate2` (1.0+), `anyhow` / `thiserror`.

## Global Constraints
- Target Language: Rust (2024 Edition)
- Terminal Graphics: 24-bit Truecolor ANSI Half-blocks (`▀`) and Unicode Braille Canvas
- Zero GPU / C library dependencies (pure Rust software rasterization)
- Maximum responsiveness: target 60 FPS in terminal during interactive rotation
- Unit test coverage for parsers, 3D math, Z-buffer, and geometry algorithms

---

### Task 1: Cargo Configuration, Dependencies, Error Handling & 3D Math Primitives

**Files:**
- Modify: `Cargo.toml`
- Create: `src/error.rs`
- Create: `src/math/mod.rs`
- Create: `src/math/vec3.rs`
- Create: `src/math/mat4.rs`
- Create: `src/math/quat.rs`
- Create: `src/math/spline.rs`
- Test: `tests/math_test.rs`

**Interfaces:**
- Produces:
  - `termpdb::error::Result<T>`, `termpdb::error::TermPdbError`
  - `termpdb::math::Vec3` (`new`, `dot`, `cross`, `norm`, `normalize`, `lerp`, operators `+`, `-`, `*`)
  - `termpdb::math::Mat4` (`identity`, `look_at`, `perspective`, `orthographic`, `transform_point`, `transform_vector`, `mul`)
  - `termpdb::math::Quat` (`identity`, `from_axis_angle`, `from_euler`, `from_drag`, `rotate_vec3`, `mul`, `to_mat4`)
  - `termpdb::math::CatmullRomSpline` (`interpolate`, `tangent`, `generate_smooth_curve`)

- [ ] **Step 1: Update `Cargo.toml` with project dependencies**

Update `Cargo.toml` with `ratatui`, `crossterm`, `clap`, `ureq`, `flate2`, `anyhow`, and `thiserror`.

- [ ] **Step 2: Write unit tests for 3D math primitives in `tests/math_test.rs`**

Include tests for vector dot/cross product, matrix-vector multiplication, quaternion rotation, and Catmull-Rom spline interpolation.

- [ ] **Step 3: Run `cargo test --test math_test` to verify failure**

Run: `cargo test --test math_test`  
Expected: FAIL (modules not found)

- [ ] **Step 4: Implement `src/error.rs` and `src/math/*` modules**

Implement `Vec3`, `Mat4`, `Quat`, and `CatmullRomSpline` with numerical stability.

- [ ] **Step 5: Run `cargo test --test math_test` to verify pass**

Run: `cargo test --test math_test`  
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml src/error.rs src/math tests/math_test.rs
git commit -m "feat(math): add 3D math vector, matrix, quaternion and spline primitives"
```

---

### Task 2: Chemical Data Model, Elements & Bond Graph

**Files:**
- Create: `src/model/elements.rs`
- Create: `src/model/atom.rs`
- Create: `src/model/residue.rs`
- Create: `src/model/chain.rs`
- Create: `src/model/bond.rs`
- Create: `src/model/mod.rs`
- Test: `tests/model_test.rs`

**Interfaces:**
- Consumes: `termpdb::math::Vec3`
- Produces:
  - `termpdb::model::Element` (atomic number, symbol, covalent radius, VDW radius, default CPK color)
  - `termpdb::model::SecondaryStructure` (`Helix`, `Sheet`, `Coil`)
  - `termpdb::model::Atom` (index, serial, name, element, pos, b_factor, res_name, res_seq, chain_id, is_hetatm)
  - `termpdb::model::Residue` (seq, name, chain_id, atoms, secondary_structure, is_amino_acid, is_nucleic)
  - `termpdb::model::Chain` (id, residues, atoms)
  - `termpdb::model::Structure` (title, chains, atoms, bonds, center_of_mass, bounding_sphere_radius, center_and_normalize)
  - `termpdb::model::BondDetector` (spatial hash grid bond detection)

- [ ] **Step 1: Write unit tests for data model and bond detection in `tests/model_test.rs`**

Test element lookup, structure centering/bounding sphere computation, residue extraction, and covalent bond detection.

- [ ] **Step 2: Run `cargo test --test model_test` to verify failure**

Run: `cargo test --test model_test`  
Expected: FAIL

- [ ] **Step 3: Implement `src/model/*` modules**

Implement `Element` periodic table database, `Atom`, `Residue`, `Chain`, `Structure`, and spatial grid `BondDetector`.

- [ ] **Step 4: Run `cargo test --test model_test` to verify pass**

Run: `cargo test --test model_test`  
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/model tests/model_test.rs
git commit -m "feat(model): add molecular data model, element table and bond detector"
```

---

### Task 3: PDB & mmCIF Parsers and RCSB Fetcher

**Files:**
- Create: `src/parser/pdb.rs`
- Create: `src/parser/cif.rs`
- Create: `src/parser/rcsb.rs`
- Create: `src/parser/mod.rs`
- Test: `tests/parser_test.rs`

**Interfaces:**
- Consumes: `termpdb::model::Structure`, `termpdb::error::Result`
- Produces:
  - `termpdb::parser::parse_pdb(input: &str) -> Result<Structure>`
  - `termpdb::parser::parse_cif(input: &str) -> Result<Structure>`
  - `termpdb::parser::load_structure(source: &str) -> Result<Structure>` (supports local path, `.gz`, or 4-letter RCSB ID)
  - `termpdb::parser::rcsb::fetch_pdb(pdb_id: &str) -> Result<String>`

- [ ] **Step 1: Write parser unit tests with embedded PDB/mmCIF snippets in `tests/parser_test.rs`**

Test ATOM, HETATM, HELIX, SHEET, CONECT parsing in PDB and `_atom_site` loop parsing in mmCIF.

- [ ] **Step 2: Run `cargo test --test parser_test` to verify failure**

Run: `cargo test --test parser_test`  
Expected: FAIL

- [ ] **Step 3: Implement `src/parser/pdb.rs`, `src/parser/cif.rs`, `src/parser/rcsb.rs`, and `src/parser/mod.rs`**

Implement robust line-based PDB parser, mmCIF loop parser, transparent `.gz` decompression, and RCSB HTTPS fetcher.

- [ ] **Step 4: Run `cargo test --test parser_test` to verify pass**

Run: `cargo test --test parser_test`  
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/parser tests/parser_test.rs
git commit -m "feat(parser): add PDB, mmCIF, gz decompression and RCSB fetcher"
```

---

### Task 4: Framebuffer, Dual-Pixel Z-Buffer, Lighting & Color Schemes

**Files:**
- Create: `src/render/buffer.rs`
- Create: `src/render/lighting.rs`
- Create: `src/render/color.rs`
- Create: `src/render/camera.rs`
- Create: `src/render/rasterizer.rs`
- Create: `src/render/mod.rs`
- Test: `tests/render_test.rs`

**Interfaces:**
- Consumes: `termpdb::math::*`, `termpdb::model::*`
- Produces:
  - `termpdb::render::Framebuffer` (`new(width, height)`, `clear`, `set_pixel`, `get_half_blocks`)
  - `termpdb::render::Camera` (`new`, `orbit`, `pan`, `zoom`, `view_matrix`, `proj_matrix`, `world_to_screen`)
  - `termpdb::render::ColorScheme` (`CPK`, `Rainbow`, `Chain`, `SecondaryStructure`, `BFactor`, `Hydrophobicity`)
  - `termpdb::render::Rasterizer` (`draw_sphere`, `draw_line_3d`, `draw_cylinder`, `draw_triangle_3d`)
  - `termpdb::render::Lighting` (`compute_shade(normal, depth, base_color)`)

- [ ] **Step 1: Write render buffer and Z-buffer tests in `tests/render_test.rs`**

Test pixel writing, depth occlusion (near pixel overrides far pixel, far pixel ignored), camera transformations, and color scheme mapping.

- [ ] **Step 2: Run `cargo test --test render_test` to verify failure**

Run: `cargo test --test render_test`  
Expected: FAIL

- [ ] **Step 3: Implement `src/render/buffer.rs`, `lighting.rs`, `color.rs`, `camera.rs`, `rasterizer.rs`, and `mod.rs`**

Implement truecolor half-block framebuffer, depth buffer, Lambertian shader with depth fog, camera orbit/arcball projection, analytical sphere rasterization, and 3D Bresenham line drawing.

- [ ] **Step 4: Run `cargo test --test render_test` to verify pass**

Run: `cargo test --test render_test`  
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/render tests/render_test.rs
git commit -m "feat(render): add truecolor framebuffer, Z-buffer, lighting and rasterizer"
```

---

### Task 5: 3D Structural Representations (Trace, Ball & Stick, Ribbon, VDW)

**Files:**
- Create: `src/render/representations/trace.rs`
- Create: `src/render/representations/ball_stick.rs`
- Create: `src/render/representations/ribbon.rs`
- Create: `src/render/representations/vdw.rs`
- Create: `src/render/representations/mod.rs`
- Test: `tests/representations_test.rs`

**Interfaces:**
- Consumes: `termpdb::render::*`, `termpdb::model::Structure`
- Produces:
  - `termpdb::render::RenderMode` (`Trace`, `BallAndStick`, `Ribbon`, `Vdw`)
  - `termpdb::render::representations::render_structure(structure, mode, color_scheme, camera, buffer)`

- [ ] **Step 1: Write representation rendering tests in `tests/representations_test.rs`**

Test that each representation (Trace, Ball & Stick, Ribbon, VDW) renders to the framebuffer without errors or panics.

- [ ] **Step 2: Run `cargo test --test representations_test` to verify failure**

Run: `cargo test --test representations_test`  
Expected: FAIL

- [ ] **Step 3: Implement representation generators in `src/render/representations/*`**

Implement Backbone Trace, Ball & Stick, Secondary Structure Cartoon Ribbon with Catmull-Rom spline, and Van der Waals space-filling sphere generator.

- [ ] **Step 4: Run `cargo test --test representations_test` to verify pass**

Run: `cargo test --test representations_test`  
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/render/representations tests/representations_test.rs
git commit -m "feat(render): implement trace, ball-and-stick, ribbon and VDW representations"
```

---

### Task 6: Ratatui TUI Application State, Widgets & Event Loop

**Files:**
- Create: `src/tui/app.rs`
- Create: `src/tui/events.rs`
- Create: `src/tui/widgets/viewport.rs`
- Create: `src/tui/widgets/hud.rs`
- Create: `src/tui/widgets/help.rs`
- Create: `src/tui/widgets/info.rs`
- Create: `src/tui/widgets/mod.rs`
- Create: `src/tui/mod.rs`
- Test: `tests/tui_test.rs`

**Interfaces:**
- Consumes: `termpdb::render::*`, `termpdb::model::Structure`
- Produces:
  - `termpdb::tui::App` (maintains structure, camera, rendering mode, color scheme, auto-spin state, overlay toggles, FPS counter)
  - `termpdb::tui::run(structure, initial_mode, initial_color, spin) -> Result<()>`
  - `termpdb::tui::widgets::ViewportWidget` (Ratatui widget for rendering the framebuffer)

- [ ] **Step 1: Write TUI state and event handling tests in `tests/tui_test.rs`**

Test app state transitions (mode switching, color scheme cycling, spin toggling, camera rotation actions, key/mouse event mapping).

- [ ] **Step 2: Run `cargo test --test tui_test` to verify failure**

Run: `cargo test --test tui_test`  
Expected: FAIL

- [ ] **Step 3: Implement `src/tui/*` modules and Ratatui widgets**

Implement full TUI application state machine, mouse drag/wheel event handling, keyboard shortcuts, HUD bars, help modal, info popup, and 60 FPS tick loop.

- [ ] **Step 4: Run `cargo test --test tui_test` to verify pass**

Run: `cargo test --test tui_test`  
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/tui tests/tui_test.rs
git commit -m "feat(tui): add interactive Ratatui TUI, viewport widget, event loop and HUD"
```

---

### Task 7: CLI Interface, Headless ANSI Export, Integration & End-to-End Tests

**Files:**
- Create: `src/cli.rs`
- Modify: `src/main.rs`
- Create: `src/lib.rs`
- Create: `tests/integration_test.rs`

**Interfaces:**
- Produces:
  - `termpdb` binary with CLI flags (`--mode`, `--color`, `--spin`, `--export-ansi`, etc.)
  - Headless ANSI snapshot export to stdout/file for scripting

- [ ] **Step 1: Implement `src/lib.rs`, `src/cli.rs`, and wire `src/main.rs`**

Expose public library API in `src/lib.rs`, define Clap CLI in `src/cli.rs`, and wire `main.rs` to parse arguments, load structure (or display spinner), and launch TUI or export ANSI.

- [ ] **Step 2: Write end-to-end integration test in `tests/integration_test.rs`**

Test loading sample structures, rendering across all 4 modes, and exporting headless ANSI strings.

- [ ] **Step 3: Run `cargo test` on entire workspace**

Run: `cargo test`  
Expected: All unit and integration tests PASS.

- [ ] **Step 4: Test CLI binary build and execution**

Run: `cargo run -- --help` and verify `--help` output formatting.

- [ ] **Step 5: Commit**

```bash
git add src/lib.rs src/cli.rs src/main.rs tests/integration_test.rs
git commit -m "feat(cli): complete CLI integration, headless ANSI export and E2E tests"
```

---

## Plan Self-Review Checklist
1. **Spec coverage**: Covers all requirements (PDB/mmCIF/GZ/RCSB parsing, 4 rendering modes, 6 color schemes, mouse/keyboard controls, HUD/help, headless export).
2. **No Placeholders**: Every task has concrete filenames, explicit interfaces, step-by-step TDD workflows, and exact commands.
3. **Type consistency**: Vector3, Matrix4, Quaternion, Structure, Atom, Residue, Framebuffer, Camera, App names match across all tasks.

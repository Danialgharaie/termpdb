# Terminal 3D Visual Enhancements Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement 8 pure-terminal visual enhancements for `termpdb` (Blinn-Phong specular highlights, silhouette depth outlines, screen-space ambient occlusion, Braille $2\times 4$ subpixel canvas, H-bond/disulfide interaction meshes, nucleic acid DNA/RNA ladders, AlphaFold pLDDT & curated palettes, and depth-of-field focal cueing).

**Architecture:** Extend `termpdb::render::lighting` with specular reflection and focal depth; introduce depth-buffer postprocessing filters (`postprocess.rs`) for SSAO and silhouette edge detection; implement `braille.rs` for subpixel Unicode rasterization; add interaction and nucleic ribbon generators in `src/model/` and `src/render/representations/`; expand `ColorScheme` and CLI/TUI controls.

**Tech Stack:** Rust 2024, `ratatui` (0.29+), `crossterm` (0.28+), `clap` (4.5+).

## Global Constraints
- Target Language: Rust (2024 Edition)
- Pure Terminal ANSI / Half-Block & Unicode Braille (zero GPU dependencies, works over tmux/SSH)
- Deterministic 3D math and analytical rasterization
- 100% test pass rate across all unit and integration test suites

---

### Task 1: Blinn-Phong Specular Lighting & Material Highlights

**Files:**
- Modify: `src/render/lighting.rs`
- Modify: `src/render/rasterizer.rs`
- Modify: `src/render/representations/ribbon.rs`
- Test: `tests/visual_lighting_test.rs`

- [x] **Step 1: Write unit tests in `tests/visual_lighting_test.rs` for specular highlights and half-vector shading**
- [x] **Step 2: Run `cargo test --test visual_lighting_test` to verify failure**
- [x] **Step 3: Update `Lighting` struct with specular intensity, shininess, and view vector support**
- [x] **Step 4: Update sphere and ribbon shaders to compute specular glint**
- [x] **Step 5: Run `cargo test --test visual_lighting_test` to verify pass**

---

### Task 2: Depth Outlines (Silhouette Edge Detection) & SSAO

**Files:**
- Create: `src/render/postprocess.rs`
- Modify: `src/render/buffer.rs`
- Modify: `src/render/mod.rs`
- Test: `tests/visual_postprocess_test.rs`

- [x] **Step 1: Write unit tests in `tests/visual_postprocess_test.rs` for silhouette edge detection and SSAO crevice darkening**
- [x] **Step 2: Run `cargo test --test visual_postprocess_test` to verify failure**
- [x] **Step 3: Implement `src/render/postprocess.rs` with edge Sobel and SSAO kernels**
- [x] **Step 4: Wire `Framebuffer::apply_postprocessing` into rendering pipeline**
- [x] **Step 5: Run `cargo test --test visual_postprocess_test` to verify pass**

---

### Task 3: AlphaFold pLDDT, Electrostatic & Curated Palettes

**Files:**
- Modify: `src/render/color.rs`
- Modify: `src/cli.rs`
- Test: `tests/visual_color_test.rs`

- [x] **Step 1: Write unit tests in `tests/visual_color_test.rs` for pLDDT thresholds, electrostatic gradient, Catppuccin, Nord, Tokyo Night, Gruvbox**
- [x] **Step 2: Run `cargo test --test visual_color_test` to verify failure**
- [x] **Step 3: Add new `ColorScheme` variants and color map lookup tables in `src/render/color.rs`**
- [x] **Step 4: Update CLI value enums and aliases in `src/cli.rs`**
- [x] **Step 5: Run `cargo test --test visual_color_test` to verify pass**

---

### Task 4: Non-Covalent Interactions & Disulfide Bridge Detector

**Files:**
- Create: `src/model/interactions.rs`
- Create: `src/render/representations/interactions.rs`
- Modify: `src/model/mod.rs`
- Modify: `src/render/representations/mod.rs`
- Test: `tests/visual_interactions_test.rs`

- [x] **Step 1: Write unit tests in `tests/visual_interactions_test.rs` for H-bond donor/acceptor pairing and disulfide `-S-S-` detection**
- [x] **Step 2: Run `cargo test --test visual_interactions_test` to verify failure**
- [x] **Step 3: Implement `InteractionDetector` and dashed 3D line segment generator**
- [x] **Step 4: Wire interactions into representation renderers**
- [x] **Step 5: Run `cargo test --test visual_interactions_test` to verify pass**

---

### Task 5: Nucleic Acid Double-Helix Ribbons & Base-Pair Slabs

**Files:**
- Create: `src/render/representations/nucleic.rs`
- Modify: `src/render/representations/ribbon.rs`
- Modify: `src/render/representations/mod.rs`
- Test: `tests/visual_nucleic_test.rs`

- [x] **Step 1: Write unit tests in `tests/visual_nucleic_test.rs` for DNA/RNA backbone spline and base slab generation on nucleic structures**
- [x] **Step 2: Run `cargo test --test visual_nucleic_test` to verify failure**
- [x] **Step 3: Implement nucleic acid ribbon & base rung generator**
- [x] **Step 4: Run `cargo test --test visual_nucleic_test` to verify pass**

---

### Task 6: Braille Subpixel Canvas ($2 \times 4$ Dot Canvas) & Wireframe Mode

**Files:**
- Create: `src/render/braille.rs`
- Modify: `src/render/representations/mod.rs`
- Modify: `src/render/mod.rs`
- Modify: `src/cli.rs`
- Test: `tests/visual_braille_test.rs`

- [x] **Step 1: Write unit tests in `tests/visual_braille_test.rs` for Braille bitmask generation and $2 \times 4$ line/dot plotting**
- [x] **Step 2: Run `cargo test --test visual_braille_test` to verify failure**
- [x] **Step 3: Implement `BrailleBuffer` and wireframe rendering mode**
- [x] **Step 4: Run `cargo test --test visual_braille_test` to verify pass**

---

### Task 7: Depth-of-Field (DoF) & Focal Plane Cueing

**Files:**
- Modify: `src/render/lighting.rs`
- Modify: `src/render/camera.rs`
- Modify: `src/render/buffer.rs`
- Test: `tests/visual_dof_test.rs`

- [x] **Step 1: Write unit tests in `tests/visual_dof_test.rs` for focal distance attenuation and fog falloff**
- [x] **Step 2: Run `cargo test --test visual_dof_test` to verify failure**
- [x] **Step 3: Implement focal plane depth cueing in `Lighting` and `Camera`**
- [x] **Step 4: Run `cargo test --test visual_dof_test` to verify pass**

---

### Task 8: TUI Keybindings, HUD Indicators, CLI Flags & Full Suite Verification

**Files:**
- Modify: `src/tui/app.rs`
- Modify: `src/tui/events.rs`
- Modify: `src/tui/widgets/hud.rs`
- Modify: `src/tui/widgets/help.rs`
- Modify: `src/main.rs`
- Modify: `src/lib.rs`
- Test: `tests/visual_integration_test.rs`

- [x] **Step 1: Write integration tests in `tests/visual_integration_test.rs` covering all new visual modes, overlays, and CLI options**
- [x] **Step 2: Run `cargo test --test visual_integration_test` to verify failure**
- [x] **Step 3: Wire TUI keybindings (`o` for outlines, `a` for SSAO, `n` for interactions, `d` for DoF, `5` for Wireframe/Braille, new color cycling)**
- [x] **Step 4: Update HUD status lines and help dialog**
- [x] **Step 5: Run all test suites across the workspace with `cargo test` and verify 100% pass**

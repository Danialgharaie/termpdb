# Scientific Tooling & Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement comprehensive scientific analysis tools for `termpdb` (Kabsch 3D structural superposition & alignment with Needleman-Wunsch sequence matching, contact geometry with angles/dihedrals/Ramachandran classification, and pure Rust DSSP secondary structure assignment).

**Architecture:** Add `src/math/kabsch.rs` (SVD & optimal rotation), `src/model/align.rs` (sequence alignment & RMSD heatmap), `src/model/geometry.rs` (angles & dihedrals), `src/model/dssp.rs` (electrostatic H-bond DSSP engine); expand `src/select.rs` to 4-atom queue; update CLI and TUI HUD status widgets.

**Tech Stack:** Rust 2024, `ratatui` (0.29+), `crossterm` (0.28+), `clap` (4.5+).

## Global Constraints
- Target Language: Rust (2024 Edition)
- Zero external C/Fortran libraries (pure software Rust implementation)
- 100% test pass rate across all unit and integration test suites with zero warnings

---

### Task 1: Mathematical Foundations & Kabsch SVD Algorithm

**Files:**
- Create: `src/math/kabsch.rs`
- Modify: `src/math/mod.rs`
- Test: `tests/scientific_kabsch_test.rs`

- [x] **Step 1: Write unit tests in `tests/scientific_kabsch_test.rs` for Kabsch rotation, translation, and identity/known RMSD alignment**
- [x] **Step 2: Run `cargo test --test scientific_kabsch_test` to verify failure**
- [x] **Step 3: Implement $3 \times 3$ SVD and `kabsch_align` in `src/math/kabsch.rs`**
- [x] **Step 4: Run `cargo test --test scientific_kabsch_test` to verify pass**

---

### Task 2: Sequence Alignment (Needleman-Wunsch) & Structural Superposition

**Files:**
- Create: `src/model/align.rs`
- Modify: `src/model/mod.rs`
- Test: `tests/scientific_align_test.rs`

- [x] **Step 1: Write unit tests in `tests/scientific_align_test.rs` for Needleman-Wunsch sequence alignment, coordinate pairing, and multi-structure superposition**
- [x] **Step 2: Run `cargo test --test scientific_align_test` to verify failure**
- [x] **Step 3: Implement BLOSUM62 alignment and `superimpose_structures` in `src/model/align.rs`**
- [x] **Step 4: Run `cargo test --test scientific_align_test` to verify pass**

---

### Task 3: Contact Geometry (Bond Angles, Dihedrals & Ramachandran Classification)

**Files:**
- Create: `src/model/geometry.rs`
- Modify: `src/model/mod.rs`
- Modify: `src/select.rs`
- Test: `tests/scientific_geometry_test.rs`

- [x] **Step 1: Write unit tests in `tests/scientific_geometry_test.rs` for 3-atom bond angle, 4-atom dihedral angle, and Ramachandran region classification**
- [x] **Step 2: Run `cargo test --test scientific_geometry_test` to verify failure**
- [x] **Step 3: Implement angle and dihedral functions in `src/model/geometry.rs`**
- [x] **Step 4: Expand `Selection` in `src/select.rs` to 4-atom FIFO queue with geometry report formatter**
- [x] **Step 5: Run `cargo test --test scientific_geometry_test` to verify pass**

---

### Task 4: Pure Rust DSSP Secondary Structure Assignment

**Files:**
- Create: `src/model/dssp.rs`
- Modify: `src/model/mod.rs`
- Modify: `src/model/structure.rs`
- Test: `tests/scientific_dssp_test.rs`

- [x] **Step 1: Write unit tests in `tests/scientific_dssp_test.rs` for electrostatic H-bond energy matrix and helical/sheet pattern recognition on Crambin (1CRN)**
- [x] **Step 2: Run `cargo test --test scientific_dssp_test` to verify failure**
- [x] **Step 3: Implement DSSP algorithm in `src/model/dssp.rs` and wire automatic fallback into structure parser**
- [x] **Step 4: Run `cargo test --test scientific_dssp_test` to verify pass**

---

### Task 5: TUI HUD Integration, CLI Multi-Structure Flags & Workspace Verification

**Files:**
- Modify: `src/cli.rs`
- Modify: `src/main.rs`
- Modify: `src/tui/app.rs`
- Modify: `src/tui/widgets/hud.rs`
- Modify: `src/tui/widgets/help.rs`
- Modify: `ROADMAP.md`
- Test: `tests/scientific_integration_test.rs`

- [x] **Step 1: Write integration tests in `tests/scientific_integration_test.rs` for CLI multi-structure superposition, angle/dihedral CLI reporting, and DSSP assignment**
- [x] **Step 2: Run `cargo test --test scientific_integration_test` to verify failure**
- [x] **Step 3: Wire CLI arguments (`--align`, `--angle`, `--dihedral`, `--dssp`) and multi-structure support into `src/cli.rs` and `src/main.rs`**
- [x] **Step 4: Update TUI HUD header and selection measurement badges**
- [x] **Step 5: Run full test suite with `cargo test` and `cargo clippy --all-targets -- -D warnings`**
- [x] **Step 6: Update `ROADMAP.md` and complete goal**

# termpdb: 3D Terminal PDB & mmCIF Structure Viewer Design Specification

**Date:** 2026-08-19  
**Status:** Approved  
**Language/Framework:** Rust (2024 Edition), `ratatui`, `crossterm`

---

## 1. Overview & Objectives
`termpdb` is a fast, lightweight, and zero-GPU dependency 3D molecular structure viewer for the terminal. It renders macromolecular structures (proteins, nucleic acids, ligands) loaded from local PDB/mmCIF files or fetched directly from RCSB PDB using high-performance software rasterization with truecolor ANSI half-blocks (`▀`) and Braille subpixel canvases.

### Key Goals:
- **Interactive 3D Navigation**: Arcball mouse rotation, scroll wheel zoom, keyboard controls, and smooth auto-spin.
- **Multiple Representations**: C-alpha/Backbone Trace, Ball & Stick, Secondary Structure Cartoon/Ribbon (Catmull-Rom spline), and Space-filling (Van der Waals) spheres.
- **Multiple Color Schemes**: CPK/Element, Rainbow (N-to-C), Chain ID, Secondary Structure (Helix/Sheet/Loop), B-factor/pLDDT, and Hydrophobicity.
- **Multi-Format Support**: Standard `.pdb`, modern `.cif` (mmCIF), `.gz` compressed files, and direct 4-letter RCSB PDB ID fetching.
- **High Performance**: Pure Rust 60+ FPS analytical software rasterizer with depth Z-buffering and Lambertian directional lighting.

---

## 2. System Architecture

The project is structured into modular Rust crates/submodules under `src/`:

```
src/
├── main.rs                 # CLI entry point, argument parsing, TUI runner
├── cli.rs                  # Clap command line definition
├── parser/
│   ├── mod.rs              # Unified structure loader (.pdb, .cif, .gz, RCSB fetch)
│   ├── pdb.rs              # PDB format parser (ATOM, HETATM, CONECT, HELIX, SHEET)
│   ├── cif.rs              # mmCIF parser (_atom_site, _struct_conf, _struct_sheet_range)
│   └── rcsb.rs             # RCSB REST API client (ureq)
├── model/
│   ├── mod.rs              # Structure data model
│   ├── atom.rs             # Atom struct (coordinates, element, residue, b-factor, radius)
│   ├── residue.rs          # Residue & secondary structure enum
│   ├── chain.rs            # Chain container
│   ├── bond.rs             # Covalent bond graph & spatial grid detection
│   └── elements.rs         # Periodic table metadata (VDW radii, covalent radii, CPK colors)
├── math/
│   ├── mod.rs              # 3D Math primitives
│   ├── vec3.rs             # Vector3 (dot, cross, norm, lerp)
│   ├── mat4.rs             # Matrix4 (perspective, orthographic, look_at, transform)
│   ├── quat.rs             # Quaternion (arcball rotation, slerp)
│   └── spline.rs           # Catmull-Rom spline interpolation for ribbons
├── render/
│   ├── mod.rs              # High-level rendering orchestrator
│   ├── camera.rs           # Orbit/Arcball camera model & viewport projection
│   ├── buffer.rs           # Dual-pixel Half-Block Framebuffer & float Z-Buffer
│   ├── rasterizer.rs       # Analytical sphere rasterizer, 3D line rasterizer, cylinder/ribbon
│   ├── lighting.rs         # Directional lighting, Lambertian diffuse & depth cueing (fog)
│   ├── color.rs            # Color palettes and mapping algorithms
│   └── representations/
│       ├── mod.rs
│       ├── trace.rs        # Backbone / CA trace generator
│       ├── ball_stick.rs   # Ball and stick generator
│       ├── ribbon.rs       # Secondary structure cartoon ribbon generator
│       └── vdw.rs          # Space-filling VDW sphere generator
└── tui/
    ├── mod.rs              # Terminal event loop & Crossterm lifecycle
    ├── app.rs              # Application state machine (mode, camera, color, spin, HUD)
    ├── widgets/
    │   ├── viewport.rs     # Ratatui custom widget rendering the framebuffer
    │   ├── hud.rs          # Header & footer status bars
    │   ├── help.rs         # Keybinding cheat sheet popup
    │   └── info.rs         # Structure metadata & chain breakdown popup
    └── events.rs           # Mouse & Keyboard event handler
```

---

## 3. Detailed Component Specifications

### 3.1 Structure Loader & Parser (`src/parser/`)
- **Format Detection**:
  - Checks file extension (`.pdb`, `.cif`, `.ent`, `.gz`).
  - If input is a 4-character alphanumeric string (e.g. `1CRN`, `7V67`), downloads `https://files.rcsb.org/download/{ID}.cif.gz` or `.pdb.gz`.
  - Automatic transparent decompression via `flate2::read::GzDecoder`.
- **PDB Parser**:
  - Parses `ATOM` / `HETATM` columns for Atom Serial, Name, AltLoc, ResName, ChainID, ResSeq, X, Y, Z, Occupancy, TempFactor (B-factor), Element.
  - Parses `HELIX` and `SHEET` records for secondary structure ranges.
  - Parses `CONECT` records for explicit bonding when present.
- **mmCIF Parser**:
  - Parses `loop_` blocks for `_atom_site.*` fields.
  - Extracts secondary structure from `_struct_conf` and `_struct_sheet_range`.

### 3.2 Molecular Data Model (`src/model/`)
- **Center of Mass & Bounding Box**:
  - Automatically calculates geometric centroid and translates model so center of mass is at $(0, 0, 0)$.
  - Calculates bounding sphere radius $R_{\text{max}}$ for automatic camera framing.
- **Bond Detection**:
  - Spatial hash grid (cell size ~ $3.0\text{ \AA}$) detecting covalent bonds where $d(A, B) \le r_{\text{cov}}(A) + r_{\text{cov}}(B) + 0.45\text{ \AA}$.

### 3.3 3D Math & Rendering Engine (`src/math/` & `src/render/`)
- **Camera**:
  - Orbit camera positioned at distance $D = 2.2 \times R_{\text{max}}$.
  - Orientation stored as a quaternion $\mathbf{q}$ updated via mouse drag delta or keyboard rotation.
  - View-projection matrix transforms world coordinates $(x,y,z)$ into screen-space pixel coordinates $(px, py)$ and depth $z_{\text{cam}}$.
- **Framebuffer & Z-Buffer (`src/render/buffer.rs`)**:
  - Framebuffer width $W = \text{terminal cols}$, height $2H = \text{terminal rows} \times 2$.
  - Depth buffer stores minimum $z$ per subpixel (depth testing).
  - Each terminal character cell renders top pixel and bottom pixel using `▀` (Upper Half Block) with FG RGB and BG RGB.
- **Lighting & Shading**:
  - Directional light vector $\vec{L} = \text{normalize}(0.3, 0.6, 1.0)$.
  - Diffuse intensity: $k_d = \max(0.15, \vec{N} \cdot \vec{L})$.
  - Depth cueing factor: $\text{fog}(z) = \text{clamp}(1.0 - (z - z_{\min}) / (z_{\max} - z_{\min}) \times 0.5, 0.5, 1.0)$.
  - Final color: $C_{\text{final}} = C_{\text{base}} \times k_d \times \text{fog}(z)$.
- **Representations**:
  1. **Trace (1)**: 3D Bresenham line rasterization connecting sequential C-alpha atoms in each chain.
  2. **Ball & Stick (2)**: Analytical sphere rasterization for atoms ($R \approx 0.35\text{ \AA}$) and 3D cylinder segments for bonds.
  3. **Cartoon / Ribbon (3)**: Catmull-Rom spline interpolation along C-alpha sequence:
     - $\alpha$-Helix: Cylindrical ribbon spiral with radius $\approx 1.2\text{ \AA}$.
     - $\beta$-Sheet: Flat planar ribbon with arrowhead at the C-terminal end.
     - Coil: Smooth thin tube ($R \approx 0.3\text{ \AA}$).
  4. **VDW Spheres (4)**: Analytical sphere rasterization with full Van der Waals radii.

### 3.4 Color Schemes (`src/render/color.rs`)
- **CPK**: C (Gray/Green), N (Blue), O (Red), S (Yellow), P (Orange), H (White), Metals (Violet).
- **Rainbow**: Linear HSV color map interpolated from Residue Index $0$ (Blue, $240^\circ$) to $N_{\text{res}}$ (Red, $0^\circ$).
- **Chain ID**: 12 distinct categorical palette colors mapped by Chain ID.
- **Secondary Structure**: Helix (Magenta/Red `#E02060`), Sheet (Gold/Yellow `#E0C020`), Coil (Cyan/Slate `#30A0B0`).
- **B-Factor / pLDDT**: Blue ($\ge 90$), Light Blue ($70-90$), Yellow ($50-70$), Orange/Red ($< 50$).
- **Hydrophobicity**: Hydrophobic (Orange `#E67E22`), Neutral (Gray `#BDC3C7`), Hydrophilic (Sky Blue `#3498DB`).

### 3.5 TUI & Interactive Controls (`src/tui/`)
- **Keybindings**:
  - `1`, `2`, `3`, `4`: Select representation mode (Trace, Ball & Stick, Ribbon, VDW).
  - `c` / `C`: Cycle color scheme forward / backward.
  - `Space`: Toggle continuous auto-spin.
  - `[` / `]`: Adjust auto-spin speed.
  - `h` / `j` / `k` / `l` or Arrow keys: Manual rotation.
  - `+` / `-` or `w` / `s`: Zoom in / out.
  - `H` / `J` / `K` / `L`: Pan camera.
  - `r`: Reset camera and zoom.
  - `i`: Toggle info overlay.
  - `?`: Toggle help modal.
  - `q` / `Esc`: Exit application.
- **Mouse Controls**:
  - Left Drag: Arcball 3D rotation.
  - Scroll Wheel: Zoom.
  - Right Drag / Shift+Drag: Pan.

---

## 4. Error Handling & Edge Cases
- **Network Failures**: When RCSB fetch fails, display clear error message with status code and guidance.
- **Missing or Non-Standard Atoms**: Ignore unrecognized atom types gracefully; parse coordinates even with non-standard formatting.
- **Window Resizing**: Dynamic buffer re-allocation on `crossterm::event::Event::Resize(w, h)`.
- **Large Macromolecules**: For files with $>50,000$ atoms, optimize representation generation to avoid frame drops.

---

## 5. Dependencies
- `ratatui` (TUI framework)
- `crossterm` (Terminal control, raw mode, mouse & keyboard events)
- `clap` (CLI argument parsing with `derive`)
- `ureq` (Lightweight synchronous HTTP client for RCSB fetch)
- `flate2` (Gzip decompression)
- `anyhow` / `thiserror` (Error handling)

---

## 6. Verification & Testing Plan
1. **Unit Tests**:
   - PDB parser parsing test with multi-chain, secondary structure, and CONECT records.
   - mmCIF parser test with loop data and coordinate extraction.
   - Vector3, Matrix4, and Quaternion transformation correctness.
   - Catmull-Rom spline interpolation continuity and smoothness.
   - Z-Buffer occlusion test (near pixel occludes far pixel).
2. **Integration Tests**:
   - Fetch and render standard benchmark structures:
     - `1CRN` (Crambin, small protein)
     - `1UBQ` (Ubiquitin, mixed alpha/beta)
     - `1BNA` (B-DNA dodecamer)
     - `1GFL` (GFP beta-barrel)
3. **Interactive Verification**:
   - Verify 60 FPS rotation and responsiveness in terminal emulator.
   - Verify all 4 rendering modes and 6 color schemes switch cleanly without visual artifacts.

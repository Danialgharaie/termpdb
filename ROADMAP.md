# TermPDB Improvement Roadmap & Feature Backlog

A living backlog of technical capabilities, improvements, and architectural milestones for [`termpdb`](Cargo.toml). Features from this document can be drawn and scoped into dedicated specifications and implementation plans.

---

## Table of Contents
1. [Visual & Rendering](#1-visual--rendering)
2. [Performance & Optimization](#2-performance--optimization)
3. [Scientific Tooling & Analysis](#3-scientific-tooling--analysis)
4. [Formats & Data Sources](#4-formats--data-sources)
5. [TUI, UX & CLI](#5-tui-ux--cli)
6. [Status Legend](#status-legend)

---

## 1. Visual & Rendering

*Enhancing visual clarity, graphical depth, and rich terminal character-cell graphics.*

### 1.1 Specular Shading & Blinn-Phong Material Highlights
- **Status:** `[Done]` | **Priority:** High
- **Description:** Add Blinn-Phong specular reflection highlights $(k_s (\vec{N} \cdot \vec{H})^\alpha)$ to atoms, bonds, and secondary structure ribbons for glossy/metallic 3D depth and curved surface glinting.
- **Approach:**
  - Update [`src/render/lighting.rs`](src/render/lighting.rs) to accept eye/view vector $\vec{V}$ and compute half-vector $\vec{H}$.
  - Apply shininess and specular intensity across sphere rasterizer and ribbon shaders.

### 1.2 Silhouette & Depth Outlining (Cartoon Cel-Shading)
- **Status:** `[Done]` | **Priority:** High
- **Description:** Add crisp dark silhouette outlines around molecular contours and depth boundaries between overlapping chains and active sites.
- **Approach:**
  - Post-process the float depth buffer ($Z$-buffer) in [`src/render/postprocess.rs`](src/render/postprocess.rs) using an adaptive Sobel/depth-difference kernel.
  - Darken silhouette pixel colors where depth jumps significantly or normal discontinuities occur.

### 1.3 Screen-Space Ambient Occlusion (SSAO / Pocket Shadowing)
- **Status:** `[Done]` | **Priority:** High
- **Description:** Softly darken deep binding pockets, clefts, and buried residue crevices to create strong spatial depth cues in low-resolution cell grids.
- **Approach:**
  - Screen-space depth sampling around each pixel in the framebuffer to compute an occlusion factor $k_{\text{ao}} \in [0.4, 1.0]$.
  - Modulate pixel RGB values by $k_{\text{ao}}$.

### 1.4 Braille Subpixel Canvas ($2 \times 4$ Dot Canvas)
- **Status:** `[Done]` | **Priority:** High
- **Description:** Provide an ultra-high resolution $2 \times 4$ subpixel mode using Unicode Braille characters (`⠁`–`⣿`) for fine wireframes, density cages, and dot-surfaces.
- **Approach:**
  - Add a Braille canvas buffer to [`src/render/braille.rs`](src/render/braille.rs) mapping 8 subpixels per character cell with 24-bit truecolor foreground styling.
  - Implement a `RenderMode::Wireframe` or Braille viewport toggle.

### 1.5 Non-Covalent Interactions (Hydrogen Bonds & Disulfide Bridges)
- **Status:** `[Done]` | **Priority:** High
- **Description:** Automatically detect and render polar donor-acceptor hydrogen bonds ($d \le 3.5\text{ \AA}$) and covalent disulfide `-S-S-` bridges with glowing dashed lines.
- **Approach:**
  - Add H-bond and disulfide detectors to [`src/model/interactions.rs`](src/model/interactions.rs).
  - Render stippled/dashed 3D lines in Ball & Stick, Trace, and Ribbon modes with distinct interaction colors.

### 1.6 Nucleic Acid Double-Helix Ribbons (DNA / RNA Ladders)
- **Status:** `[Done]` | **Priority:** High
- **Description:** Render publication-quality DNA and RNA double helices with smooth sugar-phosphate backbone ribbons and rectangular planar base-pair rungs.
- **Approach:**
  - Detect nucleic acid chains (phosphate $P$, ribose $C4'/C1'$, nitrogenous bases $N1/N9$).
  - Draw planar base-pair slabs colored by base identity (Adenine=Green, Thymine/Uracil=Red, Guanine=Yellow, Cytosine=Cyan).

### 1.7 AlphaFold pLDDT, Electrostatic, and Curated Color Palettes
- **Status:** `[Done]` | **Priority:** High
- **Description:** Expand [`ColorScheme`](src/render/color.rs) with official AlphaFold pLDDT confidence coloring, Electrostatic Potential (Red-White-Blue), and popular terminal themes (Catppuccin, Nord, Tokyo Night, Gruvbox).
- **Approach:**
  - Add `ColorScheme::Plddt`, `ColorScheme::Electrostatic`, and curated palette variants (`Catppuccin`, `Nord`, `TokyoNight`, `Gruvbox`).

### 1.8 Depth-of-Field (DoF) & Focal Plane Cueing
- **Status:** `[Done]` | **Priority:** Medium
- **Description:** Cinematic focal distance depth cueing to keep active sites crisp while softly fading or desaturating background loops.
- **Approach:**
  - Add focal plane parameter $Z_{\text{focus}}$ and focal depth range to camera and lighting pipeline.

---

## 2. Performance & Optimization

*Maximizing frame rates, multi-core scaling, and rendering massive macromolecular complexes (>100k atoms).*

### 2.1 Vectorized & Analytical Frustum-Culled Rasterizer
- **Status:** `[Done]` | **Priority:** High
- **Description:** Accelerate per-pixel sphere ray-intersection and cylinder testing with early screen-space frustum bounding culling, analytical depth calculations, and branchless $Z$-buffer tests.
- **Approach:**
  - Screen bounding box rejection prior to pixel loops in [`src/render/rasterizer.rs`](src/render/rasterizer.rs).
  - Fast single-precision inverse-radius and early nearest-$Z$ rejection.

### 2.2 Multi-Threaded Tile/Chunk Software Rasterizer (Rayon)
- **Status:** `[Done]` | **Priority:** High
- **Description:** Parallelize software rendering and post-processing across CPU cores for high-resolution displays and complex supramolecular assemblies.
- **Approach:**
  - Split the framebuffer into horizontal row bands ([`FramebufferBand`](src/render/buffer.rs)) and dispatch rasterization and SSAO post-processing across cores with `rayon`.

### 2.3 Spatial Grid Indexing & Occlusion Culling
- **Status:** `[Done]` | **Priority:** Medium
- **Description:** 3D uniform spatial hash grid ([`SpatialGrid`](src/model/spatial.rs)) for $O(1)$ voxel neighbor searches and buried core atom detection.
- **Approach:**
  - Uniform cubic spatial hashing for neighbor queries and buried-atom culling.

---

## 3. Scientific Tooling & Analysis

*Quantitative structural biology analysis, measurements, and structural comparisons.*

### 3.1 Structural Alignment & Superposition (Kabsch RMSD)
- **Status:** `[Done]` | **Priority:** High
- **Description:** Align two or more structures/chains in 3D space to compare wild-type vs. mutants or distinct conformational states.
- **Approach:**
  - Implemented the **Kabsch algorithm** (optimal rotation matrix via $3 \times 3$ SVD Jacobi eigen-decomposition in [`src/math/kabsch.rs`](src/math/kabsch.rs)).
  - Needleman-Wunsch dynamic programming sequence alignment for automatic $C\alpha$ pairing across mutations/insertions in [`src/model/align.rs`](src/model/align.rs).
  - CLI `--align` support with global RMSD and per-residue coordinate deviation calculations.

### 3.2 Advanced Angle, Dihedral & Contact Geometry Measurements
- **Status:** `[Done]` | **Priority:** Medium
- **Description:** 1-to-4 atom selection queue measuring 2-atom distance, 3-atom bond angle ($\theta^\circ$), 4-atom dihedral/torsion angle ($\phi/\psi/\chi^\circ$), and Ramachandran quadrant classification in [`src/model/geometry.rs`](src/model/geometry.rs).
- **Approach:**
  - Generalized 4-atom selection FIFO in [`src/select.rs`](src/select.rs).
  - Real-time HUD status and CLI reporting flags (`--angle`, `--dihedral`).

### 3.3 Secondary Structure Assignment (DSSP Algorithm)
- **Status:** `[Done]` | **Priority:** Medium
- **Description:** Calculates secondary structure on-the-fly when files lack `HELIX`/`SHEET` records (e.g. raw MD snapshots, plain PDBs, or AlphaFold predictions).
- **Approach:**
  - Implemented pure Rust DSSP electrostatic backbone hydrogen-bonding energy calculator ($E = q_1 q_2 [1/r_{ON} + 1/r_{CH} - 1/r_{OH} - 1/r_{CN}] \cdot 332 < -0.5\text{ kcal/mol}$) in [`src/model/dssp.rs`](src/model/dssp.rs).
  - Automatic fallback on parse + explicit `--dssp` CLI flag and recalculation.

---

## 4. Formats & Data Sources

*Expanding supported chemical and structural file formats and database integrations.*

### 4.1 Trajectory & Multi-Model MD Playback (DCD / XTC / NetCDF)
- **Status:** `[Backlog]` | **Priority:** High
- **Description:** Play molecular dynamics trajectories directly inside the terminal with timeline scrubbing and frame controls.
- **Approach:**
  - Support multi-model PDB/mmCIF trajectory stepping and binary MD formats (`.dcd`, `.xtc`).
  - Add playback controls (play, pause, step forward/back, FPS speed, looping).

### 4.2 Cryo-EM & X-Ray Density Maps (CCP4 / MRC Format)
- **Status:** `[Backlog]` | **Priority:** Medium
- **Description:** Parse electron density maps (`.mrc`, `.ccp4`, `.map`) and render 3D density isomesh wireframes or contours aligned with the coordinate model.
- **Approach:**
  - Stream header and 3D voxel grid from MRC/CCP4 files.
  - Generate wireframe contours at adjustable $\sigma$ threshold levels.

### 4.3 Direct AlphaFold DB & EMDB Fetching
- **Status:** `[Backlog]` | **Priority:** Low
- **Description:** Fetch predicted structures from AlphaFold DB via UniProt accession IDs (`AF-P00533-F1`) or EMDB density volumes via EMD IDs (`EMD-1234`).
- **Approach:**
  - Expand [`src/parser/rcsb.rs`](src/parser/rcsb.rs) into a general online fetcher supporting RCSB, AlphaFold DB API, and EMDB REST endpoints.

---

## 5. TUI, UX & CLI

*Improving ergonomics, terminal scripting, and interactive exploration.*

### 5.1 Interactive In-TUI Command Console (PyMOL-style `:` Bar)
- **Status:** `[Backlog]` | **Priority:** High
- **Description:** Provide an interactive command line inside the TUI for fine-grained selections and operations.
- **Approach:**
  - Implement command bar opened with `:` key.
  - Support commands such as:
    - `select <name> <expression>` (e.g., `select active :A and resi 50-75 and not name H*`)
    - `color <scheme|hex> <selection>`
    - `show <mode> <selection>` / `hide <selection>`
    - `center <selection>` / `zoom <selection>`

### 5.2 Session State Saving & Restoration (`.tpdb` / JSON)
- **Status:** `[Backlog]` | **Priority:** Low
- **Description:** Save and reload viewer states (camera orientation, representations, color schemes, custom selections, and measurements).
- **Approach:**
  - Serialize app state into JSON/TOML configuration files for reproducible scientific presentation and sharing.

---

## Status Legend
- `[Backlog]`: Scoped and ready to be turned into a design spec and implementation plan.
- `[In Progress]`: Active implementation or partial functionality present in tree.
- `[Done]`: Fully implemented, tested, and integrated.

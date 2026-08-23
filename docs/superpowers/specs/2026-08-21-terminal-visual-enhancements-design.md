# Terminal 3D Visual Enhancements Design Specification

**Date:** 2026-08-21  
**Status:** Approved  
**Language/Framework:** Rust (2024 Edition), `ratatui`, `crossterm`

---

## 1. Overview & Objectives

This specification defines the architecture, data models, algorithms, and test suites for 8 major visual enhancements designed specifically for **pure ANSI/Unicode terminal character-cell environments** (tmux, herdr, standard SSH):

1. **Specular Shading & Blinn-Phong Highlights**: Dynamic specular reflections with variable shininess and light half-vector computation on spheres, cylinders, and ribbon quads.
2. **Silhouette & Depth Outlining (Cartoon Cel-Shading)**: Screen-space depth-discontinuity and normal-edge filter on the Z-buffer to draw crisp 1-pixel dark borders around overlapping chains and molecular contours.
3. **Screen-Space Ambient Occlusion (SSAO)**: Depth-difference crevice darkening to make binding pockets, active site clefts, and buried residue contacts pop with realistic contact shadows.
4. **Braille Subpixel Canvas ($2 \times 4$ Dot Canvas)**: 8-subpixel-per-cell Unicode Braille canvas for ultra-crisp wireframe rendering, C-alpha backbones, and dot-surface representations.
5. **Non-Covalent Interactions**: Automatic detection and glowing dashed/stippled 3D line rendering for hydrogen bonds ($D-\text{H}\cdots A \le 3.5\text{ \AA}$) and covalent disulfide bridges (`-S-S-`).
6. **Nucleic Acid Double-Helix Ribbons**: Dedicated cartoon ladder representation for DNA and RNA: smooth sugar-phosphate backbone ribbon and rectangular planar base-pair slabs with standard nucleotide coloring (A, C, G, T, U).
7. **AlphaFold pLDDT, Electrostatic & Curated Palettes**: Official pLDDT confidence color map, Red-White-Blue electrostatic ramp, and curated terminal themes (Catppuccin, Nord, Tokyo Night, Gruvbox).
8. **Depth-of-Field (DoF) & Focal Plane Cueing**: Distance-weighted focal depth falloff that keeps selected residues/active sites crisp while softly cueing distant regions.

---

## 2. Component Design & Interfaces

### 2.1 Lighting & Shading (`src/render/lighting.rs`)
- Add specular reflection:
  $$\vec{H} = \text{normalize}(\vec{L} + \vec{V})$$
  $$I_{\text{spec}} = k_s \cdot \max(0, \vec{N} \cdot \vec{H})^\alpha$$
  $$C_{\text{final}} = \text{clamp}(C_{\text{diffuse}} + C_{\text{spec}} \cdot I_{\text{spec}}, 0, 255) \times \text{fog}(z) \times \text{ao}(x, y)$$
- `Lighting` struct maintains configurable `specular_intensity`, `shininess`, `ambient_occlusion_enabled`, `outline_enabled`, and `dof_focus_depth`.

### 2.2 Depth Outlines & SSAO Post-Processing (`src/render/buffer.rs` / `src/render/postprocess.rs`)
- `Framebuffer::apply_postprocessing(&mut self, lighting: &Lighting)`:
  - **Edge / Silhouette detection**: 4-neighbor / 8-neighbor depth difference $\Delta z = \sum |z_{i,j} - z_{neighbor}| / z_{i,j}$. If $\Delta z > \text{threshold}$, scale pixel luminance by $0.25$ (dark outline).
  - **SSAO factor**: Sample circular screen-space disc (radius 2-3 pixels) comparing center depth with sample depth to calculate local occlusion factor $k_{\text{ao}} \in [0.4, 1.0]$.

### 2.3 Braille Subpixel Canvas (`src/render/braille.rs`)
- `BrailleBuffer` struct mapping $(W, H)$ character cells to $(2W, 4H)$ binary/colored subpixels.
- Unicode offset formula: `0x2800 | bitmask` where bits `0..7` map to positions:
  $$\begin{bmatrix} (0,0): 0x01 & (1,0): 0x08 \\ (0,1): 0x02 & (1,1): 0x10 \\ (0,2): 0x04 & (1,2): 0x20 \\ (0,3): 0x40 & (1,3): 0x80 \end{bmatrix}$$
- Integrates with `RenderMode::Wireframe` and standalone Braille string generation.

### 2.4 Non-Covalent Interactions (`src/model/interactions.rs` & `src/render/representations/interactions.rs`)
- `InteractionDetector`:
  - Disulfide bonds: Cys-SG to Cys-SG ($d \le 2.2\text{ \AA}$).
  - Hydrogen bonds: Donor (N, O) to Acceptor (N, O) distance $\le 3.5\text{ \AA}$.
- Rendered as dashed 3D line segments with periodic blank gaps along the line interpolation.

### 2.5 Nucleic Acid Ladder Generator (`src/render/representations/nucleic.rs`)
- Recognizes purines (A, G) and pyrimidines (C, T, U, DA, DC, DG, DT).
- Traces sugar-phosphate backbone using Catmull-Rom spline on $P$ and $C4'$ atoms.
- Generates planar base slabs between backbone and base ring centroid.

### 2.6 Color Schemes (`src/render/color.rs`)
- Add `ColorScheme::Plddt` (AlphaFold confidence: Blue $\ge 90$, Cyan $70-90$, Yellow $50-70$, Orange $< 50$).
- Add `ColorScheme::Electrostatic` (Red-White-Blue gradient).
- Add curated themes: `ColorScheme::Catppuccin`, `ColorScheme::Nord`, `ColorScheme::TokyoNight`, `ColorScheme::Gruvbox`.

---

## 3. Verification & Testing

- Unit tests for Blinn-Phong math, half-vector computation, and specular saturation.
- Unit tests for Sobel/depth-difference outline detector and SSAO kernel.
- Unit tests for Braille bitmask generation and $2 \times 4$ coordinate mapping.
- Unit tests for H-bond and disulfide detection on real standard PDB structures (`1CRN`, `1BNA`).
- Unit tests for DNA/RNA ladder geometry generation on nucleic structures.
- Unit tests for new color palettes and CLI flag arguments.

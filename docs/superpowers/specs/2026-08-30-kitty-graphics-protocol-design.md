# Specification: Native Kitty Graphics Protocol Integration for TermPDB

**Date:** 2026-08-30  
**Status:** Approved  
**Author:** Antigravity & User  

---

## 1. Overview & Goals

TermPDB currently renders 3D macromolecular structures in terminal character cells using ANSI 24-bit half-blocks (`▀` / `▄`), where each terminal cell represents $1 \times 2$ pixels.

Modern terminal emulators (Kitty, Ghostty, WezTerm, Konsole, Alacritty with patches) implement the **Kitty Graphics Protocol**, enabling true pixel graphics to be drawn directly onto the terminal canvas.

This specification defines the native integration of the Kitty Graphics Protocol into TermPDB, providing:
1. **High-Resolution 3D Molecular Rasterization**: True 1:1 screen pixel rendering for spheres, ribbons, cylinders, non-covalent bonds, Blinn-Phong specular highlights, cel outlines, and SSAO pocket shadowing.
2. **Interactive TUI Graphics with Overlay Support**: The Kitty graphic renders directly within the central 3D viewport behind Ratatui text widgets (Header, Footer, HUD, Modal Dialogs).
3. **Runtime Graphics Backend Toggling**: Seamless switching between ANSI half-block mode and Kitty pixel mode with a single keypress (`K` or `g`).
4. **Calibrated Mouse Interaction**: High-resolution raycast atom picking mapped from terminal character cell mouse events.
5. **Headless Scriptable Export**: `--export-kitty` flag for embedding high-resolution molecular images into terminal scripts, CLI tools, and `fzf` previews.

---

## 2. Architecture & Components

```
+-------------------------------------------------------------------------+
|                                TermPDB                                  |
+-------------------------------------------------------------------------+
                                    |
            +-----------------------+-----------------------+
            |                                               |
  [ GraphicsBackend::HalfBlock ]                 [ GraphicsBackend::Kitty ]
            |                                               |
   Framebuffer: (Cols x Rows*2)                  Framebuffer: (Cols*cell_w x Rows*cell_h)
            |                                               |
  Software 3D Rasterizer                          Software 3D Rasterizer
            |                                               |
  Ratatui Buffer (Half-Block Cells)               Kitty Protocol In-band RGBA Base64
            |                                               |
            +-----------------------+-----------------------+
                                    |
                            Terminal Display
```

### 2.1 Graphics Backend Mode
A new enum `GraphicsBackend` will be introduced in `src/render/mod.rs` (or `src/render/graphics.rs`):
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GraphicsBackend {
    #[default]
    HalfBlock,
    Kitty,
}
```

- **CLI Flag**: `--kitty` / `--graphics <auto|kitty|halfblock>` in `src/cli.rs`.
- **Runtime Toggle**: Handled in `src/tui/events.rs` and `src/tui/app.rs` via `KeyCode::Char('K')` or `KeyCode::Char('g')`.
- **State Cleanup**: Whenever toggling away from Kitty mode or on application termination, a Kitty delete escape sequence (`\x1b_Ga=d,d=a\x1b\\` or specific image ID deletion `\x1b_Ga=d,d=i,i=1\x1b\\`) is sent to the terminal to remove lingering graphics.

---

## 3. Resolution, Pixel Geometry & Framebuffer Sizing

### 3.1 Terminal Window & Cell Pixel Querying
1. TermPDB queries the terminal window pixel dimensions using `crossterm::terminal::window_size()` (or fallback `ioctl(TIOCGWINSZ)`).
2. Given `window_pixel_width`, `window_pixel_height`, `window_cols`, and `window_rows`:
   $$\text{cell\_width\_px} = \max\left(1, \left\lfloor \frac{\text{window\_pixel\_width}}{\text{window\_cols}} \right\rfloor\right)$$
   $$\text{cell\_height\_px} = \max\left(1, \left\lfloor \frac{\text{window\_pixel\_height}}{\text{window\_rows}} \right\rfloor\right)$$
3. If the terminal does not report pixel dimensions (e.g. returns 0 or an error), a sensible default cell geometry of $10 \times 20$ pixels is used.

### 3.2 Viewport Framebuffer Allocation
When rendering the central 3D viewport area of size `Rect { x, y, width: viewport_cols, height: viewport_rows }`:
- **In `HalfBlock` Mode**:
  $$\text{fb\_width} = \text{viewport\_cols}$$
  $$\text{fb\_height} = \text{viewport\_rows} \times 2$$
- **In `Kitty` Mode**:
  $$\text{fb\_width} = \text{viewport\_cols} \times \text{cell\_width\_px}$$
  $$\text{fb\_height} = \text{viewport\_rows} \times \text{cell\_height\_px}$$

The software rasterizer (`rasterize_scene`, `render_structure_ctx`, etc.) operates transparently on `fb_width` and `fb_height`, resulting in crisp high-DPI output without altering the 3D projection math.

---

## 4. Kitty Protocol Transmission & In-TUI Placement

### 4.1 In-Band RGBA Stream Protocol
TermPDB transmits raw 32-bit RGBA pixel buffers directly via in-band terminal escape sequences:
- **Format**: `f=32` (32-bit RGBA).
- **Action & Replacement**: `a=T` (transmit and display immediately) with image ID `i=1` and quiet mode `q=2` (suppress terminal response messages).
- **Placement**:
  - `c=<viewport_cols>` (column count).
  - `r=<viewport_rows>` (row count).
  - `X=<area.x>` and `Y=<area.y>` (cursor cell position placement).
  - `z=-1` (negative z-index so text cells render on top).
- **Chunking**: Data is base64-encoded and transmitted in chunks of up to 4096 bytes:
  - First chunk: `\x1b_Ga=T,f=32,s=<w>,v=<h>,c=<cols>,r=<rows>,z=-1,i=1,q=2,m=1;<chunk>\x1b\\`
  - Intermediate chunks: `\x1b_Gm=1;<chunk>\x1b\\`
  - Final chunk: `\x1b_Gm=0;<chunk>\x1b\\`

### 4.2 Framebuffer to RGBA Conversion
`Framebuffer` provides a high-performance method to dump its pixel colors into a flat 32-bit RGBA byte slice (`Vec<u8>` or reusable buffer) for Base64 encoding.

### 4.3 TUI Overlay Integration
- In `tui/widgets/viewport.rs`, when `backend == GraphicsBackend::Kitty`, character cells inside the viewport area are filled with transparent/empty cells or spaces to reserve the screen region.
- Ratatui's Header, Footer, HUD, Modal Help, and Modal Info widgets continue rendering standard text cells.
- The Kitty graphic escape sequence is flushed to stdout immediately following the frame buffer draw, positioning the 3D graphic precisely inside the viewport rect.

---

## 5. Mouse Interaction & Atom Picking

In interactive mode, mouse click events report character cell coordinates `(click_col, click_row)`.

When `backend == GraphicsBackend::Kitty`:
$$\text{px} = (\text{click\_col} - \text{viewport\_rect.x}) \times \text{cell\_width\_px} + \frac{\text{cell\_width\_px}}{2}$$
$$\text{py} = (\text{click\_row} - \text{viewport\_rect.y}) \times \text{cell\_height\_px} + \frac{\text{cell\_height\_px}}{2}$$

`pick_atom_at_screen(px, py, framebuffer, structure, camera)` performs the exact depth-buffer raycasting at `(px, py)` in the high-res framebuffer, providing pixel-accurate atom selection.

---

## 6. Headless Export: `--export-kitty`

TermPDB adds `--export-kitty [FILE|-]` CLI flag in `src/cli.rs` and `src/render/export.rs`:
- Renders the structure using the active or specified `--mode`, `--color`, lighting, and SSAO.
- Formats the resulting framebuffer into the Kitty graphics escape sequence (`f=32,a=T`).
- Writes to stdout or the designated file for headless terminal automation, piping, and documentation generators.

---

## 7. Testing & Validation Plan

1. **Unit Tests**:
   - `test_kitty_escape_generation`: Verify correct formatting of `\x1b_G` escape codes, chunking with `m=1`/`m=0`, and Base64 encoding.
   - `test_cell_pixel_calculation`: Verify correct cell geometry derivation and fallback handling.
   - `test_mouse_coordinate_mapping`: Verify screen pixel conversion in Kitty vs. HalfBlock modes.
2. **Integration Tests**:
   - Verify `--export-kitty -` produces valid Kitty protocol data on stdout for PDB/mmCIF inputs.
   - Verify `--kitty` launches without panic and cleans up Kitty image IDs on exit.
3. **Interactive Verification**:
   - Test rotation, spin, zoom, and pan in Kitty mode.
   - Test switching modes dynamically (`K` / `g`).
   - Verify modals (Help `?`, Info `i`, `/` atom pick prompt) overlay properly above the 3D graphics canvas.

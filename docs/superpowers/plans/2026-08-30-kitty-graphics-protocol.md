# Kitty Graphics Protocol Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement native Kitty Graphics Protocol support in TermPDB for high-resolution 3D macromolecular rendering, interactive TUI graphics, runtime backend toggling (`K`), calibrated mouse picking, and headless `--export-kitty` export.

**Architecture:** Add `GraphicsBackend` (`HalfBlock` and `Kitty`) to the rendering pipeline. In Kitty mode, compute the native terminal cell pixel dimensions and allocate a full-resolution pixel `Framebuffer`. Encode the framebuffer directly into chunked in-band 32-bit RGBA Base64 escape sequences (`\x1b_Ga=T,f=32...`), positioned behind Ratatui text widgets.

**Tech Stack:** Rust (2024 edition), Ratatui, Crossterm, Base64, Rayon.

## Global Constraints
- Pure Rust / standard library + existing project dependencies (no external C libraries).
- In-band RGBA Base64 stream (`f=32`), 4096-byte chunking (`m=1`/`m=0`), quiet mode `q=2`, negative z-index `z=-1`.
- Graceful fallback to default cell dimensions ($10 \times 20$ px) when terminal window pixel size is unsupported.
- Clean terminal state restoration (emit Kitty delete sequence on exit or mode toggle).

---

### Task 1: Framebuffer RGBA Pixel Extraction

**Files:**
- Modify: `src/render/buffer.rs`
- Test: `tests/test_buffer_rgba.rs`

**Interfaces:**
- Consumes: `Framebuffer` and `PixelColor` from `src/render/buffer.rs`
- Produces: `Framebuffer::to_rgba_bytes(&self) -> Vec<u8>` and `Framebuffer::write_rgba_bytes(&self, out: &mut Vec<u8>)`

- [ ] **Step 1: Write the failing test**

Create `tests/test_buffer_rgba.rs`:
```rust
use termpdb::render::{Framebuffer, PixelColor};

#[test]
fn test_framebuffer_to_rgba_bytes() {
    let mut fb = Framebuffer::new(2, 2);
    fb.set(0, 0, 1.0, PixelColor(255, 0, 0));
    fb.set(1, 0, 1.0, PixelColor(0, 255, 0));
    fb.set(0, 1, 1.0, PixelColor(0, 0, 255));
    fb.set(1, 1, 1.0, PixelColor(255, 255, 255));

    let rgba = fb.to_rgba_bytes();
    assert_eq!(rgba.len(), 2 * 2 * 4);
    // (0,0): R=255, G=0, B=0, A=255
    assert_eq!(&rgba[0..4], &[255, 0, 0, 255]);
    // (1,0): R=0, G=255, B=0, A=255
    assert_eq!(&rgba[4..8], &[0, 255, 0, 255]);
    // (0,1): R=0, G=0, B=255, A=255
    assert_eq!(&rgba[8..12], &[0, 0, 255, 255]);
    // (1,1): R=255, G=255, B=255, A=255
    assert_eq!(&rgba[12..16], &[255, 255, 255, 255]);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test test_buffer_rgba`  
Expected: FAIL with "no method named `to_rgba_bytes` found for struct `Framebuffer`"

- [ ] **Step 3: Implement `to_rgba_bytes` and `write_rgba_bytes` on `Framebuffer`**

In `src/render/buffer.rs`, add:
```rust
impl Framebuffer {
    /// Copies the framebuffer pixel colors into a flat 32-bit RGBA byte vector.
    pub fn to_rgba_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(self.width * self.height * 4);
        self.write_rgba_bytes(&mut bytes);
        bytes
    }

    /// Appends the framebuffer pixel colors as 32-bit RGBA (4 bytes per pixel, A=255) into the provided buffer.
    pub fn write_rgba_bytes(&self, out: &mut Vec<u8>) {
        out.clear();
        out.reserve(self.width * self.height * 4);
        for pixel in &self.color {
            out.push(pixel.0);
            out.push(pixel.1);
            out.push(pixel.2);
            out.push(255);
        }
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test test_buffer_rgba`  
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/render/buffer.rs tests/test_buffer_rgba.rs
git commit -m "feat(render): add RGBA byte serialization to Framebuffer"
```

---

### Task 2: Kitty Graphics Protocol Escape Sequence Encoder & Cell Geometry

**Files:**
- Create: `src/render/kitty.rs`
- Modify: `src/render/mod.rs`
- Test: `tests/test_kitty_protocol.rs`

**Interfaces:**
- Consumes: Raw RGBA byte slices and terminal geometry
- Produces:
  - `GraphicsBackend` enum (`HalfBlock`, `Kitty`)
  - `encode_kitty_graphics_rgba(w: u32, h: u32, cols: u16, rows: u16, x: u16, y: u16, z_index: i32, image_id: u32, rgba: &[u8]) -> String`
  - `encode_kitty_delete(image_id: Option<u32>) -> &'static str`
  - `get_terminal_cell_size() -> (u32, u32)`

- [ ] **Step 1: Write the failing test**

Create `tests/test_kitty_protocol.rs`:
```rust
use termpdb::render::{GraphicsBackend, encode_kitty_delete, encode_kitty_graphics_rgba};

#[test]
fn test_encode_kitty_graphics_single_chunk() {
    let rgba = vec![255, 0, 0, 255, 0, 255, 0, 255]; // 2 pixels
    let seq = encode_kitty_graphics_rgba(2, 1, 10, 5, 0, 1, -1, 1, &rgba);

    assert!(seq.starts_with("\x1b_G"));
    assert!(seq.contains("a=T"));
    assert!(seq.contains("f=32"));
    assert!(seq.contains("s=2,v=1"));
    assert!(seq.contains("c=10,r=5"));
    assert!(seq.contains("X=0,Y=1"));
    assert!(seq.contains("z=-1"));
    assert!(seq.contains("i=1"));
    assert!(seq.contains("q=2"));
    assert!(seq.ends_with("\x1b\\"));
}

#[test]
fn test_encode_kitty_delete() {
    assert_eq!(encode_kitty_delete(Some(1)), "\x1b_Ga=d,d=i,i=1,q=2\x1b\\");
    assert_eq!(encode_kitty_delete(None), "\x1b_Ga=d,d=a,q=2\x1b\\");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test test_kitty_protocol`  
Expected: FAIL with unresolved symbols `encode_kitty_graphics_rgba` and `encode_kitty_delete`

- [ ] **Step 3: Implement `src/render/kitty.rs` and update `src/render/mod.rs`**

Create `src/render/kitty.rs`:
```rust
//! Kitty Graphics Protocol escape sequence generator and terminal geometry utilities.

use std::fmt::Write;

/// Active graphics rendering backend.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GraphicsBackend {
    /// ANSI 24-bit half-block characters (▀/▄) mapped to terminal cells.
    #[default]
    HalfBlock,
    /// High-resolution true-pixel rasterization via Kitty Graphics Protocol.
    Kitty,
}

impl GraphicsBackend {
    pub fn is_kitty(&self) -> bool {
        matches!(self, Self::Kitty)
    }

    pub fn toggle(&mut self) {
        *self = match self {
            Self::HalfBlock => Self::Kitty,
            Self::Kitty => Self::HalfBlock,
        };
    }
}

/// Fallback character cell pixel dimensions (width, height) when terminal query is unavailable.
pub const DEFAULT_CELL_PIXEL_WIDTH: u32 = 10;
pub const DEFAULT_CELL_PIXEL_HEIGHT: u32 = 20;
pub const KITTY_CHUNK_SIZE: usize = 4096;

/// Queries the terminal for the pixel size of a single character cell.
pub fn get_terminal_cell_size() -> (u32, u32) {
    if let Ok(size) = crossterm::terminal::window_size() {
        if size.width > 0 && size.height > 0 && size.columns > 0 && size.rows > 0 {
            let cell_w = (size.width as u32 / size.columns as u32).max(1);
            let cell_h = (size.height as u32 / size.rows as u32).max(1);
            return (cell_w, cell_h);
        }
    }
    (DEFAULT_CELL_PIXEL_WIDTH, DEFAULT_CELL_PIXEL_HEIGHT)
}

/// Generates Kitty Graphics Protocol escape sequences for displaying an RGBA buffer.
pub fn encode_kitty_graphics_rgba(
    width: u32,
    height: u32,
    cols: u16,
    rows: u16,
    x: u16,
    y: u16,
    z_index: i32,
    image_id: u32,
    rgba: &[u8],
) -> String {
    use base64::Engine;
    let b64 = base64::engine::general_purpose::STANDARD.encode(rgba);
    let bytes = b64.as_bytes();
    let total_len = bytes.len();
    let mut out = String::with_capacity(total_len + 256);

    let mut offset = 0;
    let mut first = true;

    while offset < total_len {
        let chunk_end = (offset + KITTY_CHUNK_SIZE).min(total_len);
        let chunk_str = std::str::from_utf8(&bytes[offset..chunk_end]).unwrap_or("");
        let has_more = chunk_end < total_len;
        let m = if has_more { 1 } else { 0 };

        if first {
            let _ = write!(
                out,
                "\x1b_Ga=T,f=32,s={width},v={height},c={cols},r={rows},X={x},Y={y},z={z_index},i={image_id},q=2,m={m};{chunk_str}\x1b\\"
            );
            first = false;
        } else {
            let _ = write!(out, "\x1b_Gm={m};{chunk_str}\x1b\\");
        }

        offset = chunk_end;
    }

    out
}

/// Generates a Kitty delete escape sequence.
pub fn encode_kitty_delete(image_id: Option<u32>) -> &'static str {
    match image_id {
        Some(1) => "\x1b_Ga=d,d=i,i=1,q=2\x1b\\",
        Some(_) => "\x1b_Ga=d,d=a,q=2\x1b\\",
        None => "\x1b_Ga=d,d=a,q=2\x1b\\",
    }
}
```

In `src/render/mod.rs`, add:
```rust
pub mod kitty;
pub use kitty::{
    DEFAULT_CELL_PIXEL_HEIGHT, DEFAULT_CELL_PIXEL_WIDTH, GraphicsBackend, encode_kitty_delete,
    encode_kitty_graphics_rgba, get_terminal_cell_size,
};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test test_kitty_protocol`  
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/render/kitty.rs src/render/mod.rs tests/test_kitty_protocol.rs
git commit -m "feat(render): add Kitty Graphics Protocol encoder and cell geometry query"
```

---

### Task 3: CLI Flags & Headless Export (`--kitty`, `--export-kitty`)

**Files:**
- Modify: `src/cli.rs`
- Modify: `src/render/export.rs`
- Modify: `src/main.rs`
- Test: `tests/test_kitty_export.rs`

**Interfaces:**
- Consumes: `render_structure_to_framebuffer` from `src/render/mod.rs` and `encode_kitty_graphics_rgba`
- Produces: `--kitty` and `--export-kitty <PATH>` CLI options; `export_kitty_frame(...)` in `src/render/export.rs`

- [ ] **Step 1: Write the failing test**

Create `tests/test_kitty_export.rs`:
```rust
use termpdb::model::Structure;
use termpdb::render::{ExportConfig, GraphicsBackend, export_kitty_frame};

#[test]
fn test_export_kitty_frame_generates_valid_sequence() {
    let structure = Structure::default();
    let config = ExportConfig::default();
    let output = export_kitty_frame(&structure, &config, 80, 40).expect("kitty export failed");

    assert!(output.starts_with("\x1b_G"));
    assert!(output.contains("a=T"));
    assert!(output.contains("f=32"));
    assert!(output.ends_with("\x1b\\"));
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test test_kitty_export`  
Expected: FAIL with unresolved function `export_kitty_frame`

- [ ] **Step 3: Implement `export_kitty_frame` and CLI flags**

In `src/render/export.rs`, add:
```rust
/// Exports a rendered structure directly as a Kitty Graphics Protocol escape sequence string.
pub fn export_kitty_frame(
    structure: &Structure,
    config: &ExportConfig,
    cols: u16,
    rows: u16,
) -> Result<String, crate::error::TermPdbError> {
    let (cell_w, cell_h) = crate::render::get_terminal_cell_size();
    let pixel_w = (cols as u32 * cell_w).max(1);
    let pixel_h = (rows as u32 * cell_h).max(1);

    let fb = render_structure_to_framebuffer(structure, config, pixel_w as usize, pixel_h as usize);
    let rgba = fb.to_rgba_bytes();
    let seq = crate::render::encode_kitty_graphics_rgba(
        pixel_w, pixel_h, cols, rows, 0, 0, 0, 1, &rgba,
    );
    Ok(seq)
}
```

In `src/cli.rs`:
Add `--kitty` and `--export-kitty`:
```rust
    /// Enable high-resolution Kitty Graphics Protocol rendering
    #[arg(long)]
    pub kitty: bool,

    /// Export rendered frame as Kitty Graphics Protocol escape sequence to file or stdout (-)
    #[arg(long, value_name = "FILE")]
    pub export_kitty: Option<PathBuf>,
```

In `src/main.rs`:
Handle `export_kitty` when present:
```rust
    if let Some(ref path) = args.export_kitty {
        let (cols, rows) = crossterm::terminal::size().unwrap_or((80, 24));
        let kitty_str = termpdb::render::export_kitty_frame(&structure, &export_config, cols, rows)?;
        if path.as_os_str() == "-" {
            print!("{}", kitty_str);
        } else {
            std::fs::write(path, kitty_str)?;
        }
        return Ok(());
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test test_kitty_export`  
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/cli.rs src/render/export.rs src/main.rs tests/test_kitty_export.rs
git commit -m "feat(cli): add --kitty and --export-kitty CLI flags and headless export"
```

---

### Task 4: Interactive TUI Integration, Viewport Sizing & Keybinding Toggle

**Files:**
- Modify: `src/tui/app.rs`
- Modify: `src/tui/events.rs`
- Modify: `src/tui/widgets/viewport.rs`
- Modify: `src/tui/widgets/help.rs`
- Test: `tests/test_tui_kitty.rs`

**Interfaces:**
- Consumes: `GraphicsBackend`, `encode_kitty_graphics_rgba`, `encode_kitty_delete`
- Produces:
  - `app.graphics_backend` state
  - Key `K` toggle action
  - Kitty sequence rendering during `app.draw()`
  - High-res calibrated mouse picking

- [ ] **Step 1: Write the failing test**

Create `tests/test_tui_kitty.rs`:
```rust
use termpdb::model::Structure;
use termpdb::render::GraphicsBackend;
use termpdb::tui::app::App;

#[test]
fn test_app_graphics_backend_initialization_and_toggle() {
    let structure = Structure::default();
    let mut app = App::new(structure, false, 1.0);
    assert_eq!(app.graphics_backend, GraphicsBackend::HalfBlock);

    app.toggle_graphics_backend();
    assert_eq!(app.graphics_backend, GraphicsBackend::Kitty);
    assert!(app.needs_rerender);

    app.toggle_graphics_backend();
    assert_eq!(app.graphics_backend, GraphicsBackend::HalfBlock);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --test test_tui_kitty`  
Expected: FAIL with "no field `graphics_backend` on `App`"

- [ ] **Step 3: Integrate `graphics_backend` into `App`, `events.rs`, and `viewport.rs`**

1. In `src/tui/events.rs`:
   - Add `AppAction::ToggleGraphicsBackend`.
   - Map `KeyCode::Char('K')` to `AppAction::ToggleGraphicsBackend`.
   - Update mouse picking calculation when `app.graphics_backend.is_kitty()`:
     ```rust
     let (cell_w, cell_h) = crate::render::get_terminal_cell_size();
     let px = (col.saturating_sub(app.viewport_area.x) as u32 * cell_w + cell_w / 2) as usize;
     let py = (row.saturating_sub(app.viewport_area.y) as u32 * cell_h + cell_h / 2) as usize;
     ```

2. In `src/tui/app.rs`:
   - Add `pub graphics_backend: GraphicsBackend` to `App`.
   - In `App::new`, initialize `graphics_backend: GraphicsBackend::HalfBlock` (or based on CLI `--kitty`).
   - Add method `pub fn toggle_graphics_backend(&mut self)`.
   - In `App::resize_framebuffer(&mut self, width: u16, height: u16)`:
     ```rust
     if self.graphics_backend.is_kitty() {
         let (cell_w, cell_h) = crate::render::get_terminal_cell_size();
         let pixel_w = (width as u32 * cell_w).max(1) as usize;
         let pixel_h = (height as u32 * cell_h).max(1) as usize;
         self.framebuffer.resize(pixel_w, pixel_h);
     } else {
         self.framebuffer.resize(width as usize, (height * 2) as usize);
     }
     ```
   - In `App::draw`:
     After Ratatui renders the terminal widgets, if `graphics_backend.is_kitty()`, construct and flush the Kitty escape sequence:
     ```rust
     if self.graphics_backend.is_kitty() && self.viewport_area.width > 0 && self.viewport_area.height > 0 {
         let rgba = self.framebuffer.to_rgba_bytes();
         let seq = crate::render::encode_kitty_graphics_rgba(
             self.framebuffer.width as u32,
             self.framebuffer.height as u32,
             self.viewport_area.width,
             self.viewport_area.height,
             self.viewport_area.x,
             self.viewport_area.y,
             -1,
             1,
             &rgba,
         );
         use std::io::Write;
         let _ = std::io::stdout().write_all(seq.as_bytes());
         let _ = std::io::stdout().flush();
     }
     ```
   - On cleanup/exit (or in `Drop` / shutdown hook), emit `encode_kitty_delete(None)` to stdout.

3. In `src/tui/widgets/viewport.rs`:
   - When `app.graphics_backend.is_kitty()`, fill viewport area with transparent/empty cells so Ratatui doesn't overwrite the terminal graphics layer with half-blocks.

4. In `src/tui/widgets/help.rs`:
   - Add `K: toggle Kitty graphics / half-block mode` to the help dialog.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test --test test_tui_kitty`  
Expected: PASS

- [ ] **Step 5: Run all test suites to ensure zero regressions**

Run: `cargo test`  
Expected: ALL PASS

- [ ] **Step 6: Commit**

```bash
git add src/tui/app.rs src/tui/events.rs src/tui/widgets/viewport.rs src/tui/widgets/help.rs tests/test_tui_kitty.rs
git commit -m "feat(tui): integrate Kitty graphics into TUI viewport with runtime 'K' toggle and mouse picking"
```

---

### Task 5: Documentation & Verification

**Files:**
- Modify: `README.md`
- Modify: `ROADMAP.md`

- [ ] **Step 1: Update README and ROADMAP**
  - Add `--kitty` and `--export-kitty` documentation in `README.md`.
  - Document `K` keybinding in `README.md` Interactive Controls table.
  - Update `ROADMAP.md` under Section 1 (Visual & Rendering) marking Kitty Graphics Protocol as `[Done]`.

- [ ] **Step 2: Verify build and clippy**
  Run: `cargo clippy --all-targets -- -D warnings && cargo test`  
  Expected: Clean build and all tests pass.

- [ ] **Step 3: Commit**
```bash
git add README.md ROADMAP.md
git commit -m "docs: document Kitty graphics protocol support and update roadmap"
```

//! Truecolor half-block framebuffer and depth buffer.
//!
//! Provides a software framebuffer that maps 2 vertical pixels per terminal character
//! using Unicode upper half block character (`▀`), where the upper half is colored using
//! the foreground color and the lower half is colored using the background color.

/// RGB color represented as `(r, g, b)` components in 0..=255.
pub type PixelColor = (u8, u8, u8);

/// Hard ceiling on framebuffer pixels (`width * height`).
///
/// A frame allocates ~7 bytes per pixel (RGB tuple + f32 depth), so this caps
/// a single framebuffer at roughly 448 MiB -- far above any legitimate
/// terminal grid or export resolution (4K @ 2x SSAA is ~33 Mpx), while
/// preventing hostile `--width`/`--height`/`--ssaa` combinations from
/// requesting terabytes of RAM and aborting the process with OOM.
pub const MAX_FRAMEBUFFER_PIXELS: usize = 64 * 1024 * 1024;

/// A software framebuffer containing truecolor RGB pixels and a floating-point depth buffer (Z-buffer).
#[derive(Debug, Clone, PartialEq)]
pub struct Framebuffer {
    /// Width in pixels (equivalent to terminal columns)
    pub width: usize,
    /// Height in pixels (equivalent to 2 * terminal rows)
    pub height: usize,
    /// Color buffer
    pub pixels: Vec<PixelColor>,
    /// Depth buffer (Z-buffer, closer objects have smaller z)
    pub depth: Vec<f32>,
}

impl Framebuffer {
    /// Creates a new `Framebuffer` with specified pixel dimensions.
    pub fn new(width: usize, height: usize) -> Self {
        let size = width * height;
        Self {
            width,
            height,
            pixels: vec![(0, 0, 0); size],
            depth: vec![f32::INFINITY; size],
        }
    }

    /// Resizes the framebuffer, clearing previous contents to black / infinite depth.
    pub fn resize(&mut self, width: usize, height: usize) {
        self.width = width;
        self.height = height;
        let size = width * height;
        self.pixels = vec![(0, 0, 0); size];
        self.depth = vec![f32::INFINITY; size];
    }

    /// Clears the color buffer to `bg_color` and depth buffer to positive infinity.
    pub fn clear(&mut self, bg_color: PixelColor) {
        self.pixels.fill(bg_color);
        self.depth.fill(f32::INFINITY);
    }

    /// Sets a pixel at integer screen coordinates `(x, y)` with depth `z`.
    ///
    /// If `(x, y)` is within bounds and `z < current_depth`, writes `color` and `z`
    /// into the buffers and returns `true`. Otherwise returns `false`.
    pub fn set_pixel(&mut self, x: i32, y: i32, z: f32, color: PixelColor) -> bool {
        if x < 0 || y < 0 {
            return false;
        }
        let ux = x as usize;
        let uy = y as usize;
        if ux >= self.width || uy >= self.height {
            return false;
        }

        let idx = uy * self.width + ux;
        if z < self.depth[idx] {
            self.depth[idx] = z;
            self.pixels[idx] = color;
            true
        } else {
            false
        }
    }

    /// Gets the pixel color at `(x, y)` if within bounds.
    pub fn get_pixel(&self, x: usize, y: usize) -> Option<PixelColor> {
        if x < self.width && y < self.height {
            Some(self.pixels[y * self.width + x])
        } else {
            None
        }
    }

    /// Gets the depth value at `(x, y)` if within bounds.
    pub fn get_depth(&self, x: usize, y: usize) -> Option<f32> {
        if x < self.width && y < self.height {
            Some(self.depth[y * self.width + x])
        } else {
            None
        }
    }

    /// Returns the terminal cell character and colors for terminal coordinates `(col, row)`.
    ///
    /// Top pixel is at `(col, row * 2)` (foreground).
    /// Bottom pixel is at `(col, row * 2 + 1)` (background).
    /// Returns `('▀', fg_color, Some(bg_color))` or `('▀', fg_color, None)` if bottom pixel is out of range.
    pub fn cell_at(&self, col: usize, row: usize) -> (char, PixelColor, Option<PixelColor>) {
        let y_top = row * 2;
        let y_bottom = row * 2 + 1;

        let fg = self.get_pixel(col, y_top).unwrap_or((0, 0, 0));
        let bg = self.get_pixel(col, y_bottom);

        ('▀', fg, bg)
    }

    /// Returns the full 2D grid of half-block cells formatted for terminal output.
    ///
    /// Dimension is `[terminal_rows][terminal_cols]`.
    pub fn get_half_blocks(&self) -> Vec<Vec<(char, PixelColor, Option<PixelColor>)>> {
        let term_rows = self.height.div_ceil(2);
        let mut rows = Vec::with_capacity(term_rows);

        for r in 0..term_rows {
            let mut row_cells = Vec::with_capacity(self.width);
            for c in 0..self.width {
                row_cells.push(self.cell_at(c, r));
            }
            rows.push(row_cells);
        }

        rows
    }

    /// Returns the terminal dimensions `(cols, rows)` represented by this framebuffer.
    pub fn terminal_size(&self) -> (usize, usize) {
        (self.width, self.height.div_ceil(2))
    }

    /// Converts the framebuffer contents to a standalone ANSI truecolor string.
    pub fn to_ansi(&self) -> String {
        let term_rows = self.height.div_ceil(2);
        if self.width == 0 || term_rows == 0 {
            return String::new();
        }

        use std::fmt::Write;
        // Worst case is roughly two SGR sequences (~20 bytes) plus one char per
        // cell; run-length encoding makes the common case (large uniform regions
        // such as the background and solid sphere interiors) far smaller.
        // Reserve lazily beyond 1 MiB: `with_capacity` pre-allocates eagerly,
        // and a worst-case reservation for a huge framebuffer would itself be
        // a multi-hundred-MB spike before a single cell is written. String
        // growth is amortized O(1), so under-reserving only costs copies.
        let worst_case = term_rows.saturating_mul(self.width).saturating_mul(24);
        let mut out = String::with_capacity(worst_case.min(1024 * 1024));

        for r in 0..term_rows {
            // Each row is terminated by an SGR reset, so the run state starts fresh.
            let mut run_fg: Option<PixelColor> = None;
            let mut run_bg: Option<Option<PixelColor>> = None;

            for c in 0..self.width {
                let (ch, fg, bg) = self.cell_at(c, r);

                // Only re-emit a color code when it changes from the previous
                // cell in this row, instead of per cell.
                if run_fg != Some(fg) {
                    let _ = write!(out, "\x1b[38;2;{};{};{}m", fg.0, fg.1, fg.2);
                    run_fg = Some(fg);
                }
                if run_bg != Some(bg) {
                    match bg {
                        Some(bg_color) => {
                            let _ = write!(
                                out,
                                "\x1b[48;2;{};{};{}m",
                                bg_color.0, bg_color.1, bg_color.2
                            );
                        }
                        None => out.push_str("\x1b[49m"),
                    }
                    run_bg = Some(bg);
                }
                out.push(ch);
            }
            out.push_str("\x1b[0m\n");
        }

        out
    }

    /// Splits the framebuffer into horizontal row bands of `band_height` rows each for parallel processing.
    pub fn par_bands_mut(&mut self, band_height: usize) -> Vec<FramebufferBand<'_>> {
        if self.width == 0 || self.height == 0 || band_height == 0 {
            return Vec::new();
        }

        let width = self.width;
        let band_size = width * band_height;

        let pixel_chunks: Vec<&mut [PixelColor]> = self.pixels.chunks_mut(band_size).collect();
        let depth_chunks: Vec<&mut [f32]> = self.depth.chunks_mut(band_size).collect();

        pixel_chunks
            .into_iter()
            .zip(depth_chunks)
            .enumerate()
            .map(|(i, (pixels, depth))| {
                let actual_height = pixels.len() / width;
                FramebufferBand {
                    width,
                    height: actual_height,
                    y_offset: i * band_height,
                    pixels,
                    depth,
                }
            })
            .collect()
    }
}

/// A horizontal slice/band of a Framebuffer for parallel multi-threaded rendering.
pub struct FramebufferBand<'a> {
    pub width: usize,
    pub height: usize,
    pub y_offset: usize,
    pub pixels: &'a mut [PixelColor],
    pub depth: &'a mut [f32],
}

impl<'a> FramebufferBand<'a> {
    /// Sets a pixel at local coordinates `(local_x, local_y)` relative to this band.
    #[inline]
    pub fn set_pixel(&mut self, local_x: i32, local_y: i32, z: f32, color: PixelColor) -> bool {
        if local_x < 0 || local_y < 0 {
            return false;
        }
        let ux = local_x as usize;
        let uy = local_y as usize;
        if ux >= self.width || uy >= self.height {
            return false;
        }

        let idx = uy * self.width + ux;
        if z < self.depth[idx] {
            self.depth[idx] = z;
            self.pixels[idx] = color;
            true
        } else {
            false
        }
    }
}

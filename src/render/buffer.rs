//! Truecolor half-block framebuffer and depth buffer.
//!
//! Provides a software framebuffer that maps 2 vertical pixels per terminal character
//! using Unicode upper half block character (`▀`), where the upper half is colored using
//! the foreground color and the lower half is colored using the background color.

/// RGB color represented as `(r, g, b)` components in 0..=255.
pub type PixelColor = (u8, u8, u8);

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
        let mut out = String::with_capacity(term_rows * self.width * 35);

        for r in 0..term_rows {
            for c in 0..self.width {
                let (ch, fg, bg) = self.cell_at(c, r);
                let _ = write!(out, "\x1b[38;2;{};{};{}m", fg.0, fg.1, fg.2);
                if let Some(bg) = bg {
                    let _ = write!(out, "\x1b[48;2;{};{};{}m", bg.0, bg.1, bg.2);
                } else {
                    out.push_str("\x1b[49m");
                }
                out.push(ch);
            }
            out.push_str("\x1b[0m\n");
        }

        out
    }
}

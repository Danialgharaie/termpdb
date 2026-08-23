//! Braille subpixel canvas (2x4 dots per character cell).
//!
//! Maps an ultra-high resolution subpixel grid into Unicode Braille patterns (0x2800..0x28FF)
//! with depth buffering and 24-bit truecolor foreground styling.

use crate::render::buffer::PixelColor;

/// Unicode Braille bitmask table for 2x4 subpixels:
/// (dx, dy) where dx in 0..2, dy in 0..4
const BRAILLE_MASKS: [[u8; 4]; 2] = [
    [0x01, 0x02, 0x04, 0x40], // column 0: dot 1, 2, 3, 7
    [0x08, 0x10, 0x20, 0x80], // column 1: dot 4, 5, 6, 8
];

/// A software buffer providing 2x4 subpixel resolution using Unicode Braille characters.
#[derive(Debug, Clone, PartialEq)]
pub struct BrailleBuffer {
    /// Subpixel width (2 * terminal columns)
    pub width: usize,
    /// Subpixel height (4 * terminal rows)
    pub height: usize,
    /// Subpixel dot presence
    pub dots: Vec<bool>,
    /// Per-subpixel color
    pub colors: Vec<PixelColor>,
    /// Per-subpixel depth (Z-buffer)
    pub depth: Vec<f32>,
}

impl BrailleBuffer {
    /// Creates a new `BrailleBuffer` with given subpixel dimensions.
    pub fn new(subpixel_width: usize, subpixel_height: usize) -> Self {
        let size = subpixel_width * subpixel_height;
        Self {
            width: subpixel_width,
            height: subpixel_height,
            dots: vec![false; size],
            colors: vec![(0, 0, 0); size],
            depth: vec![f32::INFINITY; size],
        }
    }

    /// Clears the buffer.
    pub fn clear(&mut self) {
        self.dots.fill(false);
        self.colors.fill((0, 0, 0));
        self.depth.fill(f32::INFINITY);
    }

    /// Sets a subpixel at integer coordinate `(x, y)` with depth `z`.
    pub fn set_subpixel(&mut self, x: i32, y: i32, z: f32, color: PixelColor) -> bool {
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
            self.dots[idx] = true;
            self.colors[idx] = color;
            true
        } else {
            false
        }
    }

    /// Draws a 3D line in subpixel coordinates with linear depth interpolation.
    pub fn draw_line_3d(&mut self, p1: (f32, f32, f32), p2: (f32, f32, f32), color: PixelColor) {
        let (x1, y1, z1) = p1;
        let (x2, y2, z2) = p2;

        let dx = x2 - x1;
        let dy = y2 - y1;
        let dz = z2 - z1;

        let steps = (dx.abs().max(dy.abs()).ceil() as usize).max(1);
        let steps_f = steps as f32;

        let x_step = dx / steps_f;
        let y_step = dy / steps_f;
        let z_step = dz / steps_f;

        for i in 0..=steps {
            let fi = i as f32;
            let x = (x1 + fi * x_step).round() as i32;
            let y = (y1 + fi * y_step).round() as i32;
            let z = z1 + fi * z_step;
            self.set_subpixel(x, y, z, color);
        }
    }

    /// Returns the character and dominant color for the cell at `(col, row)`.
    #[allow(clippy::needless_range_loop)]
    pub fn cell_at(&self, col: usize, row: usize) -> (char, PixelColor) {
        let start_x = col * 2;
        let start_y = row * 4;

        let mut mask: u8 = 0;
        let mut min_z = f32::INFINITY;
        let mut best_color = (180, 180, 180);

        for dx in 0..2 {
            for dy in 0..4 {
                let sx = start_x + dx;
                let sy = start_y + dy;
                if sx < self.width && sy < self.height {
                    let idx = sy * self.width + sx;
                    if self.dots[idx] {
                        mask |= BRAILLE_MASKS[dx][dy];
                        if self.depth[idx] < min_z {
                            min_z = self.depth[idx];
                            best_color = self.colors[idx];
                        }
                    }
                }
            }
        }

        let ch = if mask == 0 {
            ' '
        } else {
            char::from_u32(0x2800 | mask as u32).unwrap_or(' ')
        };

        (ch, best_color)
    }

    /// Formats the BrailleBuffer into an ANSI string with 24-bit truecolor foreground colors.
    pub fn to_ansi(&self) -> String {
        let term_cols = self.width.div_ceil(2);
        let term_rows = self.height.div_ceil(4);
        if term_cols == 0 || term_rows == 0 {
            return String::new();
        }

        use std::fmt::Write;
        let mut out = String::with_capacity(term_rows * term_cols * 20);

        for r in 0..term_rows {
            let mut last_color: Option<PixelColor> = None;
            for c in 0..term_cols {
                let (ch, col) = self.cell_at(c, r);
                if ch == ' ' {
                    out.push(' ');
                    continue;
                }
                if last_color != Some(col) {
                    let _ = write!(out, "\x1b[38;2;{};{};{}m", col.0, col.1, col.2);
                    last_color = Some(col);
                }
                out.push(ch);
            }
            out.push_str("\x1b[0m\n");
        }

        out
    }
}

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

/// Fallback character cell pixel dimensions (width, height) for crisp high-DPI rendering when terminal query is unavailable.
pub const DEFAULT_CELL_PIXEL_WIDTH: u32 = 16;
pub const DEFAULT_CELL_PIXEL_HEIGHT: u32 = 32;
pub const KITTY_CHUNK_SIZE: usize = 4096;

/// Queries the terminal for the pixel size of a single character cell, scaled by `scale_factor`.
pub fn get_terminal_cell_size_scaled(scale_factor: f32) -> (u32, u32) {
    let scale = scale_factor.clamp(0.25, 8.0);
    if let Ok(size) = crossterm::terminal::window_size()
        && size.width > 0
        && size.height > 0
        && size.columns > 0
        && size.rows > 0
    {
        let cell_w = ((size.width as f32 / size.columns as f32) * scale).round() as u32;
        let cell_h = ((size.height as f32 / size.rows as f32) * scale).round() as u32;
        return (cell_w.max(2), cell_h.max(4));
    }
    (
        ((DEFAULT_CELL_PIXEL_WIDTH as f32) * scale).round().max(2.0) as u32,
        ((DEFAULT_CELL_PIXEL_HEIGHT as f32) * scale)
            .round()
            .max(4.0) as u32,
    )
}

/// Queries the terminal for the pixel size of a single character cell.
pub fn get_terminal_cell_size() -> (u32, u32) {
    get_terminal_cell_size_scaled(1.0)
}

/// Generates Kitty Graphics Protocol escape sequences for displaying an RGBA buffer compressed as fast in-memory PNG (`f=100`).
#[allow(clippy::too_many_arguments)]
pub fn encode_kitty_graphics_png(
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
    if rgba.is_empty() || width == 0 || height == 0 {
        return String::new();
    }

    let mut png_bytes = Vec::with_capacity(64 * 1024);
    {
        let mut encoder = png::Encoder::new(&mut png_bytes, width, height);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        encoder.set_compression(png::Compression::Fast);
        if let Ok(mut writer) = encoder.write_header() {
            if writer.write_image_data(rgba).is_err() {
                return encode_kitty_graphics_rgba(
                    width, height, cols, rows, x, y, z_index, image_id, rgba,
                );
            }
        } else {
            return encode_kitty_graphics_rgba(
                width, height, cols, rows, x, y, z_index, image_id, rgba,
            );
        }
    }

    use base64::Engine;
    let b64 = base64::engine::general_purpose::STANDARD.encode(&png_bytes);
    let bytes = b64.as_bytes();
    let total_len = bytes.len();
    let mut out = String::with_capacity(total_len + 256);

    // Position cursor at target cell (1-indexed row and col)
    let _ = write!(out, "\x1b[{};{}H", y + 1, x + 1);

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
                "\x1b_Ga=T,f=100,c={cols},r={rows},z={z_index},i={image_id},q=2,m={m};{chunk_str}\x1b\\"
            );
            first = false;
        } else {
            let _ = write!(out, "\x1b_Gm={m};{chunk_str}\x1b\\");
        }

        offset = chunk_end;
    }

    out
}

/// Generates Kitty Graphics Protocol escape sequences for displaying an RGBA buffer.
#[allow(clippy::too_many_arguments)]
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
    if rgba.is_empty() {
        return String::new();
    }
    let b64 = base64::engine::general_purpose::STANDARD.encode(rgba);
    let bytes = b64.as_bytes();
    let total_len = bytes.len();
    let mut out = String::with_capacity(total_len + 256);

    // Position cursor at target cell (1-indexed row and col)
    let _ = write!(out, "\x1b[{};{}H", y + 1, x + 1);

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
                "\x1b_Ga=T,f=32,s={width},v={height},c={cols},r={rows},z={z_index},i={image_id},q=2,m={m};{chunk_str}\x1b\\"
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

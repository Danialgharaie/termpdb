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
    if let Ok(size) = crossterm::terminal::window_size()
        && size.width > 0
        && size.height > 0
        && size.columns > 0
        && size.rows > 0
    {
        let cell_w = (size.width as u32 / size.columns as u32).max(1);
        let cell_h = (size.height as u32 / size.rows as u32).max(1);
        return (cell_w, cell_h);
    }
    (DEFAULT_CELL_PIXEL_WIDTH, DEFAULT_CELL_PIXEL_HEIGHT)
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

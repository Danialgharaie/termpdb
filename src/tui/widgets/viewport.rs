//! Viewport widget for rendering 3D framebuffer into Ratatui buffer.

use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::Color;
use ratatui::widgets::Widget;

use crate::render::Framebuffer;

/// Ratatui widget that displays a 3D software framebuffer using half-block Unicode characters.
pub struct ViewportWidget<'a> {
    framebuffer: &'a Framebuffer,
}

impl<'a> ViewportWidget<'a> {
    /// Creates a new `ViewportWidget` from a reference to a `Framebuffer`.
    pub fn new(framebuffer: &'a Framebuffer) -> Self {
        Self { framebuffer }
    }
}

impl Widget for ViewportWidget<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.width == 0 || area.height == 0 {
            return;
        }

        let term_cols = self.framebuffer.width.min(area.width as usize);
        let term_rows = self
            .framebuffer
            .height
            .div_ceil(2)
            .min(area.height as usize);

        for r in 0..term_rows {
            let y = area.y + r as u16;
            for c in 0..term_cols {
                let x = area.x + c as u16;
                let (ch, fg, bg) = self.framebuffer.cell_at(c, r);
                let cell = &mut buf[(x, y)];
                cell.set_char(ch);
                cell.set_fg(Color::Rgb(fg.0, fg.1, fg.2));
                if let Some(bg_color) = bg {
                    cell.set_bg(Color::Rgb(bg_color.0, bg_color.1, bg_color.2));
                } else {
                    cell.set_bg(Color::Reset);
                }
            }
        }
    }
}

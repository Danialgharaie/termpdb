//! HUD Header and Footer status bars for TermPDB.

use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Paragraph, Widget};

use crate::model::Structure;
use crate::render::{ColorScheme, RenderMode};

/// Header status widget displaying macromolecule identity and summary statistics.
pub struct HeaderWidget<'a> {
    structure: &'a Structure,
}

impl<'a> HeaderWidget<'a> {
    /// Creates a new `HeaderWidget` referencing the given structure.
    pub fn new(structure: &'a Structure) -> Self {
        Self { structure }
    }
}

impl Widget for HeaderWidget<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.width == 0 || area.height == 0 {
            return;
        }

        let id_str = self
            .structure
            .id_code
            .as_deref()
            .unwrap_or("termpdb");
        let title_str = if self.structure.title.is_empty() {
            "Macromolecule"
        } else {
            self.structure.title.as_str()
        };

        let mut spans = vec![
            Span::styled(
                format!(" [{id_str}] "),
                Style::default()
                    .fg(Color::Black)
                    .bg(Color::Cyan)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::raw(" "),
            Span::styled(
                title_str,
                Style::default()
                    .fg(Color::White)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(
                format!(
                    " | Chains: {} | Residues: {} | Atoms: {}",
                    self.structure.chain_count(),
                    self.structure.residue_count(),
                    self.structure.atom_count()
                ),
                Style::default().fg(Color::DarkGray),
            ),
        ];

        let resolution = self
            .structure
            .metadata
            .get("resolution")
            .or_else(|| self.structure.metadata.get("RESOLUTION"));

        if let Some(res) = resolution {
            spans.push(Span::styled(
                format!(" | Res: {res} Å"),
                Style::default().fg(Color::Green),
            ));
        }

        let line = Line::from(spans);
        let paragraph = Paragraph::new(line).style(Style::default().bg(Color::Black));
        paragraph.render(area, buf);
    }
}

/// Footer status widget displaying active mode, color scheme, spin status, and keybindings.
pub struct FooterWidget {
    mode: RenderMode,
    color_scheme: ColorScheme,
    auto_spin: bool,
    fps: f32,
}

impl FooterWidget {
    /// Creates a new `FooterWidget` with rendering state and FPS.
    pub fn new(
        mode: RenderMode,
        color_scheme: ColorScheme,
        auto_spin: bool,
        fps: f32,
    ) -> Self {
        Self {
            mode,
            color_scheme,
            auto_spin,
            fps,
        }
    }
}

impl Widget for FooterWidget {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.width == 0 || area.height == 0 {
            return;
        }

        let spin_status = if self.auto_spin { "ON" } else { "OFF" };
        let spin_color = if self.auto_spin {
            Color::Green
        } else {
            Color::DarkGray
        };

        let spans = vec![
            Span::styled(" [1-4] ", Style::default().fg(Color::Yellow)),
            Span::styled(
                format!("Mode: {} ", self.mode.name()),
                Style::default().fg(Color::White),
            ),
            Span::styled(" [c] ", Style::default().fg(Color::Yellow)),
            Span::styled(
                format!("Color: {} ", self.color_scheme.name()),
                Style::default().fg(Color::White),
            ),
            Span::styled(" [Space] ", Style::default().fg(Color::Yellow)),
            Span::styled(
                format!("Spin: {spin_status} "),
                Style::default()
                    .fg(spin_color)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(" [?] ", Style::default().fg(Color::Yellow)),
            Span::styled("Help ", Style::default().fg(Color::White)),
            Span::styled(" [i] ", Style::default().fg(Color::Yellow)),
            Span::styled("Info ", Style::default().fg(Color::White)),
            Span::styled(" [q] ", Style::default().fg(Color::Yellow)),
            Span::styled("Quit ", Style::default().fg(Color::White)),
            Span::styled(
                format!("| {:>4.1} FPS", self.fps),
                Style::default().fg(Color::Cyan),
            ),
        ];

        let line = Line::from(spans);
        let paragraph = Paragraph::new(line).style(Style::default().bg(Color::Black));
        paragraph.render(area, buf);
    }
}

//! Help modal popup widget.

use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, BorderType, Borders, Clear, Paragraph, Widget};

/// Centered help popup modal listing all keybindings and mouse controls.
#[derive(Debug, Clone, Copy, Default)]
pub struct HelpWidget;

impl HelpWidget {
    /// Creates a new `HelpWidget`.
    pub fn new() -> Self {
        Self
    }
}

/// Helper function to create a centered rectangle within the given area.
pub fn centered_rect(percent_x: u16, percent_y: u16, r: Rect) -> Rect {
    let popup_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - percent_y) / 2),
            Constraint::Percentage(percent_y),
            Constraint::Percentage((100 - percent_y) / 2),
        ])
        .split(r);

    Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage((100 - percent_x) / 2),
            Constraint::Percentage(percent_x),
            Constraint::Percentage((100 - percent_x) / 2),
        ])
        .split(popup_layout[1])[1]
}

impl Widget for HelpWidget {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let modal_area = centered_rect(75, 75, area);
        Clear.render(modal_area, buf);

        let block = Block::default()
            .title(Span::styled(
                " Controls & Keybindings (? / Esc to close) ",
                Style::default()
                    .fg(Color::Cyan)
                    .add_modifier(Modifier::BOLD),
            ))
            .borders(Borders::ALL)
            .border_type(BorderType::Rounded)
            .border_style(Style::default().fg(Color::Cyan))
            .style(Style::default().bg(Color::Black));

        let inner = block.inner(modal_area);
        block.render(modal_area, buf);

        let lines = vec![
            Line::from(vec![Span::styled(
                "  Representation Modes:",
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled(
                    "    1 / 2 / 3 / 4 / 5",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("  Trace / Ball & Stick / Ribbon / VDW / Wireframe"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    m / M        ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Next / Previous Representation"),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "  Color Schemes & Visuals:",
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled(
                    "    c / C        ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Cycle 12 Color Schemes (CPK, Rainbow, Chain, SS, pLDDT, Themes)"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    O            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Toggle Silhouette Depth Outlines"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    k            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Toggle Screen-Space Ambient Occlusion (SSAO)"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    e            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Toggle Non-Covalent Interactions (H-Bonds & Disulfides)"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    f            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Toggle Depth-of-Field (DoF) Focal Cueing"),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "  Animation & View:",
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled(
                    "    Space        ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Toggle Auto-Spin"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    + / - (= / _)",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Increase / Decrease Spin Speed"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    r            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Reset Camera View & Zoom"),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "  Camera Controls:",
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled(
                    "    Left Click + Drag / Arrows ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw(" Orbit / Rotate 3D camera"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    Right Click + Drag / WASD  ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw(" Pan camera target"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    Scroll Wheel / [ / ]       ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw(" Zoom In / Out"),
            ]),
            Line::from(""),
            Line::from(vec![Span::styled(
                "  General:",
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![Span::styled(
                "  Selection & Measure:",
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            )]),
            Line::from(vec![
                Span::styled(
                    "    /            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Pick atom by id (A:12 or A:12:CA), Enter to confirm"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    Click        ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Pick nearest atom (drag still orbits)"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    x            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Clear selection (2: dist, 3: angle, 4: dihedral/Ramachandran)"),
            ]),
            Line::from(""),
            Line::from(vec![
                Span::styled(
                    "    n / p        ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Next / Previous model (wrap)"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    l / L        ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    LOD Auto / Full / Backbone / CA (Auto: 25k / 80k atoms)"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    b / B        ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Next / Previous biological assembly (ASU ↔ copies)"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    o            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Toggle waters / solvent (hidden by default)"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    u            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Toggle hydrogens"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    i            ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Toggle Structure Details / Info"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    ? / h        ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Toggle this Help modal"),
            ]),
            Line::from(vec![
                Span::styled(
                    "    q / Esc      ",
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::raw("    Quit application"),
            ]),
        ];

        let paragraph = Paragraph::new(lines);
        paragraph.render(inner, buf);
    }
}

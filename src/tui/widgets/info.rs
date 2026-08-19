//! Structure Information modal popup widget.

use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, BorderType, Borders, Clear, Paragraph, Widget};

use crate::model::Structure;
use crate::tui::widgets::help::centered_rect;

/// Centered popup modal displaying detailed structure statistics.
pub struct InfoWidget<'a> {
    structure: &'a Structure,
}

impl<'a> InfoWidget<'a> {
    /// Creates a new `InfoWidget` referencing the given structure.
    pub fn new(structure: &'a Structure) -> Self {
        Self { structure }
    }
}

impl Widget for InfoWidget<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let modal_area = centered_rect(75, 75, area);
        Clear.render(modal_area, buf);

        let block = Block::default()
            .title(Span::styled(
                " Structure Details (i / Esc to close) ",
                Style::default()
                    .fg(Color::Green)
                    .add_modifier(Modifier::BOLD),
            ))
            .borders(Borders::ALL)
            .border_type(BorderType::Rounded)
            .border_style(Style::default().fg(Color::Green))
            .style(Style::default().bg(Color::Black));

        let inner = block.inner(modal_area);
        block.render(modal_area, buf);

        let id_code = self.structure.id_code.as_deref().unwrap_or("N/A");
        let title = if self.structure.title.is_empty() {
            "N/A"
        } else {
            self.structure.title.as_str()
        };

        let chain_ids: Vec<String> = self.structure.chains.iter().map(|c| c.id.clone()).collect();
        let chain_str = if chain_ids.is_empty() {
            "None".to_string()
        } else {
            chain_ids.join(", ")
        };

        // Secondary structure counts
        let mut helix_count = 0;
        let mut sheet_count = 0;
        let mut coil_count = 0;

        for chain in &self.structure.chains {
            for res in &chain.residues {
                match res.secondary_structure {
                    crate::model::SecondaryStructure::Helix => helix_count += 1,
                    crate::model::SecondaryStructure::Sheet => sheet_count += 1,
                    crate::model::SecondaryStructure::Coil => coil_count += 1,
                }
            }
        }

        let (min_b, max_b) = self.structure.b_factor_range();
        let hetatm_count = self.structure.atoms.iter().filter(|a| a.is_hetatm).count();
        let heavy_atom_count = self
            .structure
            .atoms
            .iter()
            .filter(|a| !a.is_hydrogen())
            .count();

        let mut lines = vec![
            Line::from(vec![
                Span::styled(
                    "  PDB ID: ",
                    Style::default()
                        .fg(Color::Yellow)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    id_code,
                    Style::default()
                        .fg(Color::White)
                        .add_modifier(Modifier::BOLD),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  Title:  ",
                    Style::default()
                        .fg(Color::Yellow)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(title, Style::default().fg(Color::White)),
            ]),
            Line::from(""),
            Line::from(vec![
                Span::styled(
                    "  Chains (",
                    Style::default()
                        .fg(Color::Cyan)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    format!("{}", self.structure.chain_count()),
                    Style::default().fg(Color::White),
                ),
                Span::styled(
                    "): ",
                    Style::default()
                        .fg(Color::Cyan)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(chain_str, Style::default().fg(Color::White)),
            ]),
            Line::from(vec![
                Span::styled(
                    "  Total Residues: ",
                    Style::default()
                        .fg(Color::Cyan)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    format!("{}", self.structure.residue_count()),
                    Style::default().fg(Color::White),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  Secondary Structure: ",
                    Style::default()
                        .fg(Color::Cyan)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    format!("{helix_count} Helix, {sheet_count} Sheet, {coil_count} Coil"),
                    Style::default().fg(Color::White),
                ),
            ]),
            Line::from(""),
            Line::from(vec![
                Span::styled(
                    "  Total Atoms:       ",
                    Style::default()
                        .fg(Color::Magenta)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    format!("{}", self.structure.atom_count()),
                    Style::default().fg(Color::White),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  Heavy Atoms:       ",
                    Style::default()
                        .fg(Color::Magenta)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    format!("{heavy_atom_count}"),
                    Style::default().fg(Color::White),
                ),
            ]),
            Line::from(vec![
                Span::styled(
                    "  Heteroatoms/Water: ",
                    Style::default()
                        .fg(Color::Magenta)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(format!("{hetatm_count}"), Style::default().fg(Color::White)),
            ]),
            Line::from(vec![
                Span::styled(
                    "  B-Factor Range:    ",
                    Style::default()
                        .fg(Color::Magenta)
                        .add_modifier(Modifier::BOLD),
                ),
                Span::styled(
                    format!("{min_b:.2} .. {max_b:.2} Å²"),
                    Style::default().fg(Color::White),
                ),
            ]),
        ];

        // Metadata additions
        if !self.structure.metadata.is_empty() {
            lines.push(Line::from(""));
            lines.push(Line::from(vec![Span::styled(
                "  Metadata:",
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            )]));
            for (k, v) in &self.structure.metadata {
                lines.push(Line::from(vec![
                    Span::styled(format!("    {k}: "), Style::default().fg(Color::DarkGray)),
                    Span::styled(v.clone(), Style::default().fg(Color::White)),
                ]));
            }
        }

        let paragraph = Paragraph::new(lines);
        paragraph.render(inner, buf);
    }
}

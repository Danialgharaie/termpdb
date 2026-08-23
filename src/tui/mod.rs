//! Interactive Terminal User Interface (TUI) for TermPDB.
//!
//! Powered by Ratatui and Crossterm, providing 60 FPS truecolor 3D rendering,
//! smooth orbit/turntable camera controls, HUD status bars, and modal overlays.

pub mod app;
pub mod events;
pub mod widgets;

use std::io::stdout;
use std::panic;
use std::time::{Duration, Instant};

use crossterm::cursor::Show;
use crossterm::event::{self, DisableMouseCapture, EnableMouseCapture, Event};
use crossterm::execute;
use crossterm::terminal::{
    EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode,
};
use ratatui::Terminal;
use ratatui::backend::CrosstermBackend;

pub use app::App;
pub use events::{AppAction, MouseState};
pub use widgets::{FooterWidget, HeaderWidget, HelpWidget, InfoWidget, ViewportWidget};

use crate::error::Result;
use crate::model::Structure;
use crate::render::{ColorScheme, LodMode, RenderMode, Visibility};

/// RAII Guard ensuring terminal raw mode and alternate screen are properly restored on exit/panic.
struct TerminalCleanupGuard;

impl Drop for TerminalCleanupGuard {
    fn drop(&mut self) {
        let _ = disable_raw_mode();
        let _ = execute!(stdout(), LeaveAlternateScreen, DisableMouseCapture, Show);
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    structure: Structure,
    initial_mode: RenderMode,
    initial_color: ColorScheme,
    auto_spin: bool,
    spin_speed: f32,
    visibility: Visibility,
    lod: LodMode,
    postprocess: crate::render::PostProcessConfig,
    show_interactions: bool,
    dof: Option<f32>,
) -> Result<()> {
    // Set panic hook to ensure terminal is cleaned up even if a panic occurs
    let original_hook = panic::take_hook();
    panic::set_hook(Box::new(move |panic_info| {
        let _ = disable_raw_mode();
        let _ = execute!(stdout(), LeaveAlternateScreen, DisableMouseCapture, Show);
        original_hook(panic_info);
    }));

    enable_raw_mode()?;
    let mut out = stdout();
    execute!(out, EnterAlternateScreen, EnableMouseCapture)?;
    let _guard = TerminalCleanupGuard;

    let backend = CrosstermBackend::new(out);
    let mut terminal = Terminal::new(backend)?;
    terminal.clear()?;

    let mut app = App::new(structure, initial_mode, initial_color, auto_spin)
        .with_spin_speed(spin_speed)
        .with_visibility(visibility)
        .with_lod(lod)
        .with_postprocess(postprocess)
        .with_interactions(show_interactions)
        .with_dof(dof);
    let mut last_frame_time = Instant::now();
    let frame_target = Duration::from_micros(16_667); // ~60 FPS

    while !app.should_quit {
        let now = Instant::now();
        let dt = now.duration_since(last_frame_time).as_secs_f32();
        last_frame_time = now;

        app.tick(dt);

        // Drain all pending input *before* drawing so a key/mouse event is
        // reflected in this frame. Non-blocking: returns immediately when empty.
        while event::poll(Duration::ZERO)? {
            match event::read()? {
                Event::Key(key) => app.handle_key(key),
                Event::Mouse(mouse) => app.handle_mouse(mouse),
                Event::Resize(_, _) => app.needs_redraw = true,
                _ => {}
            }
        }

        // Only repaint when something actually changed (input, spin, FPS readout
        // tick, resize). While idle the TUI does NOT redraw 60 fps, so CPU sits
        // near 0 instead of constantly diffing an unchanged buffer.
        if app.needs_redraw {
            terminal.draw(|f| app.render_ui(f))?;
            app.needs_redraw = false;
        }

        // Block for the remainder of the frame budget, or until input arrives.
        // poll returns immediately when an event is ready, so input wakes the
        // loop without waiting out the timeout.
        let poll_timeout = frame_target.saturating_sub(now.elapsed());
        if !poll_timeout.is_zero() {
            let _ = event::poll(poll_timeout)?;
        }
    }

    Ok(())
}

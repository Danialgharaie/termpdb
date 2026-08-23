use ratatui::Terminal;
use ratatui::backend::TestBackend;
use termpdb::parser::parse_pdb;
use termpdb::render::{ColorScheme, RenderMode};
use termpdb::tui::app::App;

fn main() {
    let pdb = std::fs::read_to_string("/tmp/4egk.pdb").expect("read");
    let structure = parse_pdb(&pdb).expect("parse");

    for &(mode, label) in &[
        (RenderMode::Ribbon, "Ribbon"),
        (RenderMode::Vdw, "VDW"),
        (RenderMode::BallAndStick, "Ball&Stick"),
        (RenderMode::Trace, "Trace"),
    ] {
        let mut app = App::new(structure.clone(), mode, ColorScheme::Cpk, true);
        let backend = TestBackend::new(160, 100);
        let mut terminal = Terminal::new(backend).unwrap();
        for _ in 0..10 {
            app.tick(0.016);
            let _ = terminal.draw(|f| app.render_ui(f));
        }
        let iters = 300;
        let t = std::time::Instant::now();
        for _ in 0..iters {
            app.tick(0.016);
            let _ = terminal.draw(|f| app.render_ui(f));
        }
        let d = t.elapsed() / iters;
        println!(
            "{:11} TUI frame (spin, 160x100): {:6.3} ms  ({:5.0} fps)",
            label,
            d.as_secs_f32() * 1000.0,
            1.0 / d.as_secs_f32()
        );
    }

    // Idle (no spin, no input) with the needs_redraw gate: the TUI does NOT
    // redraw every frame, so per-iteration cost is just tick (no draw, no raster).
    let mut app = App::new(
        structure.clone(),
        RenderMode::Ribbon,
        ColorScheme::Cpk,
        false,
    );
    let backend = TestBackend::new(160, 100);
    let mut terminal = Terminal::new(backend).unwrap();
    terminal.draw(|f| app.render_ui(f)).unwrap();
    app.needs_redraw = false;
    let iters = 2000;
    let t = std::time::Instant::now();
    let mut draws = 0u32;
    for _ in 0..iters {
        app.tick(0.016);
        if app.needs_redraw {
            terminal.draw(|f| app.render_ui(f)).unwrap();
            app.needs_redraw = false;
            draws += 1;
        }
    }
    let d = t.elapsed() / iters;
    println!(
        "{:11} TUI idle (gated, 160x100): {:6.4} ms/iter  (draws={}/{} -> {:.0}% fewer draws vs 60fps)",
        "Ribbon",
        d.as_secs_f32() * 1000.0,
        draws,
        iters,
        100.0 * (1.0 - draws as f32 / iters as f32)
    );
}

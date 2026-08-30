use crossterm::event::{KeyCode, KeyEvent, KeyEventKind, KeyEventState, KeyModifiers};
use ratatui::Terminal;
use ratatui::backend::TestBackend;
use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::widgets::Widget;
use termpdb::math::Vec3;
use termpdb::model::{Atom, Chain, Element, Residue, Structure};
use termpdb::render::{ColorScheme, Framebuffer, GraphicsBackend, RenderMode};
use termpdb::tui::app::App;
use termpdb::tui::events::{AppAction, handle_key_event};
use termpdb::tui::widgets::ViewportWidget;

fn create_test_structure() -> Structure {
    let mut structure = Structure::with_id("1CRN", "TEST PROTEIN");
    let mut chain = Chain::new("A");
    let mut r1 = Residue::new(1, "ALA", "A");

    let c_elem = Element {
        atomic_number: 6,
        symbol: "C",
        name: "Carbon",
        covalent_radius: 0.77,
        vdw_radius: 1.70,
        cpk_color: (144, 144, 144),
    };

    let a1 = Atom::new(
        0,
        1,
        "CA",
        c_elem,
        Vec3::new(0.0, 0.0, 0.0),
        20.0,
        "ALA",
        1,
        "A",
        false,
    );

    let idx1 = structure.add_atom(a1);
    r1.atom_indices.push(idx1);
    chain.residues.push(r1);
    structure.add_chain(chain);
    structure.build_bonds();
    structure
}

#[test]
fn test_app_graphics_backend_initialization_and_toggle() {
    let structure = Structure::default();
    let mut app = App::new(structure, RenderMode::Ribbon, ColorScheme::Cpk, false);
    assert_eq!(app.graphics_backend, GraphicsBackend::HalfBlock);

    app.toggle_graphics_backend();
    assert_eq!(app.graphics_backend, GraphicsBackend::Kitty);
    assert!(app.needs_rerender);
    assert!(app.needs_redraw);

    app.toggle_graphics_backend();
    assert_eq!(app.graphics_backend, GraphicsBackend::HalfBlock);
    assert!(app.needs_rerender);
    assert!(app.needs_redraw);
}

#[test]
fn test_key_event_k_toggles_graphics_backend() {
    let key_upper = KeyEvent {
        code: KeyCode::Char('K'),
        modifiers: KeyModifiers::SHIFT,
        kind: KeyEventKind::Press,
        state: KeyEventState::empty(),
    };
    let action_upper = handle_key_event(key_upper);
    assert_eq!(action_upper, AppAction::ToggleGraphicsBackend);

    let key_g = KeyEvent {
        code: KeyCode::Char('g'),
        modifiers: KeyModifiers::empty(),
        kind: KeyEventKind::Press,
        state: KeyEventState::empty(),
    };
    let action_g = handle_key_event(key_g);
    assert_eq!(action_g, AppAction::ToggleGraphicsBackend);

    let key_g_upper = KeyEvent {
        code: KeyCode::Char('G'),
        modifiers: KeyModifiers::SHIFT,
        kind: KeyEventKind::Press,
        state: KeyEventState::empty(),
    };
    let action_g_upper = handle_key_event(key_g_upper);
    assert_eq!(action_g_upper, AppAction::ToggleGraphicsBackend);

    let key_k_lower = KeyEvent {
        code: KeyCode::Char('k'),
        modifiers: KeyModifiers::empty(),
        kind: KeyEventKind::Press,
        state: KeyEventState::empty(),
    };
    let action_k_lower = handle_key_event(key_k_lower);
    assert_eq!(action_k_lower, AppAction::ToggleSsao);
}

#[test]
fn test_apply_action_toggle_graphics_backend() {
    let structure = Structure::default();
    let mut app = App::new(structure, RenderMode::Ribbon, ColorScheme::Cpk, false);
    assert_eq!(app.graphics_backend, GraphicsBackend::HalfBlock);

    app.apply_action(AppAction::ToggleGraphicsBackend);
    assert_eq!(app.graphics_backend, GraphicsBackend::Kitty);

    app.apply_action(AppAction::ToggleGraphicsBackend);
    assert_eq!(app.graphics_backend, GraphicsBackend::HalfBlock);
}

#[test]
fn test_resize_framebuffer_kitty_mode() {
    let structure = Structure::default();
    let mut app = App::new(structure, RenderMode::Ribbon, ColorScheme::Cpk, false);
    app.resize_framebuffer(80, 24);
    assert_eq!(app.framebuffer.width, 80);
    assert_eq!(app.framebuffer.height, 48);

    app.toggle_graphics_backend();
    app.resize_framebuffer(80, 24);
    let (cell_w, cell_h) = termpdb::render::get_terminal_cell_size();
    assert_eq!(app.framebuffer.width, (80 * cell_w) as usize);
    assert_eq!(app.framebuffer.height, (24 * cell_h) as usize);
}

#[test]
fn test_viewport_widget_kitty_mode() {
    let fb = Framebuffer::new(80, 48);
    let area = Rect::new(0, 0, 10, 5);
    let mut buf = Buffer::empty(area);

    let widget = ViewportWidget::new(&fb).with_backend(GraphicsBackend::Kitty);
    widget.render(area, &mut buf);

    for y in 0..5 {
        for x in 0..10 {
            let cell = &buf[(x, y)];
            assert_eq!(cell.symbol(), " ");
        }
    }
}

#[test]
fn test_app_render_ui_kitty_mode() {
    let backend = TestBackend::new(80, 24);
    let mut terminal = Terminal::new(backend).unwrap();

    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::BallAndStick, ColorScheme::Cpk, false)
        .with_graphics_backend(GraphicsBackend::Kitty);

    terminal.draw(|f| app.render_ui(f)).unwrap();

    let (cell_w, cell_h) = termpdb::render::get_terminal_cell_size();
    assert_eq!(
        app.framebuffer.width,
        (app.viewport_area.width as u32 * cell_w) as usize
    );
    assert_eq!(
        app.framebuffer.height,
        (app.viewport_area.height as u32 * cell_h) as usize
    );
}

#[test]
fn test_mouse_pick_at_cell_kitty_mode() {
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::BallAndStick, ColorScheme::Cpk, false)
        .with_graphics_backend(GraphicsBackend::Kitty);

    let backend = TestBackend::new(80, 24);
    let mut terminal = Terminal::new(backend).unwrap();
    terminal.draw(|f| app.render_ui(f)).unwrap();

    let center_col = app.viewport_area.x + app.viewport_area.width / 2;
    let center_row = app.viewport_area.y + app.viewport_area.height / 2;

    app.apply_action(AppAction::PickAt {
        col: center_col,
        row: center_row,
    });

    assert_eq!(app.selection.atoms(), &[0]);
}

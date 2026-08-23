use crossterm::event::{
    KeyCode, KeyEvent, KeyEventKind, KeyEventState, KeyModifiers, MouseButton, MouseEvent,
    MouseEventKind,
};
use ratatui::Terminal;
use ratatui::backend::TestBackend;
use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use termpdb::math::Vec3;
use termpdb::model::{Atom, Chain, Element, Residue, SecondaryStructure, Structure};
use termpdb::render::{ColorScheme, Framebuffer, RenderMode};
use termpdb::tui::app::App;
use termpdb::tui::events::{AppAction, MouseState, handle_key_event, handle_mouse_event};
use termpdb::tui::widgets::{FooterWidget, HeaderWidget, HelpWidget, InfoWidget, ViewportWidget};

fn create_test_structure() -> Structure {
    let mut structure = Structure::with_id("1CRN", "WATER PENTAPEPTIDE");
    structure
        .metadata
        .insert("resolution".to_string(), "1.50".to_string());
    structure
        .metadata
        .insert("method".to_string(), "X-RAY DIFFRACTION".to_string());

    let mut chain = Chain::new("A");
    let mut r1 = Residue::new(1, "THR", "A");
    r1.secondary_structure = SecondaryStructure::Helix;
    let mut r2 = Residue::new(2, "THR", "A");
    r2.secondary_structure = SecondaryStructure::Helix;
    let mut r3 = Residue::new(3, "CYS", "A");
    r3.secondary_structure = SecondaryStructure::Sheet;
    let mut r4 = Residue::new(4, "PRO", "A");
    r4.secondary_structure = SecondaryStructure::Coil;

    let c_elem = Element {
        atomic_number: 6,
        symbol: "C",
        name: "Carbon",
        covalent_radius: 0.77,
        vdw_radius: 1.70,
        cpk_color: (144, 144, 144),
    };
    let o_elem = Element {
        atomic_number: 8,
        symbol: "O",
        name: "Oxygen",
        covalent_radius: 0.73,
        vdw_radius: 1.52,
        cpk_color: (255, 13, 13),
    };

    let a1 = Atom::new(
        0,
        1,
        "CA",
        c_elem,
        Vec3::new(0.0, 0.0, 0.0),
        15.0,
        "THR",
        1,
        "A",
        false,
    );
    let a2 = Atom::new(
        1,
        2,
        "CA",
        c_elem,
        Vec3::new(1.0, 2.0, 3.0),
        20.0,
        "THR",
        2,
        "A",
        false,
    );
    let a3 = Atom::new(
        2,
        3,
        "CA",
        c_elem,
        Vec3::new(2.0, 4.0, 6.0),
        25.0,
        "CYS",
        3,
        "A",
        false,
    );
    let a4 = Atom::new(
        3,
        4,
        "CA",
        c_elem,
        Vec3::new(3.0, 6.0, 9.0),
        30.0,
        "PRO",
        4,
        "A",
        false,
    );
    let a5 = Atom::new(
        4,
        5,
        "O",
        o_elem,
        Vec3::new(4.0, 8.0, 12.0),
        10.0,
        "HOH",
        5,
        "A",
        true,
    );

    let idx1 = structure.add_atom(a1);
    let idx2 = structure.add_atom(a2);
    let idx3 = structure.add_atom(a3);
    let idx4 = structure.add_atom(a4);
    let idx5 = structure.add_atom(a5);

    r1.atom_indices.push(idx1);
    r2.atom_indices.push(idx2);
    r3.atom_indices.push(idx3);
    r4.atom_indices.push(idx4);

    let mut r5 = Residue::new(5, "HOH", "A");
    r5.atom_indices.push(idx5);

    chain.residues.push(r1);
    chain.residues.push(r2);
    chain.residues.push(r3);
    chain.residues.push(r4);
    chain.residues.push(r5);

    structure.add_chain(chain);
    structure.build_bonds();

    structure
}

#[test]
fn test_app_initialization() {
    let structure = create_test_structure();
    let app = App::new(
        structure.clone(),
        RenderMode::Ribbon,
        ColorScheme::Cpk,
        true,
    );

    assert_eq!(app.render_mode, RenderMode::Ribbon);
    assert_eq!(app.color_scheme, ColorScheme::Cpk);
    assert!(app.auto_spin);
    assert!(!app.show_help);
    assert!(!app.show_info);
    assert!(!app.should_quit);
    assert_eq!(app.structure.title, "WATER PENTAPEPTIDE");
    assert_eq!(app.structure.id_code.as_deref(), Some("1CRN"));
}

#[test]
fn test_app_mode_switching() {
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::Trace, ColorScheme::Cpk, false);

    assert_eq!(app.render_mode, RenderMode::Trace);
    app.next_mode();
    assert_eq!(app.render_mode, RenderMode::BallAndStick);
    app.next_mode();
    assert_eq!(app.render_mode, RenderMode::Ribbon);
    app.next_mode();
    assert_eq!(app.render_mode, RenderMode::Vdw);
    app.next_mode();
    assert_eq!(app.render_mode, RenderMode::Wireframe);
    app.next_mode();
    assert_eq!(app.render_mode, RenderMode::Trace);

    app.prev_mode();
    assert_eq!(app.render_mode, RenderMode::Wireframe);

    app.set_mode(RenderMode::Ribbon);
    assert_eq!(app.render_mode, RenderMode::Ribbon);
}

#[test]
fn test_app_color_scheme_switching() {
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::Ribbon, ColorScheme::Cpk, false);

    assert_eq!(app.color_scheme, ColorScheme::Cpk);
    app.next_color_scheme();
    assert_eq!(app.color_scheme, ColorScheme::Rainbow);
    app.next_color_scheme();
    assert_eq!(app.color_scheme, ColorScheme::Chain);
    app.prev_color_scheme();
    assert_eq!(app.color_scheme, ColorScheme::Rainbow);
}

#[test]
fn test_app_spin_toggles_and_speed() {
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::Ribbon, ColorScheme::Cpk, true);

    assert!(app.auto_spin);
    app.toggle_spin();
    assert!(!app.auto_spin);
    app.toggle_spin();
    assert!(app.auto_spin);

    let initial_speed = app.spin_speed;
    app.increase_spin_speed();
    assert!(app.spin_speed > initial_speed);
    app.decrease_spin_speed();
    assert_eq!(app.spin_speed, initial_speed);

    // Test with_spin_speed and set_spin_speed
    let app2 = App::new(
        create_test_structure(),
        RenderMode::Ribbon,
        ColorScheme::Cpk,
        true,
    )
    .with_spin_speed(3.5);
    assert_eq!(app2.spin_speed, 3.5);

    let mut app3 = App::new(
        create_test_structure(),
        RenderMode::Ribbon,
        ColorScheme::Cpk,
        true,
    );
    app3.set_spin_speed(2.0);
    assert_eq!(app3.spin_speed, 2.0);
}

#[test]
fn test_app_modal_toggles() {
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::Ribbon, ColorScheme::Cpk, false);

    assert!(!app.show_help);
    app.toggle_help();
    assert!(app.show_help);
    app.toggle_help();
    assert!(!app.show_help);

    assert!(!app.show_info);
    app.toggle_info();
    assert!(app.show_info);
    app.toggle_info();
    assert!(!app.show_info);
}

#[test]
fn test_app_camera_reset() {
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::Ribbon, ColorScheme::Cpk, false);

    let initial_distance = app.camera.distance;
    app.camera.zoom(2.0);
    app.camera.orbit(10.0, 5.0);
    assert_ne!(app.camera.distance, initial_distance);

    app.reset_camera();
    assert_eq!(app.camera.distance, initial_distance);
}

#[test]
fn test_key_event_mapping() {
    let make_key = |code: KeyCode, modifiers: KeyModifiers| KeyEvent {
        code,
        modifiers,
        kind: KeyEventKind::Press,
        state: KeyEventState::NONE,
    };

    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('q'), KeyModifiers::NONE)),
        AppAction::Quit
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Esc, KeyModifiers::NONE)),
        AppAction::Quit
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('c'), KeyModifiers::CONTROL)),
        AppAction::Quit
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char(' '), KeyModifiers::NONE)),
        AppAction::ToggleSpin
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('1'), KeyModifiers::NONE)),
        AppAction::SetRenderMode(RenderMode::Trace)
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('2'), KeyModifiers::NONE)),
        AppAction::SetRenderMode(RenderMode::BallAndStick)
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('3'), KeyModifiers::NONE)),
        AppAction::SetRenderMode(RenderMode::Ribbon)
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('4'), KeyModifiers::NONE)),
        AppAction::SetRenderMode(RenderMode::Vdw)
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('m'), KeyModifiers::NONE)),
        AppAction::NextRenderMode
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('M'), KeyModifiers::SHIFT)),
        AppAction::PrevRenderMode
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('c'), KeyModifiers::NONE)),
        AppAction::NextColorScheme
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('C'), KeyModifiers::SHIFT)),
        AppAction::PrevColorScheme
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('r'), KeyModifiers::NONE)),
        AppAction::ResetCamera
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('?'), KeyModifiers::NONE)),
        AppAction::ToggleHelp
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('h'), KeyModifiers::NONE)),
        AppAction::ToggleHelp
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('i'), KeyModifiers::NONE)),
        AppAction::ToggleInfo
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('n'), KeyModifiers::NONE)),
        AppAction::NextModel
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('p'), KeyModifiers::NONE)),
        AppAction::PrevModel
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('o'), KeyModifiers::NONE)),
        AppAction::ToggleWaters
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('u'), KeyModifiers::NONE)),
        AppAction::ToggleHydrogens
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('/'), KeyModifiers::NONE)),
        AppAction::StartPickPrompt
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('x'), KeyModifiers::NONE)),
        AppAction::ClearSelection
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('b'), KeyModifiers::NONE)),
        AppAction::NextAssembly
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('B'), KeyModifiers::SHIFT)),
        AppAction::PrevAssembly
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('l'), KeyModifiers::NONE)),
        AppAction::NextLod
    );
    assert_eq!(
        handle_key_event(make_key(KeyCode::Char('L'), KeyModifiers::SHIFT)),
        AppAction::PrevLod
    );
}

#[test]
fn test_mouse_event_mapping() {
    let mut mouse_state = MouseState::default();

    // Mouse down at (10, 10)
    let down_event = MouseEvent {
        kind: MouseEventKind::Down(MouseButton::Left),
        column: 10,
        row: 10,
        modifiers: KeyModifiers::NONE,
    };
    let action = handle_mouse_event(down_event, &mut mouse_state);
    assert_eq!(action, AppAction::None);
    assert!(mouse_state.is_left_down);
    assert_eq!(mouse_state.last_pos, Some((10, 10)));

    // Drag left click to (15, 12) -> dx = 5, dy = 2 -> Orbit
    let drag_event = MouseEvent {
        kind: MouseEventKind::Drag(MouseButton::Left),
        column: 15,
        row: 12,
        modifiers: KeyModifiers::NONE,
    };
    let action = handle_mouse_event(drag_event, &mut mouse_state);
    match action {
        AppAction::Orbit { dx, dy } => {
            assert_eq!(dx, 5.0);
            assert_eq!(dy, 2.0);
        }
        other => panic!("Expected Orbit, got {:?}", other),
    }
    assert_eq!(mouse_state.last_pos, Some((15, 12)));

    // Drag with Shift + Left click -> Pan
    let drag_shift_event = MouseEvent {
        kind: MouseEventKind::Drag(MouseButton::Left),
        column: 18,
        row: 12,
        modifiers: KeyModifiers::SHIFT,
    };
    let action = handle_mouse_event(drag_shift_event, &mut mouse_state);
    match action {
        AppAction::Pan { dx, dy } => {
            assert_eq!(dx, 3.0);
            assert_eq!(dy, 0.0);
        }
        other => panic!("Expected Pan, got {:?}", other),
    }

    // Scroll up / down -> Zoom
    let scroll_up = MouseEvent {
        kind: MouseEventKind::ScrollUp,
        column: 18,
        row: 12,
        modifiers: KeyModifiers::NONE,
    };
    assert_eq!(
        handle_mouse_event(scroll_up, &mut mouse_state),
        AppAction::Zoom { delta: 1.0 }
    );

    let scroll_down = MouseEvent {
        kind: MouseEventKind::ScrollDown,
        column: 18,
        row: 12,
        modifiers: KeyModifiers::NONE,
    };
    assert_eq!(
        handle_mouse_event(scroll_down, &mut mouse_state),
        AppAction::Zoom { delta: -1.0 }
    );

    let mut click_state = MouseState::default();
    let down = MouseEvent {
        kind: MouseEventKind::Down(MouseButton::Left),
        column: 20,
        row: 10,
        modifiers: KeyModifiers::NONE,
    };
    assert_eq!(handle_mouse_event(down, &mut click_state), AppAction::None);
    let up = MouseEvent {
        kind: MouseEventKind::Up(MouseButton::Left),
        column: 20,
        row: 10,
        modifiers: KeyModifiers::NONE,
    };
    assert_eq!(
        handle_mouse_event(up, &mut click_state),
        AppAction::PickAt { col: 20, row: 10 }
    );
}

#[test]
fn test_viewport_widget_render() {
    let mut fb = Framebuffer::new(4, 4); // 4 columns, 2 terminal rows
    fb.set_pixel(0, 0, 1.0, (255, 0, 0)); // Top-left foreground
    fb.set_pixel(0, 1, 1.0, (0, 255, 0)); // Top-left background

    let area = Rect::new(0, 0, 4, 2);
    let mut buffer = Buffer::empty(area);

    let widget = ViewportWidget::new(&fb);
    ratatui::widgets::Widget::render(widget, area, &mut buffer);

    let cell = &buffer[(0, 0)];
    assert_eq!(cell.symbol(), "▀");
    assert_eq!(cell.fg, ratatui::style::Color::Rgb(255, 0, 0));
    assert_eq!(cell.bg, ratatui::style::Color::Rgb(0, 255, 0));
}

#[test]
fn test_hud_header_widget_render() {
    let structure = create_test_structure();
    let area = Rect::new(0, 0, 80, 1);
    let mut buffer = Buffer::empty(area);

    let widget = HeaderWidget::new(&structure);
    ratatui::widgets::Widget::render(widget, area, &mut buffer);

    // Header should contain structure title or ID
    let text: String = (0..80).map(|x| buffer[(x, 0)].symbol()).collect();
    assert!(text.contains("1CRN") || text.contains("WATER PENTAPEPTIDE"));
}

#[test]
fn test_hud_footer_widget_render() {
    let area = Rect::new(0, 0, 80, 1);
    let mut buffer = Buffer::empty(area);

    let widget = FooterWidget::new(
        RenderMode::Ribbon,
        ColorScheme::Cpk,
        true,
        60.0,
        termpdb::render::Visibility::default(),
        termpdb::render::LodMode::Auto,
        10,
    );
    ratatui::widgets::Widget::render(widget, area, &mut buffer);

    let text: String = (0..80).map(|x| buffer[(x, 0)].symbol()).collect();
    assert!(text.contains("Ribbon"));
    assert!(text.contains("CPK"));
    assert!(text.contains("ON") || text.contains("Spin"));
}

#[test]
fn test_help_widget_render() {
    let area = Rect::new(0, 0, 80, 24);
    let mut buffer = Buffer::empty(area);

    let widget = HelpWidget::new();
    ratatui::widgets::Widget::render(widget, area, &mut buffer);

    let mut found_help_text = false;
    for y in 0..24 {
        let line: String = (0..80).map(|x| buffer[(x, y)].symbol()).collect();
        if line.contains("Help") || line.contains("Controls") || line.contains("Keybindings") {
            found_help_text = true;
            break;
        }
    }
    assert!(found_help_text, "Help modal should render title/controls");
}

#[test]
fn test_info_widget_render() {
    let structure = create_test_structure();
    let area = Rect::new(0, 0, 80, 24);
    let mut buffer = Buffer::empty(area);

    let widget = InfoWidget::new(&structure);
    ratatui::widgets::Widget::render(widget, area, &mut buffer);

    let mut found_info = false;
    for y in 0..24 {
        let line: String = (0..80).map(|x| buffer[(x, y)].symbol()).collect();
        if line.contains("1CRN") || line.contains("Residues") || line.contains("Chains") {
            found_info = true;
            break;
        }
    }
    assert!(found_info, "Info modal should render structure stats");
}

#[test]
fn test_app_render_ui_integration() {
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::Ribbon, ColorScheme::Cpk, false);

    let backend = TestBackend::new(80, 24);
    let mut terminal = Terminal::new(backend).unwrap();

    let res = terminal.draw(|f| {
        app.render_ui(f);
    });
    assert!(res.is_ok());

    // Toggle help and render
    app.toggle_help();
    let res_help = terminal.draw(|f| {
        app.render_ui(f);
    });
    assert!(res_help.is_ok());

    // Toggle info and render
    app.toggle_help();
    app.toggle_info();
    let res_info = terminal.draw(|f| {
        app.render_ui(f);
    });
    assert!(res_info.is_ok());
}

#[test]
fn test_app_model_step_keeps_camera() {
    let pdb = r#"MODEL        1
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ENDMDL
MODEL        2
ATOM      1  CA  ALA A   1       8.000   0.000   0.000  1.00 10.00           C
ENDMDL
END
"#;
    let structure = termpdb::parser::parse_pdb(pdb).unwrap();
    let mut app = App::new(structure, RenderMode::Trace, ColorScheme::Cpk, false);

    app.camera.orbit(12.0, -3.0);
    app.camera.zoom(1.0);
    let orientation = app.camera.orientation;
    let distance = app.camera.distance;
    let target = app.camera.target;

    app.apply_action(AppAction::NextModel);
    assert_eq!(app.structure.active_model_serial(), 2);
    assert!((app.structure.atoms()[0].pos.x - 8.0).abs() < 1e-4);
    assert_eq!(app.camera.orientation, orientation);
    assert_eq!(app.camera.distance, distance);
    assert_eq!(app.camera.target, target);

    app.apply_action(AppAction::PrevModel);
    assert_eq!(app.structure.active_model_serial(), 1);
    assert_eq!(app.camera.orientation, orientation);
}

#[test]
fn test_header_shows_model_only_when_multiple() {
    let single = create_test_structure();
    let area = Rect::new(0, 0, 100, 1);
    let mut buffer = Buffer::empty(area);
    ratatui::widgets::Widget::render(HeaderWidget::new(&single), area, &mut buffer);
    let text: String = (0..100).map(|x| buffer[(x, 0)].symbol()).collect();
    assert!(
        !text.contains("Model"),
        "single-model HUD must not show Model: {text}"
    );

    let pdb = r#"MODEL        1
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ENDMDL
MODEL        2
ATOM      1  CA  ALA A   1       8.000   0.000   0.000  1.00 10.00           C
ENDMDL
END
"#;
    let multi = termpdb::parser::parse_pdb(pdb).unwrap();
    let mut buffer = Buffer::empty(area);
    ratatui::widgets::Widget::render(HeaderWidget::new(&multi), area, &mut buffer);
    let text: String = (0..100).map(|x| buffer[(x, 0)].symbol()).collect();
    assert!(text.contains("Model 1/2"), "expected Model 1/2 in {text}");
}

#[test]
fn test_app_visibility_defaults_and_toggles() {
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::Ribbon, ColorScheme::Cpk, false);
    assert!(!app.visibility.show_waters);
    assert!(app.visibility.show_hydrogens);

    app.apply_action(AppAction::ToggleWaters);
    assert!(app.visibility.show_waters);
    app.apply_action(AppAction::ToggleWaters);
    assert!(!app.visibility.show_waters);

    app.apply_action(AppAction::ToggleHydrogens);
    assert!(!app.visibility.show_hydrogens);
    app.apply_action(AppAction::ToggleHydrogens);
    assert!(app.visibility.show_hydrogens);

    assert_eq!(app.lod, termpdb::render::LodMode::Auto);
    app.apply_action(AppAction::NextLod);
    assert_eq!(app.lod, termpdb::render::LodMode::Full);
    app.apply_action(AppAction::PrevLod);
    assert_eq!(app.lod, termpdb::render::LodMode::Auto);
}

#[test]
fn test_app_selection_and_pick_prompt() {
    let mut app = App::new(
        create_test_structure(),
        RenderMode::Ribbon,
        ColorScheme::Cpk,
        false,
    );
    assert!(app.selection.is_empty());

    app.apply_action(AppAction::PickAtom(0));
    app.apply_action(AppAction::PickAtom(1));
    let line = app.selection.status_line(&app.structure).unwrap();
    assert!(line.contains('Å'), "{line}");

    app.apply_action(AppAction::ClearSelection);
    assert!(app.selection.is_empty());

    let make_key = |code: KeyCode| KeyEvent {
        code,
        modifiers: KeyModifiers::NONE,
        kind: KeyEventKind::Press,
        state: KeyEventState::NONE,
    };
    app.handle_key(make_key(KeyCode::Char('/')));
    assert_eq!(app.pick_prompt.as_deref(), Some(""));
    app.handle_key(make_key(KeyCode::Char('A')));
    app.handle_key(make_key(KeyCode::Char(':')));
    app.handle_key(make_key(KeyCode::Char('1')));
    app.handle_key(make_key(KeyCode::Enter));
    assert!(app.pick_prompt.is_none());
    assert_eq!(app.selection.atoms().len(), 1);
    assert_eq!(app.structure.atoms()[app.selection.atoms()[0]].res_seq, 1);

    app.handle_key(make_key(KeyCode::Char('/')));
    app.handle_key(make_key(KeyCode::Esc));
    assert!(app.pick_prompt.is_none());
}

#[test]
fn test_app_assembly_cycle_clears_selection_and_fits() {
    let pdb = r#"REMARK 350 BIOMOLECULE: 1
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000       10.00000
REMARK 350   BIOMT2   2  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
END
"#;
    let structure = termpdb::parser::parse_pdb(pdb).unwrap();
    let mut app = App::new(structure, RenderMode::Trace, ColorScheme::Cpk, false);
    app.apply_action(AppAction::PickAtom(0));
    assert!(!app.selection.is_empty());
    let asu_distance = app.camera.distance;
    app.apply_action(AppAction::NextAssembly);
    assert_eq!(app.structure.active_assembly_id(), Some("1"));
    assert!(app.selection.is_empty());
    assert_eq!(app.structure.atom_count(), 2);
    assert!(app.camera.distance >= asu_distance);
}

#[test]
fn test_switching_into_ribbon_mode_renders_geometry() {
    // Regression: mode changes never flagged the render cache dirty, so the
    // cached ribbon geometry stayed empty and the viewport rendered blank
    // after switching into Ribbon mode mid-session.
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::Trace, ColorScheme::Cpk, false);
    app.render_scene(80, 48);
    let trace_lit = lit_pixel_count(&app);
    assert!(trace_lit > 0, "trace baseline frame must not be blank");

    app.apply_action(AppAction::SetRenderMode(RenderMode::Ribbon));
    assert_eq!(app.render_mode, RenderMode::Ribbon);
    app.render_scene(80, 48);

    let ribbon_lit = lit_pixel_count(&app);
    assert!(
        ribbon_lit > 100,
        "ribbon viewport blank after direct mode switch: {ribbon_lit} lit pixels"
    );
}

#[test]
fn test_cycling_into_ribbon_mode_renders_geometry() {
    // Same regression via the cycle-key path: BallAndStick -> next == Ribbon.
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::BallAndStick, ColorScheme::Cpk, false);
    app.apply_action(AppAction::NextRenderMode);
    assert_eq!(app.render_mode, RenderMode::Ribbon);
    app.render_scene(80, 48);

    let ribbon_lit = lit_pixel_count(&app);
    assert!(
        ribbon_lit > 100,
        "ribbon viewport blank after cycling into it: {ribbon_lit} lit pixels"
    );
}

#[test]
fn test_interactions_cache_rasterizes_identically_across_frames() {
    // The interaction list is now cached across frames; the second rasterized
    // frame must reuse the cache and produce a pixel-identical image.
    let structure = create_test_structure();
    let mut app = App::new(structure, RenderMode::BallAndStick, ColorScheme::Cpk, false)
        .with_interactions(true);

    app.render_scene(80, 48);
    let first_pass = app.framebuffer.pixels.clone();
    app.render_scene(80, 48);
    assert_eq!(
        first_pass, app.framebuffer.pixels,
        "cached interactions must rasterize identically on reuse"
    );
}

/// Counts non-background pixels in the app's framebuffer.
fn lit_pixel_count(app: &App) -> usize {
    app.framebuffer
        .pixels
        .iter()
        .filter(|&&p| p != (0, 0, 0))
        .count()
}

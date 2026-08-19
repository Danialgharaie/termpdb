//! TUI Event handling and crossterm input mapping.

use crossterm::event::{
    KeyCode, KeyEvent, KeyEventKind, KeyModifiers, MouseButton, MouseEvent, MouseEventKind,
};

use crate::render::RenderMode;

/// Actions generated from user input events.
#[derive(Debug, Clone, PartialEq)]
pub enum AppAction {
    /// Exit the application
    Quit,
    /// Toggle automatic turntable spin
    ToggleSpin,
    /// Increase spin rotation speed
    IncreaseSpinSpeed,
    /// Decrease spin rotation speed
    DecreaseSpinSpeed,
    /// Set a specific molecular representation mode
    SetRenderMode(RenderMode),
    /// Cycle to next representation mode
    NextRenderMode,
    /// Cycle to previous representation mode
    PrevRenderMode,
    /// Cycle to next color scheme
    NextColorScheme,
    /// Cycle to previous color scheme
    PrevColorScheme,
    /// Reset camera view and zoom
    ResetCamera,
    /// Toggle help modal popup
    ToggleHelp,
    /// Toggle structure information popup
    ToggleInfo,
    /// Orbit camera around target by (dx, dy)
    Orbit { dx: f32, dy: f32 },
    /// Pan camera target by (dx, dy)
    Pan { dx: f32, dy: f32 },
    /// Zoom camera by delta
    Zoom { delta: f32 },
    /// No action
    None,
}

/// Mouse interaction tracking state.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MouseState {
    /// Last recorded cursor column and row
    pub last_pos: Option<(u16, u16)>,
    /// Whether left mouse button is currently held down
    pub is_left_down: bool,
    /// Whether right mouse button is currently held down
    pub is_right_down: bool,
}

/// Maps a crossterm `KeyEvent` to an `AppAction`.
pub fn handle_key_event(key: KeyEvent) -> AppAction {
    if key.kind != KeyEventKind::Press {
        return AppAction::None;
    }

    if key.modifiers.contains(KeyModifiers::CONTROL) && key.code == KeyCode::Char('c') {
        return AppAction::Quit;
    }

    match key.code {
        KeyCode::Char('q') | KeyCode::Esc => AppAction::Quit,
        KeyCode::Char(' ') => AppAction::ToggleSpin,
        KeyCode::Char('+') | KeyCode::Char('=') => AppAction::IncreaseSpinSpeed,
        KeyCode::Char('-') | KeyCode::Char('_') => AppAction::DecreaseSpinSpeed,
        KeyCode::Char('1') => AppAction::SetRenderMode(RenderMode::Trace),
        KeyCode::Char('2') => AppAction::SetRenderMode(RenderMode::BallAndStick),
        KeyCode::Char('3') => AppAction::SetRenderMode(RenderMode::Ribbon),
        KeyCode::Char('4') => AppAction::SetRenderMode(RenderMode::Vdw),
        KeyCode::Char('m') => AppAction::NextRenderMode,
        KeyCode::Char('M') => AppAction::PrevRenderMode,
        KeyCode::Char('c') => {
            if key.modifiers.contains(KeyModifiers::SHIFT) {
                AppAction::PrevColorScheme
            } else {
                AppAction::NextColorScheme
            }
        }
        KeyCode::Char('C') => AppAction::PrevColorScheme,
        KeyCode::Char('r') | KeyCode::Char('R') => AppAction::ResetCamera,
        KeyCode::Char('?') | KeyCode::Char('h') | KeyCode::Char('H') => AppAction::ToggleHelp,
        KeyCode::Char('i') | KeyCode::Char('I') => AppAction::ToggleInfo,
        KeyCode::Left => AppAction::Orbit { dx: -5.0, dy: 0.0 },
        KeyCode::Right => AppAction::Orbit { dx: 5.0, dy: 0.0 },
        KeyCode::Up => AppAction::Orbit { dx: 0.0, dy: -5.0 },
        KeyCode::Down => AppAction::Orbit { dx: 0.0, dy: 5.0 },
        KeyCode::Char('w') => AppAction::Pan { dx: 0.0, dy: 1.0 },
        KeyCode::Char('s') => AppAction::Pan { dx: 0.0, dy: -1.0 },
        KeyCode::Char('a') => AppAction::Pan { dx: -1.0, dy: 0.0 },
        KeyCode::Char('d') => AppAction::Pan { dx: 1.0, dy: 0.0 },
        KeyCode::Char('[') => AppAction::Zoom { delta: -1.0 },
        KeyCode::Char(']') => AppAction::Zoom { delta: 1.0 },
        _ => AppAction::None,
    }
}

/// Maps a crossterm `MouseEvent` to an `AppAction`, updating tracking state.
pub fn handle_mouse_event(mouse: MouseEvent, state: &mut MouseState) -> AppAction {
    match mouse.kind {
        MouseEventKind::Down(MouseButton::Left) => {
            state.is_left_down = true;
            state.last_pos = Some((mouse.column, mouse.row));
            AppAction::None
        }
        MouseEventKind::Down(MouseButton::Right) => {
            state.is_right_down = true;
            state.last_pos = Some((mouse.column, mouse.row));
            AppAction::None
        }
        MouseEventKind::Up(MouseButton::Left) => {
            state.is_left_down = false;
            state.last_pos = Some((mouse.column, mouse.row));
            AppAction::None
        }
        MouseEventKind::Up(MouseButton::Right) => {
            state.is_right_down = false;
            state.last_pos = Some((mouse.column, mouse.row));
            AppAction::None
        }
        MouseEventKind::Drag(MouseButton::Left) => {
            let (dx, dy) = if let Some((lx, ly)) = state.last_pos {
                (
                    (mouse.column as f32) - (lx as f32),
                    (mouse.row as f32) - (ly as f32),
                )
            } else {
                (0.0, 0.0)
            };
            state.last_pos = Some((mouse.column, mouse.row));

            if mouse.modifiers.contains(KeyModifiers::SHIFT) {
                AppAction::Pan { dx, dy }
            } else {
                AppAction::Orbit { dx, dy }
            }
        }
        MouseEventKind::Drag(MouseButton::Right) | MouseEventKind::Drag(MouseButton::Middle) => {
            let (dx, dy) = if let Some((lx, ly)) = state.last_pos {
                (
                    (mouse.column as f32) - (lx as f32),
                    (mouse.row as f32) - (ly as f32),
                )
            } else {
                (0.0, 0.0)
            };
            state.last_pos = Some((mouse.column, mouse.row));
            AppAction::Pan { dx, dy }
        }
        MouseEventKind::ScrollUp => AppAction::Zoom { delta: 1.0 },
        MouseEventKind::ScrollDown => AppAction::Zoom { delta: -1.0 },
        MouseEventKind::Moved => {
            state.last_pos = Some((mouse.column, mouse.row));
            AppAction::None
        }
        _ => AppAction::None,
    }
}

//! Main TUI Application state machine and rendering orchestrator.

use std::time::Instant;

use crossterm::event::{KeyEvent, MouseEvent};
use ratatui::Frame;
use ratatui::layout::{Constraint, Direction, Layout};

use crate::model::Structure;
use crate::render::{Camera, ColorScheme, Framebuffer, Lighting, RenderMode, render_structure};
use crate::tui::events::{AppAction, MouseState, handle_key_event, handle_mouse_event};
use crate::tui::widgets::{FooterWidget, HeaderWidget, HelpWidget, InfoWidget, ViewportWidget};

/// Main application state for the interactive 3D terminal molecular viewer.
pub struct App {
    /// Loaded molecular structure
    pub structure: Structure,
    /// 3D orbit camera
    pub camera: Camera,
    /// Software framebuffer and depth buffer
    pub framebuffer: Framebuffer,
    /// Directional lighting and depth cueing
    pub lighting: Lighting,
    /// Active molecular representation mode
    pub render_mode: RenderMode,
    /// Active color scheme
    pub color_scheme: ColorScheme,
    /// Whether automatic turntable spin is active
    pub auto_spin: bool,
    /// Spin rotation speed factor
    pub spin_speed: f32,
    /// Whether help overlay modal is open
    pub show_help: bool,
    /// Whether structure information modal is open
    pub show_info: bool,
    /// Current calculated frames per second
    pub fps: f32,
    /// Total rendered frame count
    pub frame_count: u64,
    /// Timestamp of the last tick
    pub last_tick: Instant,
    /// Timestamp of last FPS counter recalculation
    pub last_fps_update: Instant,
    /// Number of frames rendered since last FPS update
    pub frames_since_fps: u32,
    /// Whether the application should exit
    pub should_quit: bool,
    /// Mouse interaction tracking state
    pub mouse_state: MouseState,
}

impl App {
    /// Creates a new `App` instance initialized with the given structure and options.
    pub fn new(
        structure: Structure,
        initial_mode: RenderMode,
        initial_color: ColorScheme,
        auto_spin: bool,
    ) -> Self {
        let mut camera = Camera::new();
        let com = structure.center_of_mass();
        let radius = structure.bounding_sphere_radius();
        camera.fit_structure(com, radius);

        let framebuffer = Framebuffer::new(80, 48);
        let lighting = Lighting::default();
        let now = Instant::now();

        Self {
            structure,
            camera,
            framebuffer,
            lighting,
            render_mode: initial_mode,
            color_scheme: initial_color,
            auto_spin,
            spin_speed: 1.0,
            show_help: false,
            show_info: false,
            fps: 0.0,
            frame_count: 0,
            last_tick: now,
            last_fps_update: now,
            frames_since_fps: 0,
            should_quit: false,
            mouse_state: MouseState::default(),
        }
    }

    /// Sets the spin speed factor using builder pattern.
    pub fn with_spin_speed(mut self, speed: f32) -> Self {
        self.spin_speed = speed;
        self
    }

    /// Sets the spin speed factor.
    pub fn set_spin_speed(&mut self, speed: f32) {
        self.spin_speed = speed;
    }

    /// Dispatches an `AppAction` to update app state or camera.
    pub fn apply_action(&mut self, action: AppAction) {
        match action {
            AppAction::Quit => {
                if self.show_help {
                    self.show_help = false;
                } else if self.show_info {
                    self.show_info = false;
                } else {
                    self.should_quit = true;
                }
            }
            AppAction::ToggleSpin => self.toggle_spin(),
            AppAction::IncreaseSpinSpeed => self.increase_spin_speed(),
            AppAction::DecreaseSpinSpeed => self.decrease_spin_speed(),
            AppAction::SetRenderMode(mode) => self.set_mode(mode),
            AppAction::NextRenderMode => self.next_mode(),
            AppAction::PrevRenderMode => self.prev_mode(),
            AppAction::NextColorScheme => self.next_color_scheme(),
            AppAction::PrevColorScheme => self.prev_color_scheme(),
            AppAction::ResetCamera => self.reset_camera(),
            AppAction::ToggleHelp => self.toggle_help(),
            AppAction::ToggleInfo => self.toggle_info(),
            AppAction::Orbit { dx, dy } => self.camera.orbit(dx, dy),
            AppAction::Pan { dx, dy } => self.camera.pan(dx * 0.2, dy * 0.2),
            AppAction::Zoom { delta } => self.camera.zoom(delta),
            AppAction::None => {}
        }
    }

    /// Toggles automatic turntable camera spinning.
    pub fn toggle_spin(&mut self) {
        self.auto_spin = !self.auto_spin;
    }

    /// Increases spin speed by 25%.
    pub fn increase_spin_speed(&mut self) {
        self.spin_speed = (self.spin_speed * 1.25).min(10.0);
    }

    /// Decreases spin speed by 20%.
    pub fn decrease_spin_speed(&mut self) {
        self.spin_speed = (self.spin_speed * 0.8).max(0.1);
    }

    /// Sets the rendering representation mode directly.
    pub fn set_mode(&mut self, mode: RenderMode) {
        self.render_mode = mode;
    }

    /// Cycles to the next rendering mode.
    pub fn next_mode(&mut self) {
        self.render_mode = self.render_mode.next();
    }

    /// Cycles to the previous rendering mode.
    pub fn prev_mode(&mut self) {
        self.render_mode = self.render_mode.prev();
    }

    /// Cycles to the next color scheme.
    pub fn next_color_scheme(&mut self) {
        self.color_scheme = self.color_scheme.next();
    }

    /// Cycles to the previous color scheme.
    pub fn prev_color_scheme(&mut self) {
        self.color_scheme = self.color_scheme.prev();
    }

    /// Resets the camera orientation, target, and framing.
    pub fn reset_camera(&mut self) {
        let com = self.structure.center_of_mass();
        let radius = self.structure.bounding_sphere_radius();
        self.camera.fit_structure(com, radius);
    }

    /// Toggles the help overlay modal.
    pub fn toggle_help(&mut self) {
        self.show_help = !self.show_help;
        if self.show_help {
            self.show_info = false;
        }
    }

    /// Toggles the structure info overlay modal.
    pub fn toggle_info(&mut self) {
        self.show_info = !self.show_info;
        if self.show_info {
            self.show_help = false;
        }
    }

    /// Handles a keyboard input event.
    pub fn handle_key(&mut self, key: KeyEvent) {
        let action = handle_key_event(key);
        self.apply_action(action);
    }

    /// Handles a mouse input event.
    pub fn handle_mouse(&mut self, mouse: MouseEvent) {
        let action = handle_mouse_event(mouse, &mut self.mouse_state);
        self.apply_action(action);
    }

    /// Advances simulation/animation time by `delta_time` seconds.
    pub fn tick(&mut self, delta_time: f32) {
        if self.auto_spin {
            self.camera.orbit(self.spin_speed * delta_time * 30.0, 0.0);
        }

        self.frame_count = self.frame_count.wrapping_add(1);
        self.frames_since_fps += 1;

        let elapsed = self.last_fps_update.elapsed().as_secs_f32();
        if elapsed >= 0.5 {
            self.fps = (self.frames_since_fps as f32) / elapsed;
            self.frames_since_fps = 0;
            self.last_fps_update = Instant::now();
        }
    }

    /// Re-renders the 3D scene into the software framebuffer matching viewport dimensions.
    pub fn render_scene(&mut self, width: usize, height: usize) {
        let pixel_width = width;
        let pixel_height = height * 2;

        if pixel_width == 0 || pixel_height == 0 {
            return;
        }

        if self.framebuffer.width != pixel_width || self.framebuffer.height != pixel_height {
            self.framebuffer.resize(pixel_width, pixel_height);
        }

        self.camera.aspect = (pixel_width as f32) / (pixel_height as f32);
        self.framebuffer.clear((0, 0, 0));

        render_structure(
            &self.structure,
            self.render_mode,
            self.color_scheme,
            &self.camera,
            &mut self.framebuffer,
            &self.lighting,
        );
    }

    /// Renders the complete UI layout (Header, Viewport, Footer, Modals) into the Ratatui frame.
    pub fn render_ui(&mut self, frame: &mut Frame) {
        let area = frame.area();
        if area.width == 0 || area.height == 0 {
            return;
        }

        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Min(1),
                Constraint::Length(1),
            ])
            .split(area);

        // Header
        frame.render_widget(HeaderWidget::new(&self.structure), chunks[0]);

        // 3D Viewport
        let v_width = chunks[1].width as usize;
        let v_height = chunks[1].height as usize;
        self.render_scene(v_width, v_height);
        frame.render_widget(ViewportWidget::new(&self.framebuffer), chunks[1]);

        // Footer
        frame.render_widget(
            FooterWidget::new(
                self.render_mode,
                self.color_scheme,
                self.auto_spin,
                self.fps,
            ),
            chunks[2],
        );

        // Overlays
        if self.show_help {
            frame.render_widget(HelpWidget::new(), area);
        } else if self.show_info {
            frame.render_widget(InfoWidget::new(&self.structure), area);
        }
    }
}

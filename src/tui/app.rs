//! Main TUI Application state machine and rendering orchestrator.

use std::time::Instant;

use crossterm::event::{KeyCode, KeyEvent, KeyEventKind, KeyModifiers, MouseEvent};
use ratatui::Frame;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::widgets::Paragraph;

use crate::math::Vec3;
use crate::model::{Interaction, Structure};
use crate::render::{
    Camera, ColorScheme, Framebuffer, GraphicsBackend, Lighting, LodMode, PixelColor,
    RenderContext, RenderMode, RibbonPrimitive, Visibility, build_render_cache,
    build_ribbon_geometry, draw_selection_markers, render_structure_ctx,
};
use crate::select::{Selection, parse_atom_spec, pick_atom_at_screen, resolve_atom};
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
    /// Active graphics backend
    pub graphics_backend: GraphicsBackend,
    /// Resolution scale multiplier for Kitty graphics mode
    pub scale: f32,
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
    /// Atom visibility filter (waters / hydrogens)
    pub visibility: Visibility,
    /// Up to four selected atoms in the active model
    pub selection: Selection,
    /// Open `/` pick prompt buffer, if any
    pub pick_prompt: Option<String>,
    /// Last pick-prompt error to show beside the query
    pub pick_error: Option<String>,
    /// Last 3D viewport rect (for mapping mouse clicks to framebuffer pixels)
    pub viewport_area: Rect,
    /// Level-of-detail for large structures
    pub lod: LodMode,
    /// Postprocessing configuration (outlines and SSAO)
    pub postprocess_config: crate::render::PostProcessConfig,
    /// Whether to render non-covalent interaction meshes (H-bonds & disulfide bridges)
    pub show_interactions: bool,
    /// Whether the 3D viewport must be re-rasterized on the next draw. Set by
    /// any state change that affects the scene (camera, mode, colors, model,
    /// assembly, LOD, visibility, selection); cleared once rasterized.
    pub needs_rerender: bool,
    /// Whether the full UI must be re-painted (ratatui draw). Broader than
    /// needs_rerender: set by any input, spin, FPS readout tick, or resize so
    /// the TUI does NOT redraw 60 fps while idle (CPU ~0 at rest).
    pub needs_redraw: bool,
    /// Per-atom render cache (colors, visibility flags, bounding sphere) rebuilt
    /// only when the structure / color scheme / visibility / LOD changes -- not
    /// when the camera moves, so orbit/spin reuse it across frames.
    render_cache: RenderCache,
    /// True when render_cache must be rebuilt before the next scene render.
    render_cache_dirty: bool,
    /// Cached, camera-independent ribbon geometry (spline + ligand primitives),
    /// rebuilt alongside render_cache when structure/color/visibility/LOD changes.
    ribbon_geometry: Vec<RibbonPrimitive>,
    /// Detected non-covalent interactions (H-bonds, disulfides) for the active
    /// structure view. Detection is O(N²), so it is cached and only recomputed
    /// after a model or assembly switch -- never per frame.
    cached_interactions: Vec<Interaction>,
    /// True when `cached_interactions` must be recomputed before drawing.
    interactions_dirty: bool,
}

/// Camera-independent per-atom data cached across frames during orbit/spin.
#[derive(Clone)]
struct RenderCache {
    colors: Vec<PixelColor>,
    visible: Vec<bool>,
    com: Vec3,
    radius: f32,
    max_vdw: f32,
}

impl Default for RenderCache {
    fn default() -> Self {
        Self {
            colors: Vec::new(),
            visible: Vec::new(),
            com: Vec3::ZERO,
            radius: 1.0,
            max_vdw: 1.5,
        }
    }
}

impl App {
    /// Creates a new `App` instance initialized with the given structure and options.
    pub fn new(
        mut structure: Structure,
        initial_mode: RenderMode,
        initial_color: ColorScheme,
        auto_spin: bool,
    ) -> Self {
        structure.ensure_bonds();
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
            graphics_backend: GraphicsBackend::HalfBlock,
            scale: 1.0,
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
            visibility: Visibility::default(),
            selection: Selection::default(),
            pick_prompt: None,
            pick_error: None,
            viewport_area: Rect::default(),
            lod: LodMode::Auto,
            postprocess_config: crate::render::PostProcessConfig::default(),
            show_interactions: false,
            needs_rerender: true,
            needs_redraw: true,
            render_cache: RenderCache::default(),
            render_cache_dirty: true,
            ribbon_geometry: Vec::new(),
            cached_interactions: Vec::new(),
            interactions_dirty: true,
        }
    }

    /// Flags the camera-independent render cache (per-atom colors, visibility,
    /// bounding sphere, ribbon geometry) as stale so it is rebuilt before the
    /// next scene raster. Required after any change to structure content,
    /// color scheme, visibility, LOD, **or render mode** -- ribbon geometry is
    /// cached alongside the atom colors and only exists in Ribbon mode.
    fn mark_render_cache_dirty(&mut self) {
        self.render_cache_dirty = true;
    }

    /// Flags every structure-derived cache (render cache + detected
    /// interactions) as stale. Required after model or assembly switches,
    /// which replace the active atom view.
    fn mark_scene_dirty(&mut self) {
        self.render_cache_dirty = true;
        self.interactions_dirty = true;
    }

    /// Sets the level-of-detail mode.
    pub fn with_lod(mut self, lod: LodMode) -> Self {
        self.lod = lod;
        self.mark_render_cache_dirty();
        self
    }

    /// Sets the atom visibility filter.
    pub fn with_visibility(mut self, visibility: Visibility) -> Self {
        self.visibility = visibility;
        self.mark_render_cache_dirty();
        self
    }

    /// Sets the spin speed factor using builder pattern.
    pub fn with_spin_speed(mut self, speed: f32) -> Self {
        self.spin_speed = speed;
        self
    }

    /// Sets the post-processing configuration (outlines and SSAO).
    pub fn with_postprocess(mut self, config: crate::render::PostProcessConfig) -> Self {
        self.postprocess_config = config;
        self
    }

    /// Sets whether non-covalent interaction meshes are rendered.
    pub fn with_interactions(mut self, show: bool) -> Self {
        self.show_interactions = show;
        self
    }

    /// Sets the Depth-of-Field focus distance.
    pub fn with_dof(mut self, focus: Option<f32>) -> Self {
        self.lighting.dof_focus = focus;
        self
    }

    /// Sets the active graphics backend.
    pub fn with_graphics_backend(mut self, backend: GraphicsBackend) -> Self {
        self.graphics_backend = backend;
        self
    }

    /// Sets the Kitty graphics resolution scale factor.
    pub fn with_scale(mut self, scale: f32) -> Self {
        self.scale = scale;
        self
    }

    /// Sets the spin speed factor.
    pub fn set_spin_speed(&mut self, speed: f32) {
        self.spin_speed = speed;
    }

    /// Sets the active graphics backend.
    pub fn set_graphics_backend(&mut self, backend: GraphicsBackend) {
        if self.graphics_backend != backend {
            let was_kitty = self.graphics_backend.is_kitty();
            self.graphics_backend = backend;
            if was_kitty && !self.graphics_backend.is_kitty() {
                use std::io::Write;
                let delete_seq = crate::render::encode_kitty_delete(None);
                let _ = std::io::stdout().write_all(delete_seq.as_bytes());
                let _ = std::io::stdout().flush();
            }
            self.needs_rerender = true;
            self.needs_redraw = true;
        }
    }

    /// Toggles between HalfBlock and Kitty Graphics Protocol rendering backends.
    pub fn toggle_graphics_backend(&mut self) {
        let was_kitty = self.graphics_backend.is_kitty();
        self.graphics_backend.toggle();
        if was_kitty && !self.graphics_backend.is_kitty() {
            use std::io::Write;
            let delete_seq = crate::render::encode_kitty_delete(None);
            let _ = std::io::stdout().write_all(delete_seq.as_bytes());
            let _ = std::io::stdout().flush();
        }
        self.needs_rerender = true;
        self.needs_redraw = true;
    }

    /// Resizes the framebuffer to match the specified terminal column and row dimensions.
    pub fn resize_framebuffer(&mut self, width: u16, height: u16) {
        if self.graphics_backend.is_kitty() {
            let (cell_w, cell_h) = crate::render::get_terminal_cell_size_scaled(self.scale);
            let pixel_w = (width as u32 * cell_w).max(1) as usize;
            let pixel_h = (height as u32 * cell_h).max(1) as usize;
            self.framebuffer.resize(pixel_w, pixel_h);
        } else {
            self.framebuffer
                .resize(width as usize, (height * 2) as usize);
        }
    }

    /// Emits the in-band Kitty graphics escape sequence for the current framebuffer placed over viewport_area.
    pub fn emit_kitty_frame(&self) {
        if self.viewport_area.width == 0
            || self.viewport_area.height == 0
            || self.show_help
            || self.show_info
        {
            return;
        }
        let rgba = self.framebuffer.to_rgba_bytes();
        let seq = crate::render::encode_kitty_graphics_png(
            self.framebuffer.width as u32,
            self.framebuffer.height as u32,
            self.viewport_area.width,
            self.viewport_area.height,
            self.viewport_area.x,
            self.viewport_area.y,
            0,
            1,
            &rgba,
        );
        use std::io::Write;
        let mut out = std::io::stdout().lock();
        let _ = out.write_all(seq.as_bytes());
        let _ = out.flush();
    }

    /// Dispatches an `AppAction` to update app state or camera.
    pub fn apply_action(&mut self, action: AppAction) {
        // Any non-trivial action requires a UI redraw (overlay / mode / scene / prompt).
        if !matches!(action, AppAction::None) {
            self.needs_redraw = true;
        }
        // Actions that change the 3D scene additionally require re-rasterization.
        // Modal / spin / quit / prompt actions only touch overlays or animation
        // state, so they are excluded to avoid wasted renders while idle.
        if !matches!(
            action,
            AppAction::None
                | AppAction::Quit
                | AppAction::ToggleSpin
                | AppAction::IncreaseSpinSpeed
                | AppAction::DecreaseSpinSpeed
                | AppAction::ToggleHelp
                | AppAction::ToggleInfo
                | AppAction::StartPickPrompt
        ) {
            self.needs_rerender = true;
        }
        match action {
            AppAction::Quit => {
                if self.pick_prompt.is_some() {
                    self.pick_prompt = None;
                    self.pick_error = None;
                } else if self.show_help {
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
            AppAction::NextColorScheme => {
                self.next_color_scheme();
                self.mark_render_cache_dirty();
            }
            AppAction::PrevColorScheme => {
                self.prev_color_scheme();
                self.mark_render_cache_dirty();
            }
            AppAction::ResetCamera => self.reset_camera(),
            AppAction::ToggleHelp => self.toggle_help(),
            AppAction::ToggleInfo => self.toggle_info(),
            AppAction::NextModel => {
                self.structure.next_model();
                self.structure.ensure_bonds();
                self.selection.clear();
                self.mark_scene_dirty();
            }
            AppAction::PrevModel => {
                self.structure.prev_model();
                self.structure.ensure_bonds();
                self.selection.clear();
                self.mark_scene_dirty();
            }
            AppAction::NextAssembly => {
                self.structure.next_assembly();
                self.structure.ensure_bonds();
                self.selection.clear();
                self.reset_camera();
                self.mark_scene_dirty();
            }
            AppAction::PrevAssembly => {
                self.structure.prev_assembly();
                self.structure.ensure_bonds();
                self.selection.clear();
                self.reset_camera();
                self.mark_scene_dirty();
            }
            AppAction::NextLod => {
                self.lod = self.lod.next();
                self.mark_render_cache_dirty();
            }
            AppAction::PrevLod => {
                self.lod = self.lod.prev();
                self.mark_render_cache_dirty();
            }
            AppAction::StartPickPrompt => {
                self.pick_prompt = Some(String::new());
                self.pick_error = None;
                self.show_help = false;
                self.show_info = false;
            }
            AppAction::ClearSelection => self.selection.clear(),
            AppAction::PickAtom(idx) => self.selection.pick(idx),
            AppAction::PickAt { col, row } => self.pick_at_cell(col, row),
            AppAction::ToggleWaters => {
                self.visibility.show_waters = !self.visibility.show_waters;
                self.mark_render_cache_dirty();
            }
            AppAction::ToggleHydrogens => {
                self.visibility.show_hydrogens = !self.visibility.show_hydrogens;
                self.mark_render_cache_dirty();
            }
            AppAction::ToggleOutline => {
                self.postprocess_config.outline = !self.postprocess_config.outline;
            }
            AppAction::ToggleSsao => {
                self.postprocess_config.ssao = !self.postprocess_config.ssao;
            }
            AppAction::ToggleGraphicsBackend => self.toggle_graphics_backend(),
            AppAction::ToggleInteractions => {
                self.show_interactions = !self.show_interactions;
            }
            AppAction::ToggleDof => {
                if self.lighting.dof_focus.is_some() {
                    self.lighting.dof_focus = None;
                } else {
                    self.lighting.dof_focus = Some(self.camera.distance);
                }
            }
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
    ///
    /// Ribbon geometry lives in the render cache (it is only built while
    /// [`RenderMode::Ribbon`] is active), so every mode change must flag the
    /// cache stale -- otherwise switching into Ribbon mode renders an empty
    /// viewport.
    pub fn set_mode(&mut self, mode: RenderMode) {
        self.render_mode = mode;
        self.mark_render_cache_dirty();
    }

    /// Cycles to the next rendering mode.
    pub fn next_mode(&mut self) {
        self.set_mode(self.render_mode.next());
    }

    /// Cycles to the previous rendering mode.
    pub fn prev_mode(&mut self) {
        self.set_mode(self.render_mode.prev());
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
        if self.pick_prompt.is_some() {
            self.handle_pick_prompt_key(key);
            return;
        }
        let action = handle_key_event(key);
        self.apply_action(action);
    }

    fn handle_pick_prompt_key(&mut self, key: KeyEvent) {
        if key.kind != KeyEventKind::Press {
            return;
        }
        // Any key in the pick prompt updates the prompt line / selection / error.
        self.needs_redraw = true;
        if key.modifiers.contains(KeyModifiers::CONTROL) && key.code == KeyCode::Char('c') {
            self.should_quit = true;
            return;
        }
        match key.code {
            KeyCode::Esc => {
                self.pick_prompt = None;
                self.pick_error = None;
            }
            KeyCode::Enter => self.submit_pick_prompt(),
            KeyCode::Backspace => {
                if let Some(buf) = &mut self.pick_prompt {
                    buf.pop();
                }
                self.pick_error = None;
            }
            KeyCode::Char(c) => {
                if let Some(buf) = &mut self.pick_prompt {
                    buf.push(c);
                }
                self.pick_error = None;
            }
            _ => {}
        }
    }

    fn submit_pick_prompt(&mut self) {
        let Some(query) = self.pick_prompt.clone() else {
            return;
        };
        let trimmed = query.trim();
        if trimmed.is_empty() {
            self.pick_prompt = None;
            return;
        }
        match parse_atom_spec(trimmed)
            .and_then(|spec| resolve_atom(&self.structure, &spec, Some(&self.visibility)))
        {
            Ok(idx) => {
                self.selection.pick(idx);
                self.pick_prompt = None;
                self.pick_error = None;
                self.needs_rerender = true;
                self.needs_redraw = true;
            }
            Err(err) => {
                self.pick_error = Some(err.to_string());
                self.needs_redraw = true;
            }
        }
    }

    fn pick_at_cell(&mut self, col: u16, row: u16) {
        let area = self.viewport_area;
        if area.width == 0 || area.height == 0 {
            return;
        }
        if col < area.x || row < area.y {
            return;
        }
        if col >= area.x + area.width || row >= area.y + area.height {
            return;
        }
        let (sx, sy) = if self.graphics_backend.is_kitty() {
            let (cell_w, cell_h) = crate::render::get_terminal_cell_size_scaled(self.scale);
            let px = (col.saturating_sub(area.x) as u32 * cell_w + cell_w / 2) as f32;
            let py = (row.saturating_sub(area.y) as u32 * cell_h + cell_h / 2) as f32;
            (px, py)
        } else {
            let sx = (col - area.x) as f32 + 0.5;
            let sy = ((row - area.y) as f32) * 2.0 + 1.0;
            (sx, sy)
        };
        let pick_radius = if self.graphics_backend.is_kitty() {
            let (cell_w, cell_h) = crate::render::get_terminal_cell_size_scaled(self.scale);
            (cell_w.max(cell_h) as f32) * 1.5
        } else {
            6.0
        };
        if let Some(idx) = pick_atom_at_screen(
            &self.structure,
            &self.camera,
            self.visibility,
            self.framebuffer.width,
            self.framebuffer.height,
            sx,
            sy,
            pick_radius,
        ) {
            self.selection.pick(idx);
        }
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
            // The camera moved this frame, so the viewport must be re-rasterized.
            self.needs_rerender = true;
            self.needs_redraw = true;
        }

        self.frame_count = self.frame_count.wrapping_add(1);
        self.frames_since_fps += 1;

        let elapsed = self.last_fps_update.elapsed().as_secs_f32();
        if elapsed >= 0.5 {
            self.fps = (self.frames_since_fps as f32) / elapsed;
            self.frames_since_fps = 0;
            self.last_fps_update = Instant::now();
            // Footer FPS readout changed -> repaint (but not re-rasterize).
            self.needs_redraw = true;
        }
    }

    /// Re-renders the 3D scene into the software framebuffer matching viewport dimensions.
    pub fn render_scene(&mut self, width: usize, height: usize) {
        let (pixel_width, pixel_height) = if self.graphics_backend.is_kitty() {
            let (cell_w, cell_h) = crate::render::get_terminal_cell_size_scaled(self.scale);
            (
                (width as u32 * cell_w).max(1) as usize,
                (height as u32 * cell_h).max(1) as usize,
            )
        } else {
            (width, height * 2)
        };

        if pixel_width == 0 || pixel_height == 0 {
            return;
        }

        if self.framebuffer.width != pixel_width || self.framebuffer.height != pixel_height {
            self.framebuffer.resize(pixel_width, pixel_height);
        }

        self.camera.aspect = (pixel_width as f32) / (pixel_height as f32);
        self.framebuffer.clear((0, 0, 0));

        // Rebuild the camera-independent per-atom cache only when the structure,
        // color scheme, visibility, or LOD changed -- not every frame.
        self.ensure_render_cache();

        // Refresh the interaction cache before `ctx` borrows scene fields.
        if self.show_interactions {
            self.ensure_interactions();
        }

        let level = self.lod.resolve(self.structure.atom_count());
        let ctx = RenderContext {
            structure: &self.structure,
            camera: &self.camera,
            mats: self.camera.matrices(),
            lighting: &self.lighting,
            visibility: self.visibility,
            lod: level,
            colors: &self.render_cache.colors,
            visible: &self.render_cache.visible,
            com: self.render_cache.com,
            radius: self.render_cache.radius,
            max_vdw: self.render_cache.max_vdw,
            ribbon_geometry: Some(&self.ribbon_geometry),
        };
        render_structure_ctx(&ctx, self.render_mode, &mut self.framebuffer);

        if self.show_interactions {
            let atoms = self.structure.atoms();
            for inter in &self.cached_interactions {
                if inter.atom1_idx < atoms.len() && inter.atom2_idx < atoms.len() {
                    let a1 = &atoms[inter.atom1_idx];
                    let a2 = &atoms[inter.atom2_idx];
                    if let (Some(p1), Some(p2)) = (
                        self.camera.project(
                            &ctx.mats,
                            a1.pos,
                            self.framebuffer.width,
                            self.framebuffer.height,
                        ),
                        self.camera.project(
                            &ctx.mats,
                            a2.pos,
                            self.framebuffer.width,
                            self.framebuffer.height,
                        ),
                    ) {
                        let color = inter.kind.default_color();
                        crate::render::rasterizer::draw_dashed_line_3d(
                            &mut self.framebuffer,
                            p1,
                            p2,
                            color,
                            4.0,
                            2.0,
                        );
                    }
                }
            }
        }

        draw_selection_markers(
            &self.structure,
            &self.camera,
            &mut self.framebuffer,
            self.selection.atoms(),
        );

        crate::render::postprocess::apply_postprocessing(
            &mut self.framebuffer,
            &self.postprocess_config,
        );
    }

    /// Lazily recomputes the non-covalent interaction list after a model or
    /// assembly switch. Detection is O(N²) over all atoms, so the result is
    /// cached and reused across frames; recomputing every rasterized frame
    /// made `--interactions` sessions crawl on large structures.
    fn ensure_interactions(&mut self) {
        if !self.interactions_dirty {
            return;
        }
        self.cached_interactions = crate::model::detect_interactions(&self.structure);
        self.interactions_dirty = false;
    }

    /// Rebuilds the per-atom render cache (colors, visibility flags, bounding
    /// sphere) if the structure, color scheme, visibility, or LOD changed since
    /// the last build. Cheap (returns immediately) when the cache is still valid,
    /// which is the common case during orbit/zoom/spin.
    fn ensure_render_cache(&mut self) {
        if !self.render_cache_dirty {
            return;
        }
        let (colors, visible, com, radius, max_vdw) = build_render_cache(
            &self.structure,
            self.color_scheme,
            self.visibility,
            self.lod,
        );
        let level = self.lod.resolve(self.structure.atom_count());
        let ribbon_geometry = if self.render_mode == RenderMode::Ribbon {
            build_ribbon_geometry(&self.structure, &colors, &visible, self.visibility, level)
        } else {
            Vec::new()
        };
        self.render_cache = RenderCache {
            colors,
            visible,
            com,
            radius,
            max_vdw,
        };
        self.ribbon_geometry = ribbon_geometry;
        self.render_cache_dirty = false;
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
        let sel_line = self.selection.status_line(&self.structure);
        let header = if let Some(ref s) = sel_line {
            HeaderWidget::new(&self.structure).with_selection(s)
        } else {
            HeaderWidget::new(&self.structure)
        };
        frame.render_widget(header, chunks[0]);

        // 3D Viewport
        self.viewport_area = chunks[1];
        let v_width = chunks[1].width as usize;
        let v_height = chunks[1].height as usize;
        let expected_size = if self.graphics_backend.is_kitty() {
            let (cell_w, cell_h) = crate::render::get_terminal_cell_size_scaled(self.scale);
            (
                (v_width as u32 * cell_w).max(1) as usize,
                (v_height as u32 * cell_h).max(1) as usize,
            )
        } else {
            (v_width, v_height * 2)
        };

        // Only re-rasterize the 3D scene when something changed (or the viewport
        // was resized); an unchanged framebuffer is painted straight from cache,
        // so the event loop stays responsive while idle.
        if self.needs_rerender
            || self.framebuffer.width != expected_size.0
            || self.framebuffer.height != expected_size.1
        {
            self.render_scene(v_width, v_height);
            self.needs_rerender = false;
        }
        frame.render_widget(
            ViewportWidget::new(&self.framebuffer).with_backend(self.graphics_backend),
            chunks[1],
        );

        // Footer or pick prompt
        if let Some(query) = &self.pick_prompt {
            let err = self
                .pick_error
                .as_deref()
                .map(|e| format!("  [{e}]"))
                .unwrap_or_default();
            let line = format!(" / {query}_{err}");
            frame.render_widget(
                Paragraph::new(line).style(
                    Style::default()
                        .fg(Color::Yellow)
                        .bg(Color::Black)
                        .add_modifier(Modifier::BOLD),
                ),
                chunks[2],
            );
        } else {
            frame.render_widget(
                FooterWidget::new(
                    self.render_mode,
                    self.color_scheme,
                    self.auto_spin,
                    self.fps,
                    self.visibility,
                    self.lod,
                    self.structure.atom_count(),
                ),
                chunks[2],
            );
        }

        // Overlays
        if self.show_help {
            frame.render_widget(HelpWidget::new(), area);
        } else if self.show_info {
            frame.render_widget(InfoWidget::new(&self.structure), area);
        }
    }
}

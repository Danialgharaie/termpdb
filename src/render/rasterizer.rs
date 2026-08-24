//! 3D software rasterization primitives.
//!
//! Provides depth-buffered rasterization for spheres (analytical ray-casting),
//! 3D lines (Bresenham/DDA with depth interpolation), thick cylinders/capsules,
//! and 3D shaded triangles with barycentric interpolation.

use crate::math::Vec3;
use crate::render::buffer::{Framebuffer, PixelColor};
use crate::render::lighting::Lighting;

/// Rasterizes an analytical 3D sphere with Lambertian shading and depth-buffer testing.
///
/// `center_screen`: `(screen_x, screen_y, view_depth_z)`
/// `radius_screen`: radius in screen pixels
pub fn draw_sphere(
    buffer: &mut Framebuffer,
    center_screen: (f32, f32, f32),
    radius_screen: f32,
    base_color: PixelColor,
    lighting: &Lighting,
) {
    if radius_screen < 0.35 {
        return;
    }

    let (cx, cy, cz) = center_screen;
    let r = radius_screen;

    if cx + r < 0.0
        || cx - r >= buffer.width as f32
        || cy + r < 0.0
        || cy - r >= buffer.height as f32
    {
        return;
    }

    let r_sq = r * r;

    let min_x = ((cx - r).floor() as i32).max(0);
    let max_x = ((cx + r).ceil() as i32).min((buffer.width as i32) - 1);
    let min_y = ((cy - r).floor() as i32).max(0);
    let max_y = ((cy + r).ceil() as i32).min((buffer.height as i32) - 1);

    if min_x > max_x || min_y > max_y {
        return;
    }

    let inv_r = 1.0 / r;
    let min_depth = cz - r;
    let max_depth = cz + r;
    // Nearest possible z of this sphere (at its center pixel). Cheap early
    // reject: a pixel whose stored depth is already <= z_near cannot be
    // overwritten by this sphere, so we skip the sqrt + shading for it -- the
    // key win where space-filling spheres overlap heavily (VDW).
    let z_near = min_depth;

    let width = buffer.width;
    let depth = buffer.depth.as_mut_ptr();
    let pixels = buffer.pixels.as_mut_ptr();

    for y in min_y..=max_y {
        let dy = (y as f32 + 0.5) - cy;
        let dy_sq = dy * dy;
        let row = (y as usize) * width;
        for x in min_x..=max_x {
            // SAFETY: x in [0, width-1] and y in [0, height-1] (clamped above),
            // so idx = y*width + x is in [0, width*height-1] = the buffer length.
            let idx = row + x as usize;
            let depth_cur = unsafe { *depth.add(idx) };
            if z_near < depth_cur {
                let dx = (x as f32 + 0.5) - cx;
                let dist_sq = dx * dx + dy_sq;
                if dist_sq <= r_sq {
                    let dz = (r_sq - dist_sq).sqrt();
                    let z = cz - dz;
                    if z < depth_cur {
                        // Normal in camera view space (+X right, +Y up, +Z towards camera)
                        let normal = Vec3::new(dx * inv_r, -dy * inv_r, dz * inv_r);
                        let lit_color = lighting.shade(normal, z, base_color, min_depth, max_depth);
                        unsafe {
                            *depth.add(idx) = z;
                            *pixels.add(idx) = lit_color;
                        }
                    }
                }
            }
        }
    }
}

/// Rasterizes an analytical 3D sphere into a horizontal FramebufferBand.
pub fn draw_sphere_band(
    band: &mut crate::render::buffer::FramebufferBand<'_>,
    center_screen: (f32, f32, f32),
    radius_screen: f32,
    base_color: PixelColor,
    lighting: &Lighting,
) {
    if radius_screen < 0.35 {
        return;
    }
    let (cx, cy, cz) = center_screen;
    let r = radius_screen;
    let local_cy = cy - band.y_offset as f32;

    if cx + r < 0.0
        || cx - r >= band.width as f32
        || local_cy + r < 0.0
        || local_cy - r >= band.height as f32
    {
        return;
    }

    let r_sq = r * r;
    let min_x = ((cx - r).floor() as i32).max(0);
    let max_x = ((cx + r).ceil() as i32).min((band.width as i32) - 1);
    let min_y = ((local_cy - r).floor() as i32).max(0);
    let max_y = ((local_cy + r).ceil() as i32).min((band.height as i32) - 1);

    if min_x > max_x || min_y > max_y {
        return;
    }

    let inv_r = 1.0 / r;
    let min_depth = cz - r;
    let max_depth = cz + r;
    let z_near = min_depth;

    let width = band.width;
    let depth = band.depth.as_mut_ptr();
    let pixels = band.pixels.as_mut_ptr();

    for y in min_y..=max_y {
        let dy = (y as f32 + 0.5) - local_cy;
        let dy_sq = dy * dy;
        let row = (y as usize) * width;
        for x in min_x..=max_x {
            let idx = row + x as usize;
            let depth_cur = unsafe { *depth.add(idx) };
            if z_near < depth_cur {
                let dx = (x as f32 + 0.5) - cx;
                let dist_sq = dx * dx + dy_sq;
                if dist_sq <= r_sq {
                    let dz = (r_sq - dist_sq).sqrt();
                    let z = cz - dz;
                    if z < depth_cur {
                        let normal = Vec3::new(dx * inv_r, -dy * inv_r, dz * inv_r);
                        let lit_color = lighting.shade(normal, z, base_color, min_depth, max_depth);
                        unsafe {
                            *depth.add(idx) = z;
                            *pixels.add(idx) = lit_color;
                        }
                    }
                }
            }
        }
    }
}

/// Draws an unlit overlay sphere that always wins the depth test (selection markers).
pub fn draw_overlay_sphere(
    buffer: &mut Framebuffer,
    center_screen: (f32, f32, f32),
    radius_screen: f32,
    color: PixelColor,
) {
    if radius_screen <= 0.0 {
        return;
    }
    let (cx, cy, _) = center_screen;
    let r = radius_screen;
    let r_sq = r * r;
    let min_x = ((cx - r).floor() as i32).max(0);
    let max_x = ((cx + r).ceil() as i32).min((buffer.width as i32) - 1);
    let min_y = ((cy - r).floor() as i32).max(0);
    let max_y = ((cy + r).ceil() as i32).min((buffer.height as i32) - 1);
    if min_x > max_x || min_y > max_y {
        return;
    }
    for y in min_y..=max_y {
        let dy = (y as f32 + 0.5) - cy;
        let dy_sq = dy * dy;
        for x in min_x..=max_x {
            let dx = (x as f32 + 0.5) - cx;
            if dx * dx + dy_sq <= r_sq {
                buffer.set_pixel(x, y, -1.0, color);
            }
        }
    }
}

/// A screen-space point with linear depth: `(x, y, z)`.
pub type ScreenPoint = (f32, f32, f32);

/// Clips the segment `p1 -> p2` (screen coordinates with linear depth) to the
/// rectangle `(0, 0)-(w, h)` using Liang-Barsky parametric clipping.
///
/// Returns the clipped endpoints (z interpolated linearly along the segment),
/// or `None` if the segment lies entirely outside the viewport.
///
/// This bounds Bresenham-style stepping for lines whose endpoints project far
/// outside the frame (e.g. a bond partner just inside the near plane projects
/// to coordinates of 10^7+ pixels): without clipping, a single
/// [`draw_line_3d`] call would iterate hundreds of millions of times while
/// writing nothing, hanging the frame. After clipping, the step count is
/// inherently bounded by the viewport diagonal.
pub fn clip_segment_to_screen(
    p1: ScreenPoint,
    p2: ScreenPoint,
    w: usize,
    h: usize,
) -> Option<(ScreenPoint, ScreenPoint)> {
    let (x1, y1, z1) = p1;
    let (x2, y2, z2) = p2;
    let dx = x2 - x1;
    let dy = y2 - y1;

    // Parametric clip against the four half-planes x >= 0, x <= w, y >= 0, y <= h.
    let mut t0 = 0.0_f32;
    let mut t1 = 1.0_f32;
    let edges = [
        (-dx, x1),
        (dx, w as f32 - x1),
        (-dy, y1),
        (dy, h as f32 - y1),
    ];
    for &(p, q) in &edges {
        if p == 0.0 {
            if q < 0.0 {
                return None;
            }
        } else {
            let r = q / p;
            if p < 0.0 {
                t0 = t0.max(r);
            } else {
                t1 = t1.min(r);
            }
            if t0 > t1 {
                return None;
            }
        }
    }

    let lerp = |a: f32, b: f32, t: f32| a + (b - a) * t;
    Some((
        (lerp(x1, x2, t0), lerp(y1, y2, t0), lerp(z1, z2, t0)),
        (lerp(x1, x2, t1), lerp(y1, y2, t1), lerp(z1, z2, t1)),
    ))
}

/// Upper bound on useful Bresenham steps for a segment inside a `w * h`
/// viewport: its visible length cannot exceed the diagonal.
fn max_line_steps(w: usize, h: usize) -> usize {
    ((w as f32).hypot(h as f32)).ceil().max(1.0) as usize
}

/// Draws an overlay line that always wins the depth test.
pub fn draw_overlay_line(
    buffer: &mut Framebuffer,
    p1_screen: (f32, f32, f32),
    p2_screen: (f32, f32, f32),
    color: PixelColor,
) {
    let Some((a, b)) = clip_segment_to_screen(p1_screen, p2_screen, buffer.width, buffer.height)
    else {
        return;
    };
    let (x1, y1, _) = a;
    let (x2, y2, _) = b;
    let dx = x2 - x1;
    let dy = y2 - y1;
    let steps = (dx.abs().max(dy.abs()).ceil() as usize)
        .min(max_line_steps(buffer.width, buffer.height))
        .max(1);
    let steps_f = steps as f32;
    let x_step = dx / steps_f;
    let y_step = dy / steps_f;
    for i in 0..=steps {
        let fi = i as f32;
        let x = (x1 + fi * x_step).round() as i32;
        let y = (y1 + fi * y_step).round() as i32;
        buffer.set_pixel(x, y, -1.0, color);
    }
}

/// Rasterizes a 3D line between two screen points with linear depth interpolation.
///
/// `p1_screen`, `p2_screen`: `(screen_x, screen_y, view_depth_z)`
pub fn draw_line_3d(
    buffer: &mut Framebuffer,
    p1_screen: (f32, f32, f32),
    p2_screen: (f32, f32, f32),
    color: PixelColor,
) {
    let Some((a, b)) = clip_segment_to_screen(p1_screen, p2_screen, buffer.width, buffer.height)
    else {
        return;
    };
    let (x1, y1, z1) = a;
    let (x2, y2, z2) = b;

    let dx = x2 - x1;
    let dy = y2 - y1;
    let dz = z2 - z1;

    let steps = (dx.abs().max(dy.abs()).ceil() as usize)
        .min(max_line_steps(buffer.width, buffer.height))
        .max(1);
    let steps_f = steps as f32;

    let x_step = dx / steps_f;
    let y_step = dy / steps_f;
    let z_step = dz / steps_f;

    for i in 0..=steps {
        let fi = i as f32;
        let x = (x1 + fi * x_step).round() as i32;
        let y = (y1 + fi * y_step).round() as i32;
        let z = z1 + fi * z_step;

        buffer.set_pixel(x, y, z, color);
    }
}

/// Rasterizes a 3D dashed line between two screen points with linear depth interpolation.
pub fn draw_dashed_line_3d(
    buffer: &mut Framebuffer,
    p1_screen: (f32, f32, f32),
    p2_screen: (f32, f32, f32),
    color: PixelColor,
    dash_px: f32,
    gap_px: f32,
) {
    let Some((a, b)) = clip_segment_to_screen(p1_screen, p2_screen, buffer.width, buffer.height)
    else {
        return;
    };
    let (x1, y1, z1) = a;
    let (x2, y2, z2) = b;

    let dx = x2 - x1;
    let dy = y2 - y1;
    let dz = z2 - z1;

    let len = (dx * dx + dy * dy).sqrt();
    let steps = (dx.abs().max(dy.abs()).ceil() as usize)
        .min(max_line_steps(buffer.width, buffer.height))
        .max(1);
    let steps_f = steps as f32;

    let x_step = dx / steps_f;
    let y_step = dy / steps_f;
    let z_step = dz / steps_f;
    let cycle = (dash_px + gap_px).max(1.0);

    for i in 0..=steps {
        let fi = i as f32;
        let dist = if steps_f > 0.0 {
            (fi / steps_f) * len
        } else {
            0.0
        };
        let phase = dist.rem_euclid(cycle);

        if phase < dash_px {
            let x = (x1 + fi * x_step).round() as i32;
            let y = (y1 + fi * y_step).round() as i32;
            let z = z1 + fi * z_step;
            buffer.set_pixel(x, y, z, color);
        }
    }
}

/// Rasterizes a thick 3D cylinder/capsule connecting two screen endpoints with analytical normal and depth shading.
///
/// `p1_screen`, `p2_screen`: `(screen_x, screen_y, view_depth_z)`
/// `radius_screen`: cylinder radius in screen pixels
pub fn draw_cylinder(
    buffer: &mut Framebuffer,
    p1_screen: (f32, f32, f32),
    p2_screen: (f32, f32, f32),
    radius_screen: f32,
    color: PixelColor,
    lighting: &Lighting,
) {
    if radius_screen <= 0.5 {
        draw_line_3d(buffer, p1_screen, p2_screen, color);
        return;
    }

    let (x1, y1, z1) = p1_screen;
    let (x2, y2, z2) = p2_screen;
    let r = radius_screen;
    let r_sq = r * r;

    let min_x = (x1.min(x2) - r).floor() as i32;
    let max_x = (x1.max(x2) + r).ceil() as i32;
    let min_y = (y1.min(y2) - r).floor() as i32;
    let max_y = (y1.max(y2) + r).ceil() as i32;

    let min_x = min_x.max(0);
    let max_x = max_x.min((buffer.width as i32) - 1);
    let min_y = min_y.max(0);
    let max_y = max_y.min((buffer.height as i32) - 1);

    if min_x > max_x || min_y > max_y {
        return;
    }

    let ab_x = x2 - x1;
    let ab_y = y2 - y1;
    let len_sq = ab_x * ab_x + ab_y * ab_y;
    let dz_axis = z2 - z1;
    let inv_r = 1.0 / r;

    let min_depth = z1.min(z2) - r;
    let max_depth = z1.max(z2) + r;
    // Nearest possible z of this capsule. Cheap early reject for occluded pixels.
    let z_near = min_depth;

    let width = buffer.width;
    let depth = buffer.depth.as_mut_ptr();
    let pixels = buffer.pixels.as_mut_ptr();

    for y in min_y..=max_y {
        let py = y as f32 + 0.5;
        let ap_y = py - y1;
        let row = (y as usize) * width;
        for x in min_x..=max_x {
            // SAFETY: x in [0, width-1] and y in [0, height-1] (clamped above).
            let idx = row + x as usize;
            let depth_cur = unsafe { *depth.add(idx) };
            if z_near < depth_cur {
                let px = x as f32 + 0.5;
                let ap_x = px - x1;

                let t = if len_sq > 1e-6 {
                    ((ap_x * ab_x + ap_y * ab_y) / len_sq).clamp(0.0, 1.0)
                } else {
                    0.0
                };

                let qx = x1 + t * ab_x;
                let qy = y1 + t * ab_y;

                let dx = px - qx;
                let dy = py - qy;
                let dist_sq = dx * dx + dy * dy;

                if dist_sq <= r_sq {
                    let dz = (r_sq - dist_sq).sqrt();
                    let z = (z1 + t * dz_axis) - dz;
                    if z < depth_cur {
                        let normal = Vec3::new(dx * inv_r, -dy * inv_r, dz * inv_r);
                        let lit_color = lighting.shade(normal, z, color, min_depth, max_depth);
                        unsafe {
                            *depth.add(idx) = z;
                            *pixels.add(idx) = lit_color;
                        }
                    }
                }
            }
        }
    }
}

/// Band-local twin of [`draw_line_3d`] for one horizontal [`FramebufferBand`].
///
/// The screen clip and the step budget are computed against the FULL frame
/// (`band.width` x `full_height`, where `full_height` is the parent
/// framebuffer's pixel height), so every sampled point -- and therefore every
/// written pixel -- matches the serial path bit-for-bit. Only the final store
/// is translated into band-local coordinates. The sampled step range is then
/// narrowed to the band's row window (rounding-safe by construction): the
/// skipped samples all fall outside the band and would have been rejected by
/// `set_pixel` anyway, so thin lines cost O(their slice of the frame) per band
/// instead of O(the whole diagonal).
pub fn draw_line_3d_band(
    band: &mut crate::render::buffer::FramebufferBand<'_>,
    p1_screen: ScreenPoint,
    p2_screen: ScreenPoint,
    color: PixelColor,
    full_height: usize,
) {
    let Some((a, b)) = clip_segment_to_screen(p1_screen, p2_screen, band.width, full_height) else {
        return;
    };
    let (x1, y1, z1) = a;
    let (x2, y2, z2) = b;

    let dx = x2 - x1;
    let dy = y2 - y1;
    let dz = z2 - z1;

    let steps = (dx.abs().max(dy.abs()).ceil() as usize)
        .min(max_line_steps(band.width, full_height))
        .max(1);
    let steps_f = steps as f32;

    let x_step = dx / steps_f;
    let y_step = dy / steps_f;
    let z_step = dz / steps_f;

    // Narrow the step loop to samples that can land on a row this band owns.
    // Sampled points move monotonically in y, so sample i touches global row
    // round(y1 + i * y_step); inverting the band's row interval with a >= 1 px
    // safety margin (rounding + floor/ceil widening) yields a conservative
    // index window. Pixels written are identical to iterating every step.
    let y_top = band.y_offset as f32;
    let y_bot = (band.y_offset + band.height) as f32;
    let (first, last) = if y_step == 0.0 {
        if y1 < y_top - 1.0 || y1 > y_bot + 1.0 {
            return;
        }
        (0, steps)
    } else {
        let s = (y_top - 2.0 - y1) / y_step;
        let e = (y_bot + 2.0 - y1) / y_step;
        let (lo, hi) = if y_step > 0.0 { (s, e) } else { (e, s) };
        let first = (lo.floor().max(0.0) as usize).min(steps);
        let last = (hi.ceil().max(0.0) as usize).min(steps);
        (first, last)
    };

    for i in first..=last {
        let fi = i as f32;
        let x = (x1 + fi * x_step).round() as i32;
        let y = (y1 + fi * y_step).round() as i32;
        let z = z1 + fi * z_step;

        band.set_pixel(x, y - band.y_offset as i32, z, color);
    }
}

/// Rasterizes a thick 3D cylinder/capsule into one horizontal [`FramebufferBand`].
///
/// Band-aware mirror of [`draw_cylinder`] for the parallel bond passes. Every
/// floating-point expression (axis projection, analytical depth, shading) is
/// computed in FULL-frame screen coordinates exactly like the serial call --
/// shifting endpoints into band space first would change the rounding of the
/// axis projection and break bit-exact parity -- and only the final store is
/// translated onto the band's slices. The row range is simply clamped to the
/// band instead of the whole frame. `full_height` is the parent framebuffer's
/// pixel height, used only by the thin-line fallback so its screen clip
/// matches [`draw_line_3d`] exactly.
pub fn draw_cylinder_band(
    band: &mut crate::render::buffer::FramebufferBand<'_>,
    p1_screen: ScreenPoint,
    p2_screen: ScreenPoint,
    radius_screen: f32,
    color: PixelColor,
    lighting: &Lighting,
    full_height: usize,
) {
    if radius_screen <= 0.5 {
        draw_line_3d_band(band, p1_screen, p2_screen, color, full_height);
        return;
    }

    let (x1, y1, z1) = p1_screen;
    let (x2, y2, z2) = p2_screen;
    let r = radius_screen;
    let r_sq = r * r;

    let min_x = (x1.min(x2) - r).floor() as i32;
    let max_x = (x1.max(x2) + r).ceil() as i32;
    let min_y = (y1.min(y2) - r).floor() as i32;
    let max_y = (y1.max(y2) + r).ceil() as i32;

    let min_x = min_x.max(0);
    let max_x = max_x.min((band.width as i32) - 1);
    // Clamp rows to this band's global y range instead of the whole frame.
    let min_y = min_y.max(band.y_offset as i32);
    let max_y = max_y.min((band.y_offset + band.height) as i32 - 1);

    if min_x > max_x || min_y > max_y {
        return;
    }

    let ab_x = x2 - x1;
    let ab_y = y2 - y1;
    let len_sq = ab_x * ab_x + ab_y * ab_y;
    let dz_axis = z2 - z1;
    let inv_r = 1.0 / r;

    let min_depth = z1.min(z2) - r;
    let max_depth = z1.max(z2) + r;
    // Nearest possible z of this capsule. Cheap early reject for occluded pixels.
    let z_near = min_depth;

    let width = band.width;
    let y_base = band.y_offset;
    let depth = band.depth.as_mut_ptr();
    let pixels = band.pixels.as_mut_ptr();

    for y in min_y..=max_y {
        let py = y as f32 + 0.5;
        let ap_y = py - y1;
        let row = (y as usize - y_base) * width;
        for x in min_x..=max_x {
            // SAFETY: x in [0, width-1] and y-y_base in [0, band.height-1]
            // (clamped above), so the band-local index is in bounds.
            let idx = row + x as usize;
            let depth_cur = unsafe { *depth.add(idx) };
            if z_near < depth_cur {
                let px = x as f32 + 0.5;
                let ap_x = px - x1;

                let t = if len_sq > 1e-6 {
                    ((ap_x * ab_x + ap_y * ab_y) / len_sq).clamp(0.0, 1.0)
                } else {
                    0.0
                };

                let qx = x1 + t * ab_x;
                let qy = y1 + t * ab_y;

                let dx = px - qx;
                let dy = py - qy;
                let dist_sq = dx * dx + dy * dy;

                if dist_sq <= r_sq {
                    let dz = (r_sq - dist_sq).sqrt();
                    let z = (z1 + t * dz_axis) - dz;
                    if z < depth_cur {
                        let normal = Vec3::new(dx * inv_r, -dy * inv_r, dz * inv_r);
                        let lit_color = lighting.shade(normal, z, color, min_depth, max_depth);
                        unsafe {
                            *depth.add(idx) = z;
                            *pixels.add(idx) = lit_color;
                        }
                    }
                }
            }
        }
    }
}

/// Rasterizes a 3D triangle with barycentric depth interpolation and Lambertian lighting.
///
/// `v1`, `v2`, `v3`: `(screen_x, screen_y, view_depth_z)`
pub fn draw_triangle_3d(
    buffer: &mut Framebuffer,
    v1: (f32, f32, f32),
    v2: (f32, f32, f32),
    v3: (f32, f32, f32),
    normal: Vec3,
    color: PixelColor,
    lighting: &Lighting,
) {
    let min_x = (v1.0.min(v2.0).min(v3.0)).floor() as i32;
    let max_x = (v1.0.max(v2.0).max(v3.0)).ceil() as i32;
    let min_y = (v1.1.min(v2.1).min(v3.1)).floor() as i32;
    let max_y = (v1.1.max(v2.1).max(v3.1)).ceil() as i32;

    let min_x = min_x.max(0);
    let max_x = max_x.min((buffer.width as i32) - 1);
    let min_y = min_y.max(0);
    let max_y = max_y.min((buffer.height as i32) - 1);

    if min_x > max_x || min_y > max_y {
        return;
    }

    let area = (v2.0 - v1.0) * (v3.1 - v1.1) - (v2.1 - v1.1) * (v3.0 - v1.0);
    if area.abs() < 1e-6 {
        return;
    }
    let inv_area = 1.0 / area;

    let min_depth = v1.2.min(v2.2).min(v3.2);
    let max_depth = v1.2.max(v2.2).max(v3.2);
    // Nearest possible z of this triangle. Cheap early reject for occluded pixels.
    let z_near = min_depth;

    let width = buffer.width;
    let depth = buffer.depth.as_mut_ptr();
    let pixels = buffer.pixels.as_mut_ptr();

    for y in min_y..=max_y {
        let py = y as f32 + 0.5;
        let row = (y as usize) * width;
        for x in min_x..=max_x {
            // SAFETY: x in [0, width-1] and y in [0, height-1] (clamped above).
            let idx = row + x as usize;
            let depth_cur = unsafe { *depth.add(idx) };
            if z_near < depth_cur {
                let px = x as f32 + 0.5;

                let w1 = ((v2.0 - px) * (v3.1 - py) - (v2.1 - py) * (v3.0 - px)) * inv_area;
                let w2 = ((v3.0 - px) * (v1.1 - py) - (v3.1 - py) * (v1.0 - px)) * inv_area;
                let w3 = 1.0 - w1 - w2;

                if w1 >= -1e-4 && w2 >= -1e-4 && w3 >= -1e-4 {
                    let z = w1 * v1.2 + w2 * v2.2 + w3 * v3.2;
                    if z < depth_cur {
                        let lit_color = lighting.shade(normal, z, color, min_depth, max_depth);
                        unsafe {
                            *depth.add(idx) = z;
                            *pixels.add(idx) = lit_color;
                        }
                    }
                }
            }
        }
    }
}

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
    if radius_screen <= 0.0 {
        return;
    }

    let (cx, cy, cz) = center_screen;
    let r = radius_screen;
    let r_sq = r * r;

    let min_x = ((cx - r).floor() as i32).max(0);
    let max_x = ((cx + r).ceil() as i32).min((buffer.width as i32) - 1);
    let min_y = ((cy - r).floor() as i32).max(0);
    let max_y = ((cy + r).ceil() as i32).min((buffer.height as i32) - 1);

    if min_x > max_x || min_y > max_y {
        return;
    }

    let min_depth = cz - r;
    let max_depth = cz + r;

    for y in min_y..=max_y {
        let dy = (y as f32 + 0.5) - cy;
        let dy_sq = dy * dy;

        for x in min_x..=max_x {
            let dx = (x as f32 + 0.5) - cx;
            let dist_sq = dx * dx + dy_sq;

            if dist_sq <= r_sq {
                let dz = (r_sq - dist_sq).sqrt();
                let z = cz - dz;

                // Normal in camera view space (+X right, +Y up, +Z towards camera)
                let normal = Vec3::new(dx / r, -dy / r, dz / r);
                let lit_color = lighting.shade(normal, z, base_color, min_depth, max_depth);

                buffer.set_pixel(x, y, z, lit_color);
            }
        }
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
    let (x1, y1, z1) = p1_screen;
    let (x2, y2, z2) = p2_screen;

    let dx = x2 - x1;
    let dy = y2 - y1;
    let dz = z2 - z1;

    let steps = (dx.abs().max(dy.abs()).ceil() as usize).max(1);
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

    let min_depth = z1.min(z2) - r;
    let max_depth = z1.max(z2) + r;

    for y in min_y..=max_y {
        let py = y as f32 + 0.5;
        let ap_y = py - y1;

        for x in min_x..=max_x {
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
                let z_axis = z1 + t * (z2 - z1);
                let z = z_axis - dz;

                let normal = Vec3::new(dx / r, -dy / r, dz / r);
                let lit_color = lighting.shade(normal, z, color, min_depth, max_depth);

                buffer.set_pixel(x, y, z, lit_color);
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

    for y in min_y..=max_y {
        let py = y as f32 + 0.5;

        for x in min_x..=max_x {
            let px = x as f32 + 0.5;

            let w1 = ((v2.0 - px) * (v3.1 - py) - (v2.1 - py) * (v3.0 - px)) * inv_area;
            let w2 = ((v3.0 - px) * (v1.1 - py) - (v3.1 - py) * (v1.0 - px)) * inv_area;
            let w3 = 1.0 - w1 - w2;

            if w1 >= -1e-4 && w2 >= -1e-4 && w3 >= -1e-4 {
                let z = w1 * v1.2 + w2 * v2.2 + w3 * v3.2;
                let lit_color = lighting.shade(normal, z, color, min_depth, max_depth);
                buffer.set_pixel(x, y, z, lit_color);
            }
        }
    }
}

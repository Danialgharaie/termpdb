//! Screen-space post-processing: silhouette depth outlining and SSAO.

use crate::render::buffer::Framebuffer;

/// Configuration for post-processing filters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PostProcessConfig {
    /// Enable dark silhouette outline along depth edges and borders.
    pub outline: bool,
    /// Enable screen-space ambient occlusion (SSAO) crevice darkening.
    pub ssao: bool,
    /// Relative depth jump threshold for detecting silhouette edges (default: 0.12).
    pub outline_threshold: f32,
    /// Sampling radius in pixels for SSAO (default: 2).
    pub ssao_radius: usize,
}

impl Default for PostProcessConfig {
    fn default() -> Self {
        Self {
            outline: true,
            ssao: true,
            outline_threshold: 0.12,
            ssao_radius: 2,
        }
    }
}

/// Applies configured post-processing passes (outline and SSAO) to the framebuffer in place.
pub fn apply_postprocessing(fb: &mut Framebuffer, config: &PostProcessConfig) {
    if !config.outline && !config.ssao {
        return;
    }

    let w = fb.width;
    let h = fb.height;
    if w < 2 || h < 2 {
        return;
    }

    let orig_pixels = fb.pixels.clone();
    let depth = &fb.depth;
    let outline = config.outline;
    let ssao = config.ssao;
    let outline_threshold = config.outline_threshold;
    let ssao_radius = config.ssao_radius;

    use rayon::prelude::*;

    fb.pixels
        .par_chunks_mut(w)
        .enumerate()
        .for_each(|(y, row)| {
            for (x, pixel) in row.iter_mut().enumerate() {
                let idx = y * w + x;
                let z = depth[idx];
                if !z.is_finite() {
                    continue;
                }

                let mut pixel_scale = 1.0f32;

                // 1. Silhouette / Depth Jump Outline
                if outline {
                    let mut is_edge = false;
                    // Check 4-connectivity
                    let neighbors = [
                        (x.wrapping_sub(1), y, x > 0),
                        (x + 1, y, x + 1 < w),
                        (x, y.wrapping_sub(1), y > 0),
                        (x, y + 1, y + 1 < h),
                    ];

                    for (nx, ny, valid) in neighbors {
                        if !valid {
                            is_edge = true;
                            break;
                        }
                        let n_idx = ny * w + nx;
                        let nz = depth[n_idx];
                        if !nz.is_finite() {
                            is_edge = true;
                            break;
                        }
                        // Depth jump test
                        let diff = (z - nz).abs();
                        let min_z = z.min(nz).max(0.01);
                        if diff / min_z > outline_threshold {
                            is_edge = true;
                            break;
                        }
                    }

                    if is_edge {
                        pixel_scale *= 0.35;
                    }
                }

                // 2. Screen-Space Ambient Occlusion (SSAO)
                if ssao {
                    let r = ssao_radius as isize;
                    let mut total_samples = 0;
                    let mut occluded = 0;

                    for dy in -r..=r {
                        for dx in -r..=r {
                            if dx == 0 && dy == 0 {
                                continue;
                            }
                            let sx = x as isize + dx;
                            let sy = y as isize + dy;
                            if sx >= 0 && sx < w as isize && sy >= 0 && sy < h as isize {
                                total_samples += 1;
                                let s_idx = (sy as usize) * w + (sx as usize);
                                let sz = depth[s_idx];
                                if sz.is_finite() && sz < z && (z - sz) < 5.0 {
                                    occluded += 1;
                                }
                            }
                        }
                    }

                    if total_samples > 0 && occluded > 0 {
                        let occ_ratio = occluded as f32 / total_samples as f32;
                        let ao_factor = (1.0 - occ_ratio * 0.55).clamp(0.45, 1.0);
                        pixel_scale *= ao_factor;
                    }
                }

                if pixel_scale < 0.999 {
                    let p = orig_pixels[idx];
                    let r = ((p.0 as f32) * pixel_scale).round().clamp(0.0, 255.0) as u8;
                    let g = ((p.1 as f32) * pixel_scale).round().clamp(0.0, 255.0) as u8;
                    let b = ((p.2 as f32) * pixel_scale).round().clamp(0.0, 255.0) as u8;
                    *pixel = (r, g, b);
                }
            }
        });
}

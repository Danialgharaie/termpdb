use crate::math::vec3::Vec3;

#[derive(Debug, Clone)]
pub struct CatmullRomSpline {
    pub points: Vec<Vec3>,
}

impl CatmullRomSpline {
    pub fn new(points: Vec<Vec3>) -> Self {
        Self { points }
    }

    /// Interpolates between `p1` and `p2` given neighbor points `p0` and `p3` for `t` in [0, 1].
    pub fn interpolate_segment(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32) -> Vec3 {
        let t2 = t * t;
        let t3 = t2 * t;

        0.5 * ((2.0 * p1)
            + (-p0 + p2) * t
            + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
            + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3)
    }

    /// Calculates tangent vector at parameter `t` in [0, 1] between `p1` and `p2`.
    pub fn tangent_segment(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32) -> Vec3 {
        let t2 = t * t;

        let raw_tangent = 0.5
            * ((-p0 + p2)
                + 2.0 * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t
                + 3.0 * (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t2);

        let n = raw_tangent.norm();
        if n > 1e-8 {
            raw_tangent / n
        } else {
            let direct = (p2 - p1).normalize();
            if direct.norm() > 1e-8 {
                direct
            } else {
                Vec3::new(0.0, 0.0, 1.0)
            }
        }
    }

    /// Evaluates the spline at global parameter `t` in [0, 1].
    pub fn interpolate(&self, t: f32) -> Vec3 {
        let n = self.points.len();
        if n == 0 {
            return Vec3::zero();
        }
        if n == 1 {
            return self.points[0];
        }
        if n == 2 {
            return self.points[0].lerp(self.points[1], t.clamp(0.0, 1.0));
        }

        let num_segments = n - 1;
        let t_clamped = t.clamp(0.0, 1.0);

        let (seg, u) = if t_clamped >= 1.0 {
            (num_segments - 1, 1.0)
        } else {
            let val = t_clamped * (num_segments as f32);
            let idx = val.floor() as usize;
            let u = val - (idx as f32);
            (idx.min(num_segments - 1), u)
        };

        let p1 = self.points[seg];
        let p2 = self.points[seg + 1];
        let p0 = if seg > 0 {
            self.points[seg - 1]
        } else {
            p1 * 2.0 - p2
        };
        let p3 = if seg + 2 < n {
            self.points[seg + 2]
        } else {
            p2 * 2.0 - p1
        };

        Self::interpolate_segment(p0, p1, p2, p3, u)
    }

    /// Evaluates the tangent at global parameter `t` in [0, 1].
    pub fn tangent(&self, t: f32) -> Vec3 {
        let n = self.points.len();
        if n == 0 {
            return Vec3::new(0.0, 0.0, 1.0);
        }
        if n == 1 {
            return Vec3::new(0.0, 0.0, 1.0);
        }
        if n == 2 {
            let dir = (self.points[1] - self.points[0]).normalize();
            return if dir.norm() > 1e-8 {
                dir
            } else {
                Vec3::new(0.0, 0.0, 1.0)
            };
        }

        let num_segments = n - 1;
        let t_clamped = t.clamp(0.0, 1.0);

        let (seg, u) = if t_clamped >= 1.0 {
            (num_segments - 1, 1.0)
        } else {
            let val = t_clamped * (num_segments as f32);
            let idx = val.floor() as usize;
            let u = val - (idx as f32);
            (idx.min(num_segments - 1), u)
        };

        let p1 = self.points[seg];
        let p2 = self.points[seg + 1];
        let p0 = if seg > 0 {
            self.points[seg - 1]
        } else {
            p1 * 2.0 - p2
        };
        let p3 = if seg + 2 < n {
            self.points[seg + 2]
        } else {
            p2 * 2.0 - p1
        };

        Self::tangent_segment(p0, p1, p2, p3, u)
    }

    /// Generates a smooth sampled curve with `samples_per_segment` samples between each pair of control points.
    pub fn generate_smooth_curve(&self, samples_per_segment: usize) -> Vec<Vec3> {
        let n = self.points.len();
        if n == 0 {
            return Vec::new();
        }
        if n == 1 {
            return vec![self.points[0]];
        }

        let samples = samples_per_segment.max(1);
        let num_segments = n - 1;
        let mut curve = Vec::with_capacity(num_segments * samples + 1);

        for seg in 0..num_segments {
            let p1 = self.points[seg];
            let p2 = self.points[seg + 1];
            let p0 = if seg > 0 {
                self.points[seg - 1]
            } else {
                p1 * 2.0 - p2
            };
            let p3 = if seg + 2 < n {
                self.points[seg + 2]
            } else {
                p2 * 2.0 - p1
            };

            for k in 0..samples {
                let u = (k as f32) / (samples as f32);
                curve.push(Self::interpolate_segment(p0, p1, p2, p3, u));
            }
        }

        curve.push(*self.points.last().unwrap());
        curve
    }

    /// Generates points and corresponding normalized tangents along the curve.
    pub fn generate_smooth_curve_with_tangents(
        &self,
        samples_per_segment: usize,
    ) -> Vec<(Vec3, Vec3)> {
        let n = self.points.len();
        if n == 0 {
            return Vec::new();
        }
        if n == 1 {
            return vec![(self.points[0], Vec3::new(0.0, 0.0, 1.0))];
        }

        let samples = samples_per_segment.max(1);
        let num_segments = n - 1;
        let mut curve = Vec::with_capacity(num_segments * samples + 1);

        for seg in 0..num_segments {
            let p1 = self.points[seg];
            let p2 = self.points[seg + 1];
            let p0 = if seg > 0 {
                self.points[seg - 1]
            } else {
                p1 * 2.0 - p2
            };
            let p3 = if seg + 2 < n {
                self.points[seg + 2]
            } else {
                p2 * 2.0 - p1
            };

            for k in 0..samples {
                let u = (k as f32) / (samples as f32);
                let pt = Self::interpolate_segment(p0, p1, p2, p3, u);
                let tan = Self::tangent_segment(p0, p1, p2, p3, u);
                curve.push((pt, tan));
            }
        }

        let last_pt = *self.points.last().unwrap();
        let last_tan = self.tangent(1.0);
        curve.push((last_pt, last_tan));

        curve
    }
}

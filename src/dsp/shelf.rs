//! The single filter primitive: a first-order pole/zero shelf in the Brown & Duda form
//! `H(s) = (w0 + az*s) / (w0 + ap*s)`, bilinear-transformed with a prewarped corner.
//!
//! Properties relied on elsewhere (and asserted in `tests/`):
//!   * DC gain is exactly 1 for every `(az, ap)` and every sample rate,
//!   * the high-frequency plateau is `az/ap`,
//!   * the pole is strictly inside the unit circle for all positive parameters.

#[derive(Clone, Copy, Debug, Default)]
pub struct Shelf1 {
    b0: f32,
    b1: f32,
    a1: f32,
    z: f32,
}

impl Shelf1 {
    /// `w0` is the analog corner in rad/s; `az`/`ap` are the zero/pole multipliers.
    pub fn set(&mut self, az: f64, ap: f64, w0: f64, fs: f64) {
        let k = 2.0 * fs;
        let w0 = k * (w0 / k).tan();
        let d = w0 + ap * k;
        self.b0 = ((w0 + az * k) / d) as f32;
        self.b1 = ((w0 - az * k) / d) as f32;
        self.a1 = ((w0 - ap * k) / d) as f32;
    }

    pub fn reset(&mut self) {
        self.z = 0.0;
    }

    /// Pole magnitude, for stability assertions.
    pub fn pole_radius(&self) -> f32 {
        self.a1.abs()
    }

    /// Coefficient access for the GUI's frequency-response drawing.
    pub fn coeffs(&self) -> (f32, f32, f32) {
        (self.b0, self.b1, self.a1)
    }

    /// Magnitude at `f`, for the editor only (never called from `process`).
    pub fn magnitude_at(&self, f: f64, fs: f64) -> f64 {
        let w = 2.0 * std::f64::consts::PI * f / fs;
        let (sn, cs) = (-w).sin_cos();
        let num = (self.b0 as f64) + (self.b1 as f64) * cs;
        let num_i = (self.b1 as f64) * sn;
        let den = 1.0 + (self.a1 as f64) * cs;
        let den_i = (self.a1 as f64) * sn;
        ((num * num + num_i * num_i) / (den * den + den_i * den_i)).sqrt()
    }

    #[inline]
    pub fn step(&mut self, x: f32) -> f32 {
        let y = self.b0 * x + self.z;
        self.z = self.b1 * x - self.a1 * y;
        y
    }
}

//! Fractional delay for the crossfeed path: an integer ring read plus a 3rd-order Thiran
//! all-pass carrying the fractional remainder.
//!
//! Why 3rd order: the all-pass must be *continuous* as the delay moves (the speaker-angle
//! knob sweeps D from ~6 to ~23 samples), and every Thiran order degenerates into a
//! marginal pole at integer delays. With `n = floor(D) - 2` and `d in [2, 3)` the
//! fractional part never sits near a degenerate point except exactly at an integer, so
//! clamping into `[2.001, 2.999]` costs at most 0.001 samples (< 0.03 us) and never
//! introduces a delay discontinuity.
//!
//! Group delay is maximally flat at DC, which is where the ITD cue lives.

pub struct FracDelay {
    buf: Box<[f32]>,
    mask: usize,
    w: usize,
    n: usize,
    // Thiran coefficients b0..b2 and a1..a3 (b3 = 1)
    b0: f32,
    b1: f32,
    b2: f32,
    a1: f32,
    a2: f32,
    a3: f32,
    s1: f32,
    s2: f32,
    s3: f32,
}

impl FracDelay {
    /// `max_delay_samples` is the largest total delay that will ever be requested.
    pub fn with_capacity(max_delay_samples: f64) -> Self {
        let size = (max_delay_samples.ceil() as usize + 4).next_power_of_two();
        Self {
            buf: vec![0.0; size].into_boxed_slice(),
            mask: size - 1,
            w: 0,
            n: 2,
            b0: 0.0,
            b1: 0.0,
            b2: 0.0,
            a1: 0.0,
            a2: 0.0,
            a3: 0.0,
            s1: 0.0,
            s2: 0.0,
            s3: 0.0,
        }
    }

    pub fn reset(&mut self) {
        self.buf.fill(0.0);
        self.w = 0;
        self.s1 = 0.0;
        self.s2 = 0.0;
        self.s3 = 0.0;
    }

    pub fn max_delay_samples(&self) -> f64 {
        (self.buf.len() - 4) as f64
    }

    /// Total group delay in samples (at DC).
    pub fn set_delay(&mut self, d_samples: f64) {
        let n = ((d_samples.floor() - 2.0).max(0.0)) as usize;
        let d = (d_samples - n as f64).clamp(2.001, 2.999);

        // Thiran 3rd order, maximally flat group delay at DC:
        //   a_k = (-1)^k C(3,k) prod_{m=0..3} (d-3+m) / (d-3+k+m)
        let a1 = -3.0 * (d - 3.0) / (d + 1.0);
        let a2 = 3.0 * (d - 3.0) * (d - 2.0) / ((d + 1.0) * (d + 2.0));
        let a3 = -(d - 3.0) * (d - 2.0) * (d - 1.0) / ((d + 1.0) * (d + 2.0) * (d + 3.0));

        self.n = n.min(self.buf.len() - 4);
        self.a1 = a1 as f32;
        self.a2 = a2 as f32;
        self.a3 = a3 as f32;
        self.b0 = a3 as f32;
        self.b1 = a2 as f32;
        self.b2 = a1 as f32;
    }

    #[inline]
    pub fn step(&mut self, x: f32) -> f32 {
        self.buf[self.w] = x;
        let xd = self.buf[(self.w + self.buf.len() - self.n) & self.mask];
        self.w = (self.w + 1) & self.mask;

        let y = self.b0 * xd + self.s1;
        self.s1 = self.b1 * xd - self.a1 * y + self.s2;
        self.s2 = self.b2 * xd - self.a2 * y + self.s3;
        self.s3 = xd - self.a3 * y;
        y
    }
}

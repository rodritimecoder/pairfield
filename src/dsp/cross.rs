//! The cross path `C = ITF(z) * z^-D`, isolated so it can be measured directly.

use super::shelf::Shelf1;
use super::thiran::FracDelay;
use super::{alpha_pair, woodworth_delay_s, HEAD_W0};

pub struct CrossPath {
    filt: Shelf1,
    delay: FracDelay,
    alpha_near: f64,
    alpha_far: f64,
    delay_samples: f64,
}

impl CrossPath {
    pub fn new(max_delay_samples: f64) -> Self {
        Self {
            filt: Shelf1::default(),
            delay: FracDelay::with_capacity(max_delay_samples),
            alpha_near: 1.0,
            alpha_far: 1.0,
            delay_samples: 0.0,
        }
    }

    pub fn reset(&mut self) {
        self.filt.reset();
        self.delay.reset();
    }

    pub fn set(&mut self, speaker_deg: f64, fs: f64) {
        let (an, af) = alpha_pair(speaker_deg);
        self.alpha_near = an;
        self.alpha_far = af;
        self.filt.set(af, an, HEAD_W0, fs);
        self.delay_samples = woodworth_delay_s(speaker_deg) * fs;
        self.delay.set_delay(self.delay_samples);
    }

    pub fn alpha_pair(&self) -> (f64, f64) {
        (self.alpha_near, self.alpha_far)
    }

    pub fn delay_samples(&self) -> f64 {
        self.delay_samples
    }

    #[inline]
    pub fn step(&mut self, x: f32) -> f32 {
        self.delay.step(self.filt.step(x))
    }
}

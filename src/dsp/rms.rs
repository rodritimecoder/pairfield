//! Program-dependent level matching so Bypass is an honest A/B.
//!
//! The analytic DC normaliser in the kernel already pins centred material to unity, so this
//! stage only has to absorb the residue: the comb's energy loss, panned content that does not
//! get the full crossfeed, and the reflection bus. It is a feed-forward measurement (the wet
//! signal is measured *before* the gain is applied), so there is no servo loop to ring.

pub struct RmsMatcher {
    det: f32,
    smooth: f32,
    p_dry: f32,
    p_wet: f32,
    gain: f32,
}

const FLOOR: f32 = 1e-10;

impl RmsMatcher {
    pub fn new(fs: f64) -> Self {
        let mut m = Self {
            det: 0.0,
            smooth: 0.0,
            p_dry: 0.0,
            p_wet: 0.0,
            gain: 1.0,
        };
        m.set_rate(fs);
        m
    }

    pub fn set_rate(&mut self, fs: f64) {
        // 300 ms power averaging, 30 ms gain smoothing.
        self.det = (1.0 - (-1.0 / (0.30 * fs)).exp()) as f32;
        self.smooth = (1.0 - (-1.0 / (0.03 * fs)).exp()) as f32;
    }

    pub fn reset(&mut self) {
        self.p_dry = 0.0;
        self.p_wet = 0.0;
        self.gain = 1.0;
    }

    pub fn gain(&self) -> f32 {
        self.gain
    }

    /// `dry_power`/`wet_power` are the mean of the squared samples across channels.
    #[inline]
    pub fn step(&mut self, dry_power: f32, wet_power: f32, enabled: bool) -> f32 {
        self.p_dry += (dry_power - self.p_dry) * self.det;
        self.p_wet += (wet_power - self.p_wet) * self.det;

        let target = if !enabled || self.p_wet < FLOOR || self.p_dry < FLOOR {
            // Hold, rather than slew to 1, so silence cannot ratchet the gain.
            self.gain
        } else {
            let g = (self.p_dry / self.p_wet).sqrt().clamp(0.3548, 2.8184); // +-9 dB
            let g = if enabled { g } else { 1.0 };
            self.gain + (g - self.gain) * self.smooth
        };

        if !enabled {
            // Auto-gain off: converge on unity so the knob click is click-free.
            self.gain += (1.0 - self.gain) * self.smooth;
            return self.gain;
        }
        self.gain = target;
        self.gain
    }
}

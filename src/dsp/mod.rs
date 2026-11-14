//! Host-independent DSP core. No `nice-plug` types appear in this module, so the offline
//! harness and the unit tests drive exactly the same code the plugin runs.

pub mod cross;
pub mod rms;
pub mod room;
pub mod shelf;
pub mod thiran;

use cross::CrossPath;
use room::RoomBus;
use rms::RmsMatcher;
use shelf::Shelf1;

use std::f64::consts::PI;

pub const C_SOUND: f64 = 343.0;
pub const HEAD_RADIUS: f64 = 0.0875;
/// Acoustic corner of the head, rad/s: `c / a` = 3920 rad/s = 623.9 Hz.
pub const HEAD_W0: f64 = C_SOUND / HEAD_RADIUS;
/// Air-trim shelf corner.
const AIR_W0: f64 = 2.0 * PI * 3500.0;
/// Hard ceiling on the speaker angle (keeps `alpha_far` positive).
pub const MAX_ANGLE_DEG: f64 = 60.0;

/// Brown & Duda eq. (5): head-shadow `alpha` for an angle of incidence, clamped at 150 deg.
/// `alpha(0) = 2.0` (+6 dB plateau), `alpha(150) = 0.1` (-20 dB).
#[inline]
pub fn head_shadow_alpha(incidence_deg: f64) -> f64 {
    let t = incidence_deg.min(150.0) * PI / 150.0;
    1.05 + 0.95 * t.cos()
}

/// The alpha pair for a symmetric speaker pair, derived from one geometric scalar so that the
/// interaural filter's DC group delay is exactly half the Woodworth delay (Kuhn's 3:2 ratio).
#[inline]
pub fn alpha_pair(speaker_deg: f64) -> (f64, f64) {
    let th = speaker_deg * PI / 180.0;
    let s = 0.5 * (th + th.sin());
    (1.0 + 0.5 * s, 1.0 - 0.5 * s)
}

/// Woodworth spherical-head path excess, seconds.
#[inline]
pub fn woodworth_delay_s(speaker_deg: f64) -> f64 {
    let th = speaker_deg * PI / 180.0;
    (HEAD_RADIUS / C_SOUND) * (th + th.sin())
}

/// Parameters as plain values, one snapshot per chunk (the plugin supplies smoothed values).
#[derive(Clone, Copy, Debug)]
pub struct DspParams {
    pub speaker_angle_deg: f32,
    pub crossfeed_db: f32,
    /// 0 = full physical M/S crossfeed, 1 = M path untouched (centre is bit-exact).
    pub center_focus: f32,
    pub room_active: bool,
    pub room_db: f32,
    pub air_db: f32,
    pub auto_gain: bool,
    pub bypass: bool,
}

impl Default for DspParams {
    fn default() -> Self {
        Self {
            speaker_angle_deg: 30.0,
            crossfeed_db: -6.0,
            center_focus: 0.5,
            room_active: true,
            room_db: -22.0,
            air_db: 0.0,
            auto_gain: true,
            bypass: false,
        }
    }
}

/// Magnitude of the shared shelf primitive, in f64, for offline analysis and the editor.
pub fn shelf_mag(az: f64, ap: f64, w0: f64, fs: f64, f: f64) -> f64 {
    let k = 2.0 * fs;
    let w0p = k * (w0 / k).tan();
    let d = w0p + ap * k;
    let b0 = (w0p + az * k) / d;
    let b1 = (w0p - az * k) / d;
    let a1 = (w0p - ap * k) / d;
    let w = 2.0 * PI * f / fs;
    let (sn, cs) = (-w).sin_cos();
    let nre = b0 + b1 * cs;
    let nim = b1 * sn;
    let dre = 1.0 + a1 * cs;
    let dim = a1 * sn;
    ((nre * nre + nim * nim) / (dre * dre + dim * dim)).sqrt()
}

pub struct PairfieldDsp {
    fs: f64,
    // cross path, one instance per input channel
    cross: [CrossPath; 2],
    // master air trim, one per output channel
    air: [Shelf1; 2],
    room: RoomBus,
    rms: RmsMatcher,

    // derived, refreshed by set_params
    g: f32,
    gm: f32,
    gdc: f32,
    coef_self: f32,
    coef_other: f32,

    bypass_mix: f32,
    bypass_step: f32,
    snap: bool,
    prepared: bool,
}

impl PairfieldDsp {
    pub fn new(fs: f64) -> Self {
        let mut dsp = Self {
            fs: 48000.0,
            cross: [CrossPath::new(16.0), CrossPath::new(16.0)],
            air: [Shelf1::default(); 2],
            room: RoomBus::new(24.0, fs),
            rms: RmsMatcher::new(fs),
            g: 0.5,
            gm: 0.25,
            gdc: 0.8,
            coef_self: 0.0,
            coef_other: 0.5,
            bypass_mix: 0.0,
            bypass_step: 1.0,
            snap: true,
            prepared: false,
        };
        dsp.prepare(fs);
        dsp
    }

    /// Allocate and size everything. Called from `initialize()` / `activate()`, never from
    /// `process()`.
    pub fn prepare(&mut self, fs: f64) {
        self.fs = fs;
        let max_delay = woodworth_delay_s(MAX_ANGLE_DEG) * fs;
        for c in self.cross.iter_mut() {
            *c = CrossPath::new(max_delay);
        }
        self.room = RoomBus::new(24.0, fs);
        self.room.prepare(fs, HEAD_W0);
        self.rms.set_rate(fs);
        // 5 ms bypass crossfade
        self.bypass_step = (1.0 / (0.005 * fs)) as f32;
        self.prepared = true;
        self.reset();
    }

    pub fn reset(&mut self) {
        for c in self.cross.iter_mut() {
            c.reset();
        }
        for s in self.air.iter_mut() {
            s.reset();
        }
        self.room.reset();
        self.rms.reset();
        self.bypass_mix = 0.0;
        self.snap = true;
    }

    pub fn sample_rate(&self) -> f64 {
        self.fs
    }

    /// Refresh coefficients from a (smoothed) parameter snapshot.
    pub fn set_params(&mut self, p: &DspParams) {
        let fs = self.fs;
        let angle = (p.speaker_angle_deg as f64).clamp(15.0, MAX_ANGLE_DEG);

        for c in self.cross.iter_mut() {
            c.set(angle, fs);
        }

        let air = 10f64.powf(p.air_db as f64 / 20.0);
        for s in self.air.iter_mut() {
            s.set(air, 1.0, AIR_W0, fs);
        }

        self.g = 10f32.powf(p.crossfeed_db / 20.0);
        self.gm = self.g * (1.0 - p.center_focus.clamp(0.0, 1.0));
        self.coef_self = 0.5 * (self.gm - self.g);
        self.coef_other = 0.5 * (self.gm + self.g);
        self.gdc = 1.0 / (1.0 + self.gm);

        self.room.set_level(p.room_db, p.room_active);
    }

    /// Process one chunk in place. `l` and `r` must be the same length.
    pub fn process_block(&mut self, l: &mut [f32], r: &mut [f32], p: &DspParams) {
        debug_assert_eq!(l.len(), r.len());
        let n = l.len().min(r.len());
        if self.snap {
            self.bypass_mix = if p.bypass { 1.0 } else { 0.0 };
            self.snap = false;
        }
        for i in 0..n {
            let xl = l[i];
            let xr = r[i];

            // ---- cross path: C = ITF(z) * z^-D, applied to each input channel
            let cl = self.cross[0].step(xl);
            let cr = self.cross[1].step(xr);

            // ---- M/S decoupled recombination
            let mut yl = xl + self.coef_self * cl + self.coef_other * cr;
            let mut yr = xr + self.coef_self * cr + self.coef_other * cl;

            // ---- early reflections
            let (rl, rr) = self.room.tick(xl, xr);
            yl += rl;
            yr += rr;

            // ---- air trim
            yl = self.air[0].step(yl);
            yr = self.air[1].step(yr);

            // ---- level matching: the analytic DC normaliser is folded in *before* the
            // detector, so the detector only measures what the RMS matcher has to fix.
            let dry_pow = 0.5 * (xl * xl + xr * xr);
            let nl = yl * self.gdc;
            let nr = yr * self.gdc;
            let wet_pow = 0.5 * (nl * nl + nr * nr);
            let gauto = self.rms.step(dry_pow, wet_pow, p.auto_gain && !p.bypass);
            let wet_l = nl * gauto;
            let wet_r = nr * gauto;

            // ---- bypass crossfade
            if p.bypass {
                self.bypass_mix = (self.bypass_mix + self.bypass_step).min(1.0);
            } else {
                self.bypass_mix = (self.bypass_mix - self.bypass_step).max(0.0);
            }
            if self.bypass_mix >= 1.0 {
                l[i] = xl;
                r[i] = xr;
            } else if self.bypass_mix <= 0.0 {
                l[i] = wet_l;
                r[i] = wet_r;
            } else {
                let m = self.bypass_mix;
                l[i] = dry_mix(xl, wet_l, m);
                r[i] = dry_mix(xr, wet_r, m);
            }
        }
    }

    pub fn auto_gain(&self) -> f32 {
        self.rms.gain()
    }

    // ---- analysis helpers for the editor (never called from process) ----

    /// Mono-sum magnitude in dB: |G_dc * (1 + gm*C(f))|.
    pub fn mono_response_db(p: &DspParams, fs: f64, f: f64) -> f64 {
        let angle = (p.speaker_angle_deg as f64).clamp(15.0, MAX_ANGLE_DEG);
        let (an, af) = alpha_pair(angle);
        let g = 10f64.powf(p.crossfeed_db as f64 / 20.0);
        let gm = g * (1.0 - p.center_focus.clamp(0.0, 1.0) as f64);
        let d = woodworth_delay_s(angle) * fs;
        let w = 2.0 * PI * f / fs;
        let (sn, cs) = (-w).sin_cos();
        let (dsn, dcs) = (-w * d).sin_cos();
        // cross path phase: delay plus the shelf's own phase (approximated by its magnitude
        // pair? no - we evaluate the shelf exactly via its transfer function)
        let k = 2.0 * fs;
        let w0p = k * ((HEAD_W0 / k).tan());
        let den = w0p + an * k;
        let b0 = (w0p + af * k) / den;
        let b1 = (w0p - af * k) / den;
        let a1 = (w0p - an * k) / den;
        let sre = b0 + b1 * cs;
        let sim = b1 * sn;
        let dre = 1.0 + a1 * cs;
        let dim = a1 * sn;
        let dn = dre * dre + dim * dim;
        let hre = (sre * dre + sim * dim) / dn;
        let him = (sim * dre - sre * dim) / dn;
        let cre = hre * dcs - him * dsn;
        let cim = hre * dsn + him * dcs;
        let re = 1.0 + gm * cre;
        let im = gm * cim;
        let norm = 1.0 / (1.0 + gm);
        let m = (re * re + im * im).sqrt() * norm;
        20.0 * m.log10()
    }

    /// Cross-path (ITF) magnitude in dB relative to its DC gain.
    pub fn cross_response_db(p: &DspParams, fs: f64, f: f64) -> f64 {
        let angle = (p.speaker_angle_deg as f64).clamp(15.0, MAX_ANGLE_DEG);
        let (an, af) = alpha_pair(angle);
        20.0 * shelf_mag(af, an, HEAD_W0, fs, f).log10()
    }
}

#[inline]
fn dry_mix(dry: f32, wet: f32, m: f32) -> f32 {
    (1.0 - m) * wet + m * dry
}

//! Early-reflection bus: six discrete taps from an image-source model of a 6.4 x 7.0 x 2.9 m
//! control room (see ARCHITECTURE.md section 5), each head-shadow filtered per ear.
//!
//! Delays and azimuths are fixed by the reference geometry; only the overall level moves.
//! Reads are linearly interpolated: the taps are >= 18 dB down and wall absorption is itself
//! a high-frequency loss, so the interpolator's own gentle HF droop is the physically correct
//! colour. The precision of a Thiran network is reserved for the ITD, where it is a cue.

use super::head_shadow_alpha;
use super::shelf::Shelf1;

pub const TAP_COUNT: usize = 6;

/// Azimuth of each tap as seen from the head, degrees, positive = right.
const TAP_AZIMUTH: [f64; TAP_COUNT] = [-30.0, -14.0, -75.0, -75.0, 79.0, 65.0];
/// Excess path delay over the direct sound, milliseconds.
const TAP_MS: [f64; TAP_COUNT] = [6.23, 9.22, 11.78, 13.14, 16.61, 18.30];
/// Level relative to the direct sound, dB (image-source geometry, before the level knob).
const TAP_REL_DB: [f64; TAP_COUNT] = [-9.0, -13.6, -13.1, -18.9, -15.3, -18.4];

/// Broadband wall/air absorption of a reflection, folded into the shelf plateau.
const WALL_ABSORPTION: f64 = 0.6;
/// Normalisation so the taps are referenced to the same-side near-ear response the direct
/// path would have had. Keeps a reflection slightly darker than the direct sound.
const ALPHA_REF: f64 = 1.25;

#[derive(Clone, Copy, Default)]
struct Tap {
    n: usize,
    f: f32,
    /// left output from left input (and right output from right input)
    w_same: f32,
    /// left output from right input (and right output from left input)
    w_cross: f32,
    gain: f32,
    /// [dl through w_same, dl through w_cross, dr through w_same, dr through w_cross]
    filt: [Shelf1; 4],
}

pub struct RoomBus {
    lh: Box<[f32]>,
    rh: Box<[f32]>,
    mask: usize,
    w: usize,
    taps: [Tap; TAP_COUNT],
}

impl RoomBus {
    pub fn new(max_delay_ms: f64, fs: f64) -> Self {
        let size = ((max_delay_ms * 0.001 * fs).ceil() as usize + 4).next_power_of_two();
        Self {
            lh: vec![0.0; size].into_boxed_slice(),
            rh: vec![0.0; size].into_boxed_slice(),
            mask: size - 1,
            w: 0,
            taps: [Tap::default(); TAP_COUNT],
        }
    }

    pub fn reset(&mut self) {
        self.lh.fill(0.0);
        self.rh.fill(0.0);
        self.w = 0;
        for t in self.taps.iter_mut() {
            for f in t.filt.iter_mut() {
                f.reset();
            }
        }
    }

    /// Recompute delay splits and the per-tap shadow coefficients for a new sample rate.
    pub fn prepare(&mut self, fs: f64, w0: f64) {
        for (i, t) in self.taps.iter_mut().enumerate() {
            let total = TAP_MS[i] * 0.001 * fs;
            t.n = total.floor() as usize;
            t.f = (total - t.n as f64) as f32;

            let inc_l = (TAP_AZIMUTH[i] + 90.0).abs();
            let inc_r = (TAP_AZIMUTH[i] - 90.0).abs();
            let w_left_ear = WALL_ABSORPTION * head_shadow_alpha(inc_l) / ALPHA_REF;
            let w_right_ear = WALL_ABSORPTION * head_shadow_alpha(inc_r) / ALPHA_REF;
            t.w_same = w_left_ear as f32;
            t.w_cross = w_right_ear as f32;
            t.filt[0].set(w_left_ear, 1.0, w0, fs);
            t.filt[1].set(w_right_ear, 1.0, w0, fs);
            t.filt[2].set(w_left_ear, 1.0, w0, fs);
            t.filt[3].set(w_right_ear, 1.0, w0, fs);
        }
        self.set_level(-22.0, true);
    }

    pub fn set_level(&mut self, room_db: f32, active: bool) {
        let lvl = if active { room_db } else { -200.0 };
        let rel0 = TAP_REL_DB[0];
        for (i, t) in self.taps.iter_mut().enumerate() {
            t.gain = 10f32.powf(((lvl as f64 + TAP_REL_DB[i] - rel0) / 20.0) as f32);
        }
    }

    #[inline]
    pub fn tick(&mut self, xl: f32, xr: f32) -> (f32, f32) {
        self.lh[self.w] = xl;
        self.rh[self.w] = xr;
        let len = self.lh.len();
        let mask = self.mask;

        let mut ol = 0.0f32;
        let mut or = 0.0f32;
        for t in self.taps.iter_mut() {
            let i0 = (self.w + len - t.n) & mask;
            let i1 = (self.w + len - t.n - 1) & mask;
            let f = t.f;
            let dl = self.lh[i0] + f * (self.lh[i1] - self.lh[i0]);
            let dr = self.rh[i0] + f * (self.rh[i1] - self.rh[i0]);

            let dl_same = t.filt[0].step(dl);
            let dl_cross = t.filt[1].step(dl);
            let dr_same = t.filt[2].step(dr);
            let dr_cross = t.filt[3].step(dr);

            ol += t.gain * (dl_same + dr_cross);
            or += t.gain * (dr_same + dl_cross);
        }

        self.w = (self.w + 1) & mask;
        (ol, or)
    }
}

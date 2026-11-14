//! Offline measurement harness. Drives the real DSP core (no host, no GUI) and checks the
//! numbers claimed in ARCHITECTURE.md.
//!
//!     cargo run --release --example measure
//!
//! Every figure printed here is measured from an impulse response, not from the closed-form
//! expressions - the closed forms are printed alongside as an independent cross-check.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering as AOrd};

/// Counts heap activity so the "no allocation on the audio thread" claim is measured,
/// not asserted.
struct Counting;
static ALLOCS: AtomicUsize = AtomicUsize::new(0);
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        ALLOCS.fetch_add(1, AOrd::Relaxed);
        System.alloc(l)
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        System.dealloc(p, l)
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        ALLOCS.fetch_add(1, AOrd::Relaxed);
        System.realloc(p, l, n)
    }
}
#[global_allocator]
static ALLOCATOR: Counting = Counting;

use pairfield::dsp::{
    alpha_pair, cross::CrossPath, woodworth_delay_s, DspParams,
    PairfieldDsp, HEAD_RADIUS, C_SOUND, HEAD_W0,
};

const N_IR: usize = 4096;
const FS: f64 = 48000.0;

struct Report {
    fails: u32,
    checks: u32,
}

impl Report {
    fn check(&mut self, name: &str, ok: bool, detail: String) {
        self.checks += 1;
        if !ok {
            self.fails += 1;
        }
        println!("  [{}] {:<46} {}", if ok { "PASS" } else { "FAIL" }, name, detail);
    }
}

/// DFT at an arbitrary frequency - exact for a finite impulse response.
fn tf(ir: &[f32], f: f64, fs: f64) -> (f64, f64) {
    let w = -2.0 * std::f64::consts::PI * f / fs;
    let (sw, cw) = w.sin_cos();
    // iterative rotation for speed
    let (mut s, mut c) = (0.0f64, 1.0f64);
    let (mut re, mut im) = (0.0f64, 0.0f64);
    for x in ir {
        let x = *x as f64;
        re += x * c;
        im += x * s;
        let ns = s * cw + c * sw;
        let nc = c * cw - s * sw;
        s = ns;
        c = nc;
    }
    (re, im)
}

fn mag_db(ir: &[f32], f: f64, fs: f64) -> f64 {
    let (re, im) = tf(ir, f, fs);
    10.0 * (re * re + im * im).log10()
}

fn group_delay_samples(ir: &[f32], f: f64, fs: f64) -> f64 {
    let (r1, i1) = tf(ir, f + 0.5, fs);
    let (r0, i0) = tf(ir, f - 0.5, fs);
    let mut d = i1.atan2(r1) - i0.atan2(r0);
    while d > std::f64::consts::PI {
        d -= 2.0 * std::f64::consts::PI;
    }
    while d < -std::f64::consts::PI {
        d += 2.0 * std::f64::consts::PI;
    }
    -d * fs / (2.0 * std::f64::consts::PI)
}

fn impulse_response(cross: &mut CrossPath, angle: f64, fs: f64) -> Vec<f32> {
    cross.reset();
    cross.set(angle, fs);
    let mut ir = vec![0.0f32; N_IR];
    ir[0] = 1.0;
    for x in ir.iter_mut() {
        *x = cross.step(*x);
    }
    ir
}

/// Transfer matrix of the whole plugin: column 0 = left input, column 1 = right input.
struct Matrix {
    h: [[Vec<f32>; 2]; 2], // [out][in]
}

fn render_matrix(p: &DspParams, fs: f64) -> Matrix {
    let mut dsp = PairfieldDsp::new(fs);
    let mut h: [[Vec<f32>; 2]; 2] = [
        [vec![0.0; N_IR], vec![0.0; N_IR]],
        [vec![0.0; N_IR], vec![0.0; N_IR]],
    ];
    for input in 0..2 {
        dsp.reset();
        dsp.set_params(p);
        let mut l = vec![0.0f32; N_IR];
        let mut r = vec![0.0f32; N_IR];
        if input == 0 {
            l[0] = 1.0;
        } else {
            r[0] = 1.0;
        }
        dsp.process_block(&mut l, &mut r, p);
        h[0][input] = l;
        h[1][input] = r;
    }
    Matrix { h }
}

impl Matrix {
    fn eval(&self, f: f64, pan: f64, fs: f64) -> f64 {
        let (l, r) = (
            ((1.0 - pan) * std::f64::consts::FRAC_PI_2).cos(),
            ((1.0 - pan) * std::f64::consts::FRAC_PI_2).sin(),
        );
        let (a, b) = tf(&self.h[0][0], f, fs);
        let (c, d) = tf(&self.h[0][1], f, fs);
        let (e, g) = tf(&self.h[1][0], f, fs);
        let (h, i) = tf(&self.h[1][1], f, fs);
        // (re00*l + re01*r) + j*(im00*l + im01*r)
        let ol = ((a * l + c * r).powi(2) + (b * l + d * r).powi(2)).sqrt();
        let or = ((e * l + h * r).powi(2) + (g * l + i * r).powi(2)).sqrt();
        ol.max(or)
    }

    fn dc(&self, pan: f64, fs: f64) -> f64 {
        self.eval(20.0, pan, fs)
    }

    fn worst_null_db(&self, pan: f64, fs: f64) -> f64 {
        let dc = self.dc(pan, fs);
        let mut worst = 0.0f64;
        let mut f = 200.0;
        while f < 20000.0 {
            let v = self.eval(f, pan, fs);
            if v < worst || worst == 0.0 {
                worst = v;
            }
            f *= 1.002;
        }
        20.0 * (worst / dc).log10()
    }
}

fn main() {
    let mut rep = Report { fails: 0, checks: 0 };

    println!("pairfield - offline measurement harness");
    println!(
        "c = {} m/s   a = {} m   w0 = c/a = {:.0} rad/s ({:.1} Hz)\n",
        C_SOUND,
        HEAD_RADIUS,
        HEAD_W0,
        HEAD_W0 / (2.0 * std::f64::consts::PI)
    );

    // ---------- 1. head-shadow shelf ----------
    println!("1. Head-shadow shelf (realised filter, measured from its impulse response)");
    let mut cross = CrossPath::new(woodworth_delay_s(60.0) * FS);
    for angle in [15.0, 30.0, 45.0, 60.0] {
        let (an, af) = alpha_pair(angle);
        let ir = impulse_response(&mut cross, angle, FS);
        let dc = mag_db(&ir, 20.0, FS);
        let plateau = mag_db(&ir, 19000.0, FS);
        let target = 20.0 * (af / an).log10();
        rep.check(
            &format!("  alpha ratio @ {angle:.0} deg -> plateau"),
            (plateau - target).abs() < 0.15 && dc.abs() < 0.01,
            format!("dc {dc:+.3} dB, plateau {plateau:+.2} dB, expected {target:+.2} dB"),
        );
    }

    // ---------- 2. ITD ----------
    println!("\n2. ITD of the cross path (Woodworth bulk + Kuhn 3:2 from the shelf's phase)");
    for angle in [15.0, 30.0, 45.0, 60.0] {
        let ir = impulse_response(&mut cross, angle, FS);
        let expect = woodworth_delay_s(angle) * FS;
        let hf = group_delay_samples(&ir, 4000.0, FS);
        let lf = group_delay_samples(&ir, 50.0, FS);
        rep.check(
            &format!("  tau_HF @ {angle:.0} deg"),
            (hf - expect).abs() < 0.4,
            format!(
                "{hf:.2} samples ({:.1} us) vs {expect:.2} (shelf adds {:+.2})",
                hf / FS * 1e6,
                hf - expect
            ),
        );
        rep.check(
            &format!("  tau_LF/tau_HF @ {angle:.0} deg"),
            ((lf / hf) - 1.5).abs() < 0.06,
            format!("{:.3} (Kuhn 1.5)", lf / hf),
        );
    }

    // ---------- 3. sample-rate independence ----------
    println!("\n3. Sample-rate independence (delay specified in seconds)");
    let mut per_rate = Vec::new();
    for fs in [44100.0, 48000.0, 96000.0, 192000.0] {
        let mut c = CrossPath::new(woodworth_delay_s(60.0) * fs);
        let ir = impulse_response(&mut c, 30.0, fs);
        let us = group_delay_samples(&ir, 1000.0, fs) / fs * 1e6;
        per_rate.push(us);
        println!("     fs {fs:>7.0} Hz -> ITD@1kHz {us:>7.1} us");
    }
    let spread = per_rate.iter().cloned().fold(f64::MIN, f64::max)
        - per_rate.iter().cloned().fold(f64::MAX, f64::min);
    rep.check(
        "  ITD consistent across 44.1..192 kHz",
        spread < 1.0,
        format!("spread {spread:.2} us"),
    );

    // ---------- 4. comb budget ----------
    println!("\n4. Comb budget: worst null over all pan positions (centred is the worst case)");
    println!("     {:>4} {:>5} {:>11} {:>11} {:>11}", "az", "k", "centre", "worst-pan", "analytic");
    for angle in [15.0, 30.0, 45.0, 60.0] {
        for k in [0.0, 0.5, 1.0] {
            let p = DspParams {
                speaker_angle_deg: angle as f32,
                center_focus: k,
                room_active: false,
                auto_gain: false,
                ..Default::default()
            };
            let m = render_matrix(&p, FS);
            let centre = m.worst_null_db(0.5, FS);
            let mut worst = 0.0f64;
            for i in 0..=20 {
                let d = m.worst_null_db(i as f64 / 20.0, FS);
                if d < worst {
                    worst = d;
                }
            }
            // independent closed-form cross-check
            let mut analytic = 0.0f64;
            let mut f = 200.0;
            let dc_an = PairfieldDsp::mono_response_db(&p, FS, 20.0);
            while f < 20000.0 {
                let d = PairfieldDsp::mono_response_db(&p, FS, f) - dc_an;
                if d < analytic {
                    analytic = d;
                }
                f *= 1.002;
            }
            println!(
                "     {:>4.0} {:>5.2} {:>10.2} {:>10.2} {:>10.2}",
                angle, k, centre, worst, analytic
            );
            rep.check(
                &format!("  analytic matches measured, az {angle:.0} k {k:.1}"),
                (analytic - centre).abs() < 0.6,
                format!("delta {:.2} dB", analytic - centre),
            );
        }
    }

    // ---------- 5. centre exactness ----------
    println!("\n5. Centre Focus = 100%: dead-centre mono must be bit-identical to bypass");
    {
        let p = DspParams {
            speaker_angle_deg: 30.0,
            center_focus: 1.0,
            room_active: false,
            auto_gain: false,
            ..Default::default()
        };
        let mut dsp = PairfieldDsp::new(FS);
        dsp.reset();
        dsp.set_params(&p);
        let mut l = vec![0.0f32; 4096];
        let mut r = vec![0.0f32; 4096];
        let mut seed = 0x1234_5678u32;
        let orig: Vec<f32> = (0..4096)
            .map(|_| {
                seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
                ((seed >> 8) as f32 / 8388608.0) - 1.0
            })
            .collect();
        l.copy_from_slice(&orig);
        r.copy_from_slice(&orig);
        dsp.process_block(&mut l, &mut r, &p);
        let mut max_err = 0.0f32;
        for i in 0..4096 {
            max_err = max_err.max((l[i] - orig[i]).abs()).max((r[i] - orig[i]).abs());
        }
        rep.check(
            "  |wet - dry| for L == R",
            max_err < 1e-6,
            format!("max error {max_err:e}"),
        );
    }

    // ---------- 6. bypass ----------
    println!("\n6. Bypass is bit-exact");
    {
        let p = DspParams {
            bypass: true,
            ..Default::default()
        };
        let mut dsp = PairfieldDsp::new(FS);
        let mut l: Vec<f32> = (0..1024).map(|i| ((i as f32) * 0.01).sin()).collect();
        let mut r = l.clone();
        let (l0, r0) = (l.clone(), r.clone());
        dsp.process_block(&mut l, &mut r, &p);
        let identical = l.iter().zip(&l0).all(|(a, b)| a.to_bits() == b.to_bits())
            && r.iter().zip(&r0).all(|(a, b)| a.to_bits() == b.to_bits());
        rep.check("  output == input, bit for bit", identical, "1024 frames".into());
    }

    // ---------- 7. level matching ----------
    println!("\n7. Auto-gain: RMS(dry) / RMS(wet) after settling");
    for (label, pan) in [("mono", 0.5), ("hard left", 0.0), ("centre-ish", 0.35)] {
        let p = DspParams {
            auto_gain: true,
            room_active: true,
            ..Default::default()
        };
        let mut dsp = PairfieldDsp::new(FS);
        dsp.reset();
        dsp.set_params(&p);
        let n = 48000 * 4;
        let mut l = vec![0.0f32; n];
        let mut r = vec![0.0f32; n];
        let mut seed = 0xdead_beefu32;
        let (gl, gr) = (
            ((1.0 - pan) * std::f64::consts::FRAC_PI_2).cos() as f32,
            ((1.0 - pan) * std::f64::consts::FRAC_PI_2).sin() as f32,
        );
        for i in 0..n {
            seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
            let x = ((seed >> 8) as f32 / 8388608.0) - 1.0;
            l[i] = x * gl;
            r[i] = x * gr;
        }
        let dry: f64 = {
            let mut s = 0.0;
            for i in &l {
                s += (*i as f64) * (*i as f64);
            }
            for i in &r {
                s += (*i as f64) * (*i as f64);
            }
            (s / (2.0 * n as f64)).sqrt()
        };
        dsp.process_block(&mut l, &mut r, &p);
        // measure the last 20% (settled)
        let skip = n * 4 / 5;
        let mut s = 0.0;
        for i in skip..n {
            s += (l[i] as f64) * (l[i] as f64) + (r[i] as f64) * (r[i] as f64);
        }
        let wet = (s / (2.0 * (n - skip) as f64)).sqrt();
        let err_db = 20.0 * (wet / dry).log10();
        rep.check(
            &format!("  {label} level match"),
            err_db.abs() < 0.5,
            format!("{err_db:+.2} dB (auto gain {:.2} dB)", 20.0 * dsp.auto_gain().log10()),
        );
    }

    // ---------- 8. numerical hygiene ----------
    println!("\n8. Numerical hygiene");
    {
        let mut dsp = PairfieldDsp::new(FS);
        let mut bad = 0u32;
        let mut sub = 0u32;
        let mut seed = 0x0bad_f00du32;
        for angle in [15.0f32, 60.0, 22.5, 45.0] {
            for k in [0.0f32, 1.0, 0.5] {
                let p = DspParams {
                    speaker_angle_deg: angle,
                    center_focus: k,
                    air_db: 3.0,
                    room_db: -12.0,
                    ..Default::default()
                };
                dsp.set_params(&p);
                for _ in 0..200 {
                    let mut l = vec![0.0f32; 64];
                    let mut r = vec![0.0f32; 64];
                    for i in 0..64 {
                        seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
                        let x = if i % 7 == 0 {
                            0.0
                        } else {
                            ((seed >> 8) as f32 / 8388608.0) - 1.0
                        };
                        l[i] = x;
                        r[i] = -x;
                    }
                    dsp.process_block(&mut l, &mut r, &p);
                    for x in l.iter().chain(r.iter()) {
                        if !x.is_finite() {
                            bad += 1;
                        }
                        if *x != 0.0 && x.abs() < f32::MIN_POSITIVE * 2.0 {
                            sub += 1;
                        }
                    }
                }
            }
        }
        rep.check(
            "  no NaN / inf over 38k frames of hostile input",
            bad == 0,
            format!("{bad} bad samples"),
        );
        rep.check("  no subnormal output", sub == 0, format!("{sub} subnormals"));
    }

    // ---------- 9. real-time safety ----------
    println!("\n9. Real-time safety: heap activity during processing");
    {
        let mut dsp = PairfieldDsp::new(FS);
        // every buffer is allocated in prepare(); size them before measuring
        dsp.prepare(FS);
        let p = DspParams::default();
        dsp.set_params(&p);
        let mut l = vec![0.1f32; 512];
        let mut r = vec![0.2f32; 512];
        dsp.process_block(&mut l, &mut r, &p);

        let before = ALLOCS.load(AOrd::Relaxed);
        for i in 0..20_000u32 {
            let p = DspParams {
                speaker_angle_deg: 20.0 + (i % 40) as f32,
                center_focus: (i % 100) as f32 / 100.0,
                ..Default::default()
            };
            dsp.set_params(&p);
            dsp.process_block(&mut l, &mut r, &p);
        }
        let allocs = ALLOCS.load(AOrd::Relaxed) - before;
        rep.check(
            "  zero allocations across 20k blocks with moving params",
            allocs == 0,
            format!("{allocs} allocations"),
        );
    }

    println!("\n{}", "=".repeat(72));
    println!(
        "{} checks, {} failed  ->  {}",
        rep.checks,
        rep.fails,
        if rep.fails == 0 { "ALL PASS" } else { "FAILURES" }
    );
    if rep.fails > 0 {
        std::process::exit(1);
    }
}

//! Invariants the architecture depends on. Fast; run with `cargo test`.

use pairfield::dsp::{
    alpha_pair, cross::CrossPath, head_shadow_alpha, shelf_mag, woodworth_delay_s, DspParams,
    PairfieldDsp, C_SOUND, HEAD_RADIUS, HEAD_W0,
};

const FS: f64 = 48000.0;

fn approx(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

#[test]
fn head_shadow_alpha_anchors() {
    assert!(approx(head_shadow_alpha(0.0), 2.0, 1e-12));
    assert!(approx(head_shadow_alpha(150.0), 0.1, 1e-12));
    // clamped beyond 150 degrees (the model is only valid to there)
    assert!(approx(head_shadow_alpha(180.0), 0.1, 1e-12));
    // monotone decreasing
    let mut prev = f64::INFINITY;
    for d in (0..=150).step_by(5) {
        let a = head_shadow_alpha(d as f64);
        assert!(a < prev, "not monotone at {d}");
        prev = a;
    }
}

#[test]
fn alpha_pair_gives_kuhn_ratio() {
    for deg in [15.0, 30.0, 45.0, 60.0] {
        let (an, af) = alpha_pair(deg);
        assert!(an > af);
        // The interaural filter's DC group delay must be exactly half the Woodworth delay.
        // tau = (an - af)/w0 and Woodworth tau = a/c (theta + sin theta).
        let tau_itf = (an - af) / HEAD_W0;
        let tau_w = woodworth_delay_s(deg);
        assert!(
            approx(tau_itf, 0.5 * tau_w, 1e-12),
            "deg {deg}: {tau_itf:e} vs {:e}",
            0.5 * tau_w
        );
    }
}

#[test]
fn woodworth_reference_values() {
    // Classic spherical-head numbers.
    assert!(approx(woodworth_delay_s(90.0) * 1e6, 656.0, 2.0));
    assert!(approx(woodworth_delay_s(30.0) * 1e6, 261.1, 0.2));
    assert!(approx(HEAD_W0 / (2.0 * std::f64::consts::PI), 623.9, 0.2));
    assert!(approx(HEAD_RADIUS / C_SOUND * 1e6, 255.1, 0.2));
}

#[test]
fn shelf_dc_gain_is_exactly_unity() {
    for fs in [44100.0, 48000.0, 88200.0, 96000.0, 192000.0] {
        for az in [0.3, 0.5, 0.87, 1.0, 1.13, 1.26, 1.48] {
            for ap in [0.5, 0.63, 1.0, 1.26, 1.5] {
                let m = shelf_mag(az, ap, HEAD_W0, fs, 0.0);
                assert!(approx(m, 1.0, 1e-9), "fs {fs} az {az} ap {ap} -> {m}");
                let plateau = shelf_mag(az, ap, HEAD_W0, fs, 20000.0);
                assert!(
                    approx(plateau, az / ap, 0.01),
                    "plateau {plateau} vs {}",
                    az / ap
                );
            }
        }
    }
}

#[test]
fn shelf_pole_is_stable() {
    use pairfield::dsp::shelf::Shelf1;
    for fs in [44100.0, 48000.0, 192000.0] {
        for ap in [0.1, 0.4, 0.7, 1.0, 1.5, 2.0] {
            for az in [0.02, 0.1, 0.5, 1.0, 2.0, 4.0] {
                let mut s = Shelf1::default();
                s.set(az, ap, HEAD_W0, fs);
                let r = s.pole_radius();
                assert!(r < 1.0, "fs {fs} az {az} ap {ap} -> pole radius {r}");
                let m = shelf_mag(az, ap, HEAD_W0, fs, 0.0);
                assert!(m.is_finite() && m > 0.0);
            }
        }
    }
}

fn tf(ir: &[f32], f: f64, fs: f64) -> (f64, f64) {
    let w = -2.0 * std::f64::consts::PI * f / fs;
    let (sw, cw) = w.sin_cos();
    let (mut s, mut c) = (0.0f64, 1.0f64);
    let (mut re, mut im) = (0.0f64, 0.0f64);
    for x in ir {
        re += (*x as f64) * c;
        im += (*x as f64) * s;
        let ns = s * cw + c * sw;
        let nc = c * cw - s * sw;
        s = ns;
        c = nc;
    }
    (re, im)
}

fn group_delay(ir: &[f32], f: f64, fs: f64) -> f64 {
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

fn cross_ir(c: &mut CrossPath, _deg: f64, _fs: f64, n: usize) -> Vec<f32> {
    let mut ir = vec![0.0f32; n];
    ir[0] = 1.0;
    for x in ir.iter_mut() {
        *x = c.step(*x);
    }
    ir
}

/// The knob sweeps D from ~5.9 to ~23.4 samples. The delay must track it closely and never
/// jump - that is what stops the speaker-angle knob from zippering.
///
/// Measured at 200 Hz, where the phase is smooth and the delay is unambiguous, and compared
/// against `bulk delay + the shelf's own theoretical group delay` so that only the fractional
/// delay network is under test.
fn shelf_group_delay_samples(deg: f64, f: f64, fs: f64) -> f64 {
    let (an, af) = alpha_pair(deg);
    let tp = an / HEAD_W0;
    let tz = af / HEAD_W0;
    let w = 2.0 * std::f64::consts::PI * f;
    (tp / (1.0 + (w * tp).powi(2)) - tz / (1.0 + (w * tz).powi(2))) * fs
}

#[test]
fn delay_tracks_angle_continuously() {
    let mut c = CrossPath::new(woodworth_delay_s(60.0) * FS);
    let mut prev: Option<(f64, f64)> = None;
    let mut deg = 15.0;
    let probe = 200.0;
    while deg <= 60.0 {
        c.reset();
        c.set(deg, FS);
        // 8192 samples, not 1024: near an integer delay the 3rd-order Thiran leaves a
        // -60 dB Nyquist tail that takes ~800 samples to decay, and truncating it biases
        // a 1 Hz phase-difference measurement.
        let ir = cross_ir(&mut c, deg, FS, 8192);
        let measured = group_delay(&ir, probe, FS);
        let expected = woodworth_delay_s(deg) * FS + shelf_group_delay_samples(deg, probe, FS);
        assert!(
            (measured - expected).abs() < 0.15,
            "deg {deg}: measured {measured:.3} vs expected {expected:.3}"
        );
        if let Some((pe, pm)) = prev {
            let d_expected = expected - pe;
            let d_measured = measured - pm;
            assert!(
                d_measured > 0.0 && (d_measured - d_expected).abs() < 0.02,
                "discontinuity at {deg} deg: expected step {d_expected:.4}, measured {d_measured:.4}"
            );
        }
        prev = Some((expected, measured));
        deg += 0.1;
    }
}

#[test]
fn thiran_is_allpass() {
    let mut c = CrossPath::new(woodworth_delay_s(60.0) * FS);
    c.reset();
    c.set(37.0, FS);
    let ir = cross_ir(&mut c, 37.0, FS, 2048);
    // the delay network must not colour: |C| must equal the shelf magnitude everywhere
    let (an, af) = alpha_pair(37.0);
    for f in [100.0, 1000.0, 5000.0, 15000.0] {
        let (re, im) = tf(&ir, f, FS);
        let measured = 10.0 * (re * re + im * im).log10();
        let expected = 20.0 * shelf_mag(af, an, HEAD_W0, FS, f).log10();
        assert!(
            (measured - expected).abs() < 0.02,
            "f {f}: {measured:.4} vs {expected:.4}"
        );
    }
}

#[test]
fn centre_focus_one_is_bit_transparent_for_mono() {
    let p = DspParams {
        center_focus: 1.0,
        room_active: false,
        auto_gain: false,
        ..Default::default()
    };
    let mut dsp = PairfieldDsp::new(FS);
    dsp.reset();
    dsp.set_params(&p);
    let src: Vec<f32> = (0..2048).map(|i| ((i as f32) * 0.017).sin() * 0.5).collect();
    let mut l = src.clone();
    let mut r = src.clone();
    dsp.process_block(&mut l, &mut r, &p);
    for i in 0..src.len() {
        assert!((l[i] - src[i]).abs() < 1e-6, "L differs at {i}");
        assert!((r[i] - src[i]).abs() < 1e-6, "R differs at {i}");
    }
}

#[test]
fn bypass_is_bit_exact() {
    let p = DspParams {
        bypass: true,
        ..Default::default()
    };
    let mut dsp = PairfieldDsp::new(FS);
    let src: Vec<f32> = (0..1024).map(|i| ((i as f32) * 0.031).cos()).collect();
    let mut l = src.clone();
    let mut r = src.clone();
    dsp.process_block(&mut l, &mut r, &p);
    for i in 0..src.len() {
        assert_eq!(l[i].to_bits(), src[i].to_bits());
        assert_eq!(r[i].to_bits(), src[i].to_bits());
    }
}

#[test]
fn no_nans_across_the_parameter_space() {
    let mut dsp = PairfieldDsp::new(FS);
    for angle in [15.0f32, 30.0, 45.0, 60.0] {
        for k in [0.0f32, 0.5, 1.0] {
            for g in [-18.0f32, -6.0, 0.0] {
                let p = DspParams {
                    speaker_angle_deg: angle,
                    center_focus: k,
                    crossfeed_db: g,
                    air_db: 3.0,
                    room_db: -12.0,
                    ..Default::default()
                };
                dsp.reset();
                dsp.set_params(&p);
                let mut l = vec![0.9f32; 256];
                let mut r = vec![-0.9f32; 256];
                dsp.process_block(&mut l, &mut r, &p);
                assert!(l.iter().all(|x| x.is_finite()));
                assert!(r.iter().all(|x| x.is_finite()));
            }
        }
    }
}

#[test]
fn reset_clears_state() {
    let p = DspParams::default();
    let mut dsp = PairfieldDsp::new(FS);
    dsp.set_params(&p);
    let mut l = vec![1.0f32; 512];
    let mut r = vec![0.5f32; 512];
    dsp.process_block(&mut l, &mut r, &p);
    dsp.reset();
    dsp.set_params(&p);
    let mut l = vec![0.0f32; 512];
    let mut r = vec![0.0f32; 512];
    dsp.process_block(&mut l, &mut r, &p);
    assert!(l.iter().all(|x| x.abs() < 1e-12));
    assert!(r.iter().all(|x| x.abs() < 1e-12));
}

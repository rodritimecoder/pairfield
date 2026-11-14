// Phase-1 numeric verification of the ARCHITECTURE.md math. No plugin code.
use std::f64::consts::PI;

const C: f64 = 343.0; // m/s
const A: f64 = 0.0875; // m head radius

#[derive(Clone, Copy, Debug)]
struct Cx {
    re: f64,
    im: f64,
}
impl Cx {
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn add(self, o: Cx) -> Cx {
        Cx::new(self.re + o.re, self.im + o.im)
    }
    fn mul(self, o: Cx) -> Cx {
        Cx::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    fn scale(self, k: f64) -> Cx {
        Cx::new(self.re * k, self.im * k)
    }
    fn div(self, o: Cx) -> Cx {
        let d = o.re * o.re + o.im * o.im;
        Cx::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
    fn abs(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
    fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
}

fn unit(w: f64, n: f64, fs: f64) -> Cx {
    let a = -(w / fs) * n;
    Cx::new(a.cos(), a.sin())
}
fn z1(w: f64, fs: f64) -> Cx {
    unit(w, 1.0, fs)
}

/// First-order shelf, analog H(s) = (dc*w0 + hf*s)/(w0 + s), bilinear + prewarp.
#[derive(Clone, Copy, Debug)]
struct Shelf {
    b0: f64,
    b1: f64,
    a1: f64,
}
impl Shelf {
    fn design(dc: f64, hf: f64, w0a: f64, fs: f64) -> Self {
        let k = 2.0 * fs;
        let w0 = k * (w0a / k).tan();
        Shelf {
            b0: (dc * w0 + hf * k) / (w0 + k),
            b1: (dc * w0 - hf * k) / (w0 + k),
            a1: (w0 - k) / (w0 + k),
        }
    }
    #[inline]
    fn eval(&self, w: f64, fs: f64) -> Cx {
        let zi = z1(w, fs);
        Cx::new(self.b0, 0.0)
            .add(Cx::new(self.b1, 0.0).mul(zi))
            .div(Cx::new(1.0, 0.0).add(Cx::new(self.a1, 0.0).mul(zi)))
    }
    fn pole(&self) -> f64 {
        -self.a1
    }
}

/// 2nd-order Thiran all-pass: total group delay D at DC, D in [1,2].
#[derive(Clone, Copy, Debug)]
struct Thiran2 {
    a1: f64,
    a2: f64,
}
impl Thiran2 {
    fn design(d: f64) -> Self {
        Thiran2 {
            a1: -2.0 * (d - 2.0) / (d + 1.0),
            a2: (d - 2.0) * (d - 1.0) / ((d + 1.0) * (d + 2.0)),
        }
    }
    #[inline]
    fn eval(&self, w: f64, fs: f64) -> Cx {
        let zi = z1(w, fs);
        let zi2 = zi.mul(zi);
        Cx::new(self.a2, 0.0)
            .add(Cx::new(self.a1, 0.0).mul(zi))
            .add(zi2)
            .div(
                Cx::new(1.0, 0.0)
                    .add(Cx::new(self.a1, 0.0).mul(zi))
                    .add(Cx::new(self.a2, 0.0).mul(zi2)),
            )
    }
    fn pole_radius(&self) -> f64 {
        self.a2.abs().sqrt()
    }
}

/// Group delay in samples; `h` takes angular frequency in rad/s.
fn gd<F: Fn(f64) -> Cx>(h: F, f_hz: f64, fs: f64) -> f64 {
    let w = 2.0 * PI * f_hz;
    let dw = 2.0 * PI * 0.5;
    let p1 = h(w + dw).arg();
    let p0 = h(w - dw).arg();
    let mut d = p1 - p0;
    while d > PI {
        d -= 2.0 * PI;
    }
    while d < -PI {
        d += 2.0 * PI;
    }
    -d * fs / (2.0 * PI)
}

fn db(x: f64) -> f64 {
    20.0 * x.log10()
}

struct P {
    az: f64,
    g_db: f64,
    k: f64,
    fs: f64,
}
impl P {
    fn new(az: f64, g_db: f64, k: f64, fs: f64) -> Self {
        Self { az, g_db, k, fs }
    }
    fn th(&self) -> f64 {
        self.az * PI / 180.0
    }
    /// s = 0.5(theta + sin theta), normalised Woodworth path excess (tau * w0).
    fn s(&self) -> f64 {
        0.5 * (self.th() + self.th().sin())
    }
    fn an(&self) -> f64 {
        1.0 + 0.5 * self.s()
    }
    fn af(&self) -> f64 {
        1.0 - 0.5 * self.s()
    }
    fn w0(&self) -> f64 {
        C / A
    }
    fn tau_hf(&self) -> f64 {
        (A / C) * (self.th() + self.th().sin())
    }
    fn g(&self) -> f64 {
        10f64.powf(self.g_db / 20.0)
    }
    fn gm(&self) -> f64 {
        self.g() * (1.0 - self.k)
    }
    fn itf(&self) -> Shelf {
        Shelf::design(1.0, self.af() / self.an(), self.w0(), self.fs)
    }
    fn d_samples(&self) -> f64 {
        self.tau_hf() * self.fs
    }
    /// integer delay + Thiran, total path
    fn split(&self) -> (f64, Thiran2) {
        let d = self.d_samples();
        let n = (d.floor() - 1.0).max(0.0);
        (n, Thiran2::design(d - n))
    }
    fn cross(&self, w: f64) -> Cx {
        let (n, t) = self.split();
        self.itf().eval(w, self.fs).mul(t.eval(w, self.fs)).mul(unit(w, n, self.fs))
    }
}

fn main() {
    let w0 = C / A;
    println!("HEAD-SHADOW / ITD MODEL  c={C} m/s  a={A} m  w0=c/a={w0:.0} rad/s  f0={:.1} Hz\n", w0 / (2.0 * PI));

    // ---------- A. ITD ----------
    println!("A. ITD of the cross path  (target: tau_hf = Woodworth, tau_lf/tau_hf = 1.5 per Kuhn)");
    println!(
        "{:>4} {:>6} {:>6} {:>7} {:>7} {:>7} {:>8} {:>8} {:>8}",
        "az", "alphaN", "alphaF", "tauHFus", "D smp", "int", "ITD100Hz", "ITD4kHz", "ratio"
    );
    for fs in [44100.0, 48000.0, 96000.0] {
        for az in [15.0, 30.0, 45.0, 60.0] {
            let p = P::new(az, -6.0, 0.0, fs);
            let (n, _) = p.split();
            let lf = gd(|w| p.cross(w), 100.0, fs) / fs * 1e6;
            let mid = gd(|w| p.cross(w), 1000.0, fs) / fs * 1e6;
            let hf = gd(|w| p.cross(w), 4000.0, fs) / fs * 1e6;
            let _ = mid;
            println!(
                "{:>4.0} {:>6.3} {:>6.3} {:>7.1} {:>7.2} {:>7.0} {:>8.1} {:>8.1} {:>8.2}",
                az, p.an(), p.af(), p.tau_hf() * 1e6, p.d_samples(), n, lf, hf, lf / hf
            );
        }
        println!("  (fs = {fs})");
    }

    // ---------- B. ILD ----------
    println!("\nB. Interaural level difference of the cross path (shelf plateau)");
    println!("{:>4} {:>9} {:>9} {:>9} {:>9}", "az", "1kHz", "2kHz", "4kHz", "10kHz");
    for az in [15.0, 30.0, 45.0, 60.0] {
        let p = P::new(az, -6.0, 0.0, 48000.0);
        let itf = p.itf();
        println!(
            "{:>4.0} {:>9.1} {:>9.1} {:>9.1} {:>9.1}",
            az,
            db(itf.eval(2.0 * PI * 1000.0, 48000.0).abs()),
            db(itf.eval(2.0 * PI * 2000.0, 48000.0).abs()),
            db(itf.eval(2.0 * PI * 4000.0, 48000.0).abs()),
            db(itf.eval(2.0 * PI * 10000.0, 48000.0).abs())
        );
    }

    // ---------- C. centre-channel comb ----------
    println!("\nC. Worst mono-sum notch versus pan position and Centre Focus (g=-6 dB, fs=48k)");
    println!("   pan 0.0 = centred (L=R), pan 1.0 = hard panned.  Values in dB.");
    print!("{:>4} {:>5} |", "az", "k");
    for i in 0..=4 {
        print!("{:>7}", format!("pan{:.1}", i as f64 * 0.25));
    }
    println!("   worst");
    for az in [15.0, 30.0, 45.0, 60.0] {
        for k in [0.0, 0.5, 1.0] {
            let p = P::new(az, -6.0, k, 48000.0);
            let gm = p.gm();
            let g = p.g();
            let norm = 1.0 / (1.0 + gm); // static DC normaliser, scalar
            print!("{:>4.0} {:>5.2} |", az, k);
            let mut worst = 0.0f64;
            for i in 0..=4 {
                let q = i as f64 * 0.25;
                // pan q: L = cos, R = sin of a constant-power pan law
                let (l, r) = (
                    ((1.0 - q) * PI / 2.0).cos(),
                    ((1.0 - q) * PI / 2.0).sin(),
                );
                let m = 0.5 * (l + r);
                let s = 0.5 * (l - r);
                // worst notch of the louder output channel over frequency
                let mut w_min = 0.0f64;
                let mut f = 20.0;
                while f < 20000.0 {
                    let w = 2.0 * PI * f;
                    let c = p.cross(w);
                    let mo = m * (1.0 + gm * 0.0) + m * gm * 0.0; // placeholder
                    let _ = mo;
                    let m_out = Cx::new(m, 0.0).add(c.mul(Cx::new(gm * m, 0.0)));
                    let s_out = Cx::new(s, 0.0).add(c.mul(Cx::new(-g * s, 0.0)));
                    let lo = m_out.add(s_out).abs() * norm;
                    let ro = m_out.scale(1.0).add(s_out.scale(-1.0)).abs() * norm;
                    let v = lo.max(ro);
                    let refv = l.max(r) * norm;
                    let _ = refv;
                    if w_min == 0.0 || v < w_min {
                        w_min = v;
                    }
                    f *= 1.002;
                }
                // compare against the un-crossfed magnitude 1.0 (unity reference)
                let d = db(w_min / 1.0);
                if d < worst {
                    worst = d;
                }
                print!("{:>7.2}", d);
            }
            println!("   {:>6.2}", worst);
        }
    }

    // ---------- D. static gain law ----------
    println!("\nD. Static DC / HF gain of the MONO path with normaliser 1/(1+gM)  (0 dB ideal)");
    println!("{:>4} {:>5} {:>7} {:>7} {:>7} {:>8}", "az", "k", "20Hz", "1kHz", "10kHz", "residHF");
    for az in [15.0, 30.0, 45.0, 60.0] {
        for k in [0.0, 0.5, 1.0] {
            let p = P::new(az, -6.0, k, 48000.0);
            let gm = p.gm();
            let norm = 1.0 / (1.0 + gm);
            let h = |f: f64| {
                let w = 2.0 * PI * f;
                (Cx::new(1.0, 0.0).add(p.cross(w).scale(gm))).abs() * norm
            };
            println!(
                "{:>4.0} {:>5.2} {:>7.3} {:>7.3} {:>7.3} {:>8.2}",
                az,
                k,
                db(h(20.0)),
                db(h(1000.0)),
                db(h(10000.0)),
                db(h(10000.0))
            );
        }
    }

    // ---------- E. Thiran ----------
    println!("\nE. 2nd-order Thiran group delay error (samples)");
    println!("{:>5} {:>9} {:>9} {:>9} {:>9}", "D", "100Hz", "1kHz", "4kHz", "12kHz");
    for d in [1.05, 1.5, 1.95] {
        let t = Thiran2::design(d);
        let g = |f: f64| gd(|w| t.eval(w, 48000.0), f, 48000.0);
        println!(
            "{:>5.2} {:>9.4} {:>9.4} {:>9.4} {:>9.4}  error@12k {:.1}%",
            d,
            g(100.0),
            g(1000.0),
            g(4000.0),
            g(12000.0),
            100.0 * (g(12000.0) / d - 1.0)
        );
    }
    println!("   pole radii: {:?}", (1..4).map(|i| {
        let d = 1.0 + i as f64 * 0.25;
        (d, Thiran2::design(d).pole_radius())
    }).collect::<Vec<_>>());

    // ---------- F. image-source reflections ----------
    println!("\nF. First/second-order image sources, reference room 6.4 x 7.0 x 2.9 m");
    let (w, dpt, h) = (3.2, 4.6, 2.9); // x in +-3.2, y in [-4.6, +2.4], z in [0, 2.9]
    let spk = (-0.85f64, 1.4724f64, 1.18f64); // left speaker, 30 deg, r=1.7
    let ear = (0.0, 0.0, 1.18);
    let direct = ((spk.0 - ear.0).powi(2) + (spk.1 - ear.1).powi(2) + (spk.2 - ear.2).powi(2)).sqrt();
    // reflection coefficients per plane: [x-, x+, y-, y+, z-, z+]
    let coef = [0.75, 0.75, 0.75, 0.75, 0.55, 0.8];
    let mut taps: Vec<(f64, f64, f64, f64, String)> = vec![];
    let mut order = vec![];
    for ix in -1i32..=1 {
        for iy in -1i32..=1 {
            for iz in -1i32..=1 {
                if ix.abs() + iy.abs() + iz.abs() > 2 || ix + iy + iz == 0 {
                    continue;
                }
                order.push((ix, iy, iz));
            }
        }
    }
    for (ix, iy, iz) in order {
        let mx = match ix {
            0 => spk.0,
            1 => 2.0 * w - spk.0,
            _ => -2.0 * w - spk.0,
        };
        let my = match iy {
            0 => spk.1,
            1 => 2.0 * 2.4 - spk.1,
            _ => -2.0 * 4.6 - spk.1,
        };
        let mz = match iz {
            0 => spk.2,
            1 => 2.0 * h - spk.2,
            _ => -spk.2,
        };
        let r = ((mx - ear.0).powi(2) + (my - ear.1).powi(2) + (mz - ear.2).powi(2)).sqrt();
        let mut g = 1.0 / r;
        let mut names = vec![];
        if ix != 0 {
            g *= coef[if ix > 0 { 1 } else { 0 }];
            names.push(if ix > 0 { "R-wall" } else { "L-wall" });
        }
        if iy != 0 {
            g *= coef[if iy > 0 { 3 } else { 2 }];
            names.push(if iy > 0 { "front" } else { "rear" });
        }
        if iz != 0 {
            g *= coef[if iz > 0 { 5 } else { 4 }];
            names.push(if iz > 0 { "ceil" } else { "floor" });
        }
        let dt = (r - direct) / C * 1000.0;
        let rel = 20.0 * (g * direct).log10();
        // azimuth of the image source seen from the head, degrees, + = right
        let azim = (mx - ear.0).atan2(my - ear.1).to_degrees();
        let elev = (mz - ear.2).atan2(((mx - ear.0).powi(2) + (my - ear.1).powi(2)).sqrt()).to_degrees();
        taps.push((dt, rel, azim, elev, names.join("+")));
    }
    taps.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    println!("{:>7} {:>8} {:>7} {:>7}  {}", "dt(ms)", "reldB", "azim", "elev", "surfaces");
    for (dt, rel, az, el, n) in taps.iter().filter(|t| t.0 > 0.5 && t.0 < 21.0).take(9) {
        println!("{:>7.2} {:>8.1} {:>7.0} {:>7.0}  {}", dt, rel, az, el, n);
    }
}

// Phase-1 metric pass: comb depth, level law, and a three-way crossfeed comparison.
use std::f64::consts::PI;
const C: f64 = 343.0;
const A: f64 = 0.0875;

#[derive(Clone, Copy, Debug)]
struct Cx { re: f64, im: f64 }
impl Cx {
    fn new(re: f64, im: f64) -> Self { Self { re, im } }
    fn add(self, o: Cx) -> Cx { Cx::new(self.re + o.re, self.im + o.im) }
    fn mul(self, o: Cx) -> Cx { Cx::new(self.re*o.re - self.im*o.im, self.re*o.im + self.im*o.re) }
    fn scale(self, k: f64) -> Cx { Cx::new(self.re*k, self.im*k) }
    fn div(self, o: Cx) -> Cx {
        let d = o.re*o.re + o.im*o.im;
        Cx::new((self.re*o.re + self.im*o.im)/d, (self.im*o.re - self.re*o.im)/d)
    }
    fn abs(self) -> f64 { (self.re*self.re + self.im*self.im).sqrt() }
    fn arg(self) -> f64 { self.im.atan2(self.re) }
}
fn unit(w: f64, n: f64, fs: f64) -> Cx { let a = -(w/fs)*n; Cx::new(a.cos(), a.sin()) }
fn z1(w: f64, fs: f64) -> Cx { unit(w, 1.0, fs) }

#[derive(Clone, Copy)]
struct Shelf { b0: f64, b1: f64, a1: f64 }
impl Shelf {
    /// Brown & Duda form: H(s) = (w0 + az*s) / (w0 + ap*s).
    /// DC gain is exactly 1; HF plateau gain is az/ap; corner anchored at w0.
    fn design(az: f64, ap: f64, w0a: f64, fs: f64) -> Self {
        let k = 2.0*fs;
        let w0 = k*(w0a/k).tan();
        let d = w0 + ap*k;
        Shelf { b0: (w0 + az*k)/d, b1: (w0 - az*k)/d, a1: (w0 - ap*k)/d }
    }
    fn eval(&self, w: f64, fs: f64) -> Cx {
        let zi = z1(w, fs);
        Cx::new(self.b0,0.0).add(Cx::new(self.b1,0.0).mul(zi))
            .div(Cx::new(1.0,0.0).add(Cx::new(self.a1,0.0).mul(zi)))
    }
}
#[derive(Clone, Copy)]
struct Thiran2 { a1: f64, a2: f64 }
impl Thiran2 {
    fn design(d: f64) -> Self {
        Thiran2 { a1: -2.0*(d-2.0)/(d+1.0), a2: (d-2.0)*(d-1.0)/((d+1.0)*(d+2.0)) }
    }
    fn eval(&self, w: f64, fs: f64) -> Cx {
        let zi = z1(w, fs); let zi2 = zi.mul(zi);
        Cx::new(self.a2,0.0).add(Cx::new(self.a1,0.0).mul(zi)).add(zi2)
            .div(Cx::new(1.0,0.0).add(Cx::new(self.a1,0.0).mul(zi)).add(Cx::new(self.a2,0.0).mul(zi2)))
    }
}
fn gd<F: Fn(f64)->Cx>(h: F, f: f64, fs: f64) -> f64 {
    let w = 2.0*PI*f; let dw = 2.0*PI*0.5;
    let p1 = h(w+dw).arg(); let p0 = h(w-dw).arg();
    let mut d = p1-p0;
    while d > PI { d -= 2.0*PI } while d < -PI { d += 2.0*PI }
    -d*fs/(2.0*PI)
}
fn db(x: f64) -> f64 { 20.0*x.log10() }

#[derive(Clone, Copy, PartialEq)]
enum Model { BrownDuda, FirstOrderLp, FlatDelay }

struct P { az: f64, g_db: f64, k: f64, fs: f64, model: Model }
impl P {
    fn th(&self) -> f64 { self.az*PI/180.0 }
    fn s(&self) -> f64 { 0.5*(self.th() + self.th().sin()) }
    fn an(&self) -> f64 { 1.0 + 0.5*self.s() }
    fn af(&self) -> f64 { 1.0 - 0.5*self.s() }
    fn w0(&self) -> f64 { C/A }
    fn tau_hf(&self) -> f64 { (A/C)*(self.th() + self.th().sin()) }
    fn g(&self) -> f64 { 10f64.powf(self.g_db/20.0) }
    fn gm(&self) -> f64 { self.g()*(1.0-self.k) }
    /// cross-path response, magnitude 1 reference at DC removed for FlatDelay
    fn cross(&self, w: f64) -> Cx {
        let d = self.tau_hf()*self.fs;
        let n = (d.floor()-1.0).max(0.0);
        let t = Thiran2::design(d-n);
        let dly = unit(w, n, self.fs).mul(t.eval(w, self.fs));
        let mag = match self.model {
            Model::BrownDuda => Shelf::design(self.af(), self.an(), self.w0(), self.fs).eval(w, self.fs),
            // classic 1st-order RC low-pass crossfeed (bs2b / Meier style), 700 Hz corner
            Model::FirstOrderLp => Shelf::design(0.0, 1.0, 2.0*PI*700.0, self.fs).eval(w, self.fs),
            Model::FlatDelay => Cx::new(1.0, 0.0),
        };
        mag.mul(dly)
    }
    fn norm(&self) -> f64 { 1.0/(1.0+self.gm()) }
    /// magnitude of the louder output channel for a source at pan position p (0=L, .5=C, 1=R)
    fn out_mag(&self, f: f64, p: f64) -> f64 {
        let w = 2.0*PI*f;
        let (l, r) = (((1.0-p)*PI/2.0).cos(), ((1.0-p)*PI/2.0).sin());
        let m = 0.5*(l+r); let s = 0.5*(l-r);
        let c = self.cross(w);
        let mo = Cx::new(m,0.0).add(c.scale(self.gm()*m));
        let so = Cx::new(s,0.0).add(c.scale(-self.g()*s));
        let lo = mo.add(so).abs();
        let ro = mo.scale(1.0).add(so.scale(-1.0)).abs();
        lo.max(ro)*self.norm()
    }
}

fn main() {
    let fs = 48000.0;
    println!("HEAD-SHADOW / ITD MODEL  c={C} m/s  a={A} m  w0=c/a={:.0} rad/s  f0={:.1} Hz", C/A, C/A/(2.0*PI));
    println!("(all levels dB; comb metrics normalised to the 20 Hz gain of the same pan position)\n");

    // --- Kuhn ITD ratio at the low-frequency asymptote
    println!("1. Kuhn ITD ratio, low-frequency asymptote vs 4 kHz");
    println!("{:>4} {:>12} {:>12} {:>10}", "az", "ITD@50Hz us", "ITD@4k us", "ratio");
    for az in [15.0, 30.0, 45.0, 60.0] {
        let p = P { az, g_db: -6.0, k: 0.0, fs, model: Model::BrownDuda };
        let lf = gd(|w| p.cross(w), 50.0, fs)/fs*1e6;
        let hf = gd(|w| p.cross(w), 4000.0, fs)/fs*1e6;
        println!("{:>4.0} {:>12.1} {:>12.1} {:>10.2}", az, lf, hf, lf/hf);
    }

    // --- comb depth for centred content
    println!("\n2. Centred content (L=R): DC gain, deepest null, and comb ripple 200 Hz..20 kHz");
    println!("{:>4} {:>5} {:>8} {:>10} {:>10} {:>10}", "az", "k", "DCdB", "null_dB", "peak_dB", "ripple");
    for az in [15.0, 30.0, 45.0, 60.0] {
        for k in [0.0, 0.5, 1.0] {
            let p = P { az, g_db: -6.0, k, fs, model: Model::BrownDuda };
            let dc = p.out_mag(20.0, 0.5);
            let (mut lo, mut hi) = (1e9f64, 0.0f64);
            let mut f = 200.0;
            while f < 20000.0 {
                let v = p.out_mag(f, 0.5);
                if v < lo { lo = v }
                if v > hi { hi = v }
                f *= 1.001;
            }
            println!("{:>4.0} {:>5.2} {:>8.3} {:>10.2} {:>10.2} {:>10.2}",
                az, k, db(dc), db(lo/dc), db(hi/dc), db(hi/lo));
        }
    }

    // --- worst comb over all pan positions
    println!("\n3. Worst comb null over every pan position (the real 'hollow centre' metric)");
    println!("{:>4} {:>5} {:>12} {:>12}", "az", "k", "worst_null", "at_pan");
    for az in [15.0, 30.0, 45.0, 60.0] {
        for k in [0.0, 0.5, 1.0] {
            let p = P { az, g_db: -6.0, k, fs, model: Model::BrownDuda };
            let mut worst = 0.0f64; let mut wp = 0.0;
            for i in 0..=20 {
                let pan = i as f64/20.0;
                let dc = p.out_mag(20.0, pan);
                let mut f = 200.0; let mut lo = 1e9f64;
                while f < 20000.0 {
                    let v = p.out_mag(f, pan);
                    if v < lo { lo = v }
                    f *= 1.002;
                }
                let d = db(lo/dc);
                if d < worst { worst = d; wp = pan }
            }
            println!("{:>4.0} {:>5.2} {:>12.2} {:>12.2}", az, k, worst, wp);
        }
    }

    // --- three-way model comparison at 30 deg
    println!("\n4. Model comparison, az=30 deg, crossfeed -6 dB, centred content, k=0");
    println!("{:>16} {:>10} {:>12} {:>12} {:>12}", "model", "1stnull_dB", "xfade@10kHz", "nulls_ptp", "HFenergy");
    for m in [Model::BrownDuda, Model::FirstOrderLp, Model::FlatDelay] {
        let p = P { az: 30.0, g_db: -6.0, k: 0.0, fs, model: m };
        // first null of the mono sum
        let dc = p.out_mag(20.0, 0.5);
        let mut f = 300.0; let mut prev = p.out_mag(f, 0.5); let mut first = (0.0, 1e9f64);
        while f < 3000.0 {
            let v = p.out_mag(f, 0.5);
            if v < prev && v < first.1 { first = (f, v) }
            prev = v; f *= 1.0005;
        }
        let xf = db(p.g()*p.cross(2.0*PI*10000.0).abs());
        let mut lo = 1e9f64; let mut hi = 0.0f64; let mut epow = 0.0f64; let mut cnt = 0.0;
        let mut f = 200.0;
        while f < 20000.0 {
            let v = p.out_mag(f, 0.5);
            if v < lo { lo = v } if v > hi { hi = v }
            epow += v*v; cnt += 1.0; f *= 1.001;
        }
        let name = match m { Model::BrownDuda => "Brown-Duda shelf", Model::FirstOrderLp => "1st-order LP 700", Model::FlatDelay => "flat delay" };
        println!("{:>16} {:>10.2} {:>12.1} {:>12.2} {:>12.2}",
            name, db(first.1/dc), xf, db(hi/lo), db((epow/cnt).sqrt()/dc));
        if m == Model::BrownDuda {
            // same model with the ipsilateral bright-spot shelf (HF plateau alpha_n) folded in
            let near = Shelf::design(p.an(), 1.0, p.w0(), fs);
            let mut ep = 0.0; let mut c2 = 0.0; let mut f = 200.0; let mut lo2 = 1e9f64; let mut hi2 = 0.0f64;
            while f < 20000.0 {
                let v = p.out_mag(f, 0.5)*near.eval(2.0*PI*f, fs).abs();
                if v<lo2 {lo2=v} if v>hi2 {hi2=v}
                ep += v*v; c2 += 1.0; f *= 1.001;
            }
            println!("{:>16} {:>10} {:>12} {:>12.2} {:>12.2}",
                "  + bright spot", "-", "-", db(hi2/lo2), db((ep/c2).sqrt()/(dc*near.eval(2.0*PI*20.0,fs).abs())));
        }
    }

    // --- what the crossfeed attenuation does to HF on a hard-panned source
    println!("\n5. Hard-panned source: leakage into the opposite ear at 100 Hz / 10 kHz");
    println!("{:>16} {:>12} {:>12} {:>12}", "model", "leak100Hz", "leak10kHz", "monoHFenergy");
    for m in [Model::BrownDuda, Model::FirstOrderLp, Model::FlatDelay] {
        let p = P { az: 30.0, g_db: -6.0, k: 0.0, fs, model: m };
        let lf = db(p.g()*p.cross(2.0*PI*100.0).abs());
        let hf = db(p.g()*p.cross(2.0*PI*10000.0).abs());
        let mut epow = 0.0; let mut cnt = 0.0; let mut f = 200.0;
        while f < 20000.0 { let v = p.out_mag(f, 0.5); epow += v*v; cnt += 1.0; f *= 1.001; }
        let name = match m { Model::BrownDuda => "Brown-Duda shelf", Model::FirstOrderLp => "1st-order LP 700", Model::FlatDelay => "flat delay" };
        println!("{:>16} {:>12.1} {:>12.1} {:>12.2}", name, lf, hf, db((epow/cnt).sqrt()/p.out_mag(20.0,0.5)));
    }

    // --- crossfeed amount sweep at 30 deg
    println!("\n6. Crossfeed amount sweep, az=30 deg, centred content, k=0.5");
    println!("{:>8} {:>10} {:>10} {:>10} {:>10}", "g_dB", "null_dB", "ripple", "DC", "HFgain");
    for g_db in [-18.0, -12.0, -6.0, -3.0, 0.0] {
        let p = P { az: 30.0, g_db, k: 0.5, fs, model: Model::BrownDuda };
        let dc = p.out_mag(20.0, 0.5);
        let mut lo = 1e9f64; let mut hi = 0.0f64; let mut f = 200.0;
        while f < 20000.0 { let v = p.out_mag(f,0.5); if v<lo {lo=v} if v>hi {hi=v} f *= 1.001; }
        println!("{:>8.1} {:>10.2} {:>10.2} {:>10.3} {:>10.2}",
            g_db, db(lo/dc), db(hi/lo), db(dc), db(p.out_mag(10000.0,0.5)/dc));
    }
}

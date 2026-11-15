use pairfield::dsp::{cross::CrossPath, woodworth_delay_s};

const FS: f64 = 48000.0;

fn tf(ir: &[f32], f: f64) -> (f64, f64) {
    let w = -2.0 * std::f64::consts::PI * f / FS;
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

fn gd(ir: &[f32], f: f64) -> f64 {
    let (r1, i1) = tf(ir, f + 0.5);
    let (r0, i0) = tf(ir, f - 0.5);
    let mut d = i1.atan2(r1) - i0.atan2(r0);
    while d > std::f64::consts::PI {
        d -= 2.0 * std::f64::consts::PI;
    }
    while d < -std::f64::consts::PI {
        d += 2.0 * std::f64::consts::PI;
    }
    -d * FS / (2.0 * std::f64::consts::PI)
}

fn main() {
    let mut c = CrossPath::new(woodworth_delay_s(60.0) * FS);
    for deg in [15.0, 15.5, 16.0, 16.4, 16.5, 17.0, 20.0, 25.0] {
        c.reset();
        c.set(deg, FS);
        let mut ir = vec![0.0f32; 1024];
        ir[0] = 1.0;
        for x in ir.iter_mut() {
            *x = c.step(*x);
        }
        let nominal = woodworth_delay_s(deg) * FS;
        println!("{deg:>5.1} deg  nominal {nominal:>7.3}  measured@16k {:>7.3}  measured@8k  {:>7.3}  peak idx {:?}",
            gd(&ir, 16000.0), gd(&ir, 8000.0),
            ir.iter().enumerate().max_by(|a,b| a.1.abs().partial_cmp(&b.1.abs()).unwrap()).map(|(i,_)| i));
    }
}

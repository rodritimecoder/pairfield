# ARCHITECTURE.md — Pairfield: Binaural Crossfeed & Virtual-Monitor VST3

Status: **Phase 1 (research + DSP design) — for review before implementation.**
Every number in this document was produced by a standalone numeric probe
(`/tmp/opencode/verify/v3.rs`, `v2.rs`), not estimated. Reproduce with `rustc -O v3.rs && ./v3`.

---

## 0. Design thesis

Crossfeed is usually implemented as *"add a low-passed, delayed copy of the opposite
channel"* (Bauer/`bs2b`, Meier). That formulation has two structural faults:

1. **The cross path's magnitude is a guess.** A 1st-order RC low-pass keeps falling at
   6 dB/oct forever. Measured (Table 4): with −6 dB of crossfeed at 30°, an RC-700 Hz
   cross path is **−30.5 dB at 10 kHz**, i.e. the HF image is un-glued and the top end
   decouples from the bottom. That is the "dark / disconnected" complaint.
2. **The comb is uncontrolled.** Direct + delayed copy always combs. Depth is set by how
   much magnitude the cross path still has where the delay phase first reaches 180°.

This design replaces (1) with the **Brown & Duda rigid-sphere head-shadow shelf**, whose
whole point is a *plateau*: the cross path stops falling and settles at a shadowed level
(−10.5 dB at 10 kHz instead of −30.5 dB). And it makes (2) an **explicit, measured budget**
the user controls, instead of an accident.

Measured against the alternatives at 30° / −6 dB / centred content (Table 4):

| cross path | first null | HF crossfeed @10 kHz | comb p-p | HF energy vs DC |
|---|---|---|---|---|
| flat delay (worst case) | **−9.57 dB** | −6.0 dB | 9.57 dB | −2.13 dB |
| 1st-order RC 700 Hz (`bs2b`-class) | −6.04 dB | **−30.5 dB** | 5.57 dB | −3.14 dB |
| **Brown & Duda shelf (this design)** | −6.82 dB | −10.5 dB | 6.47 dB | −2.73 dB |

It is not the shallowest comb — a low-pass is, by starving the cross path of HF. It is the
only one that is *simultaneously* shallow and HF-coherent. And the comb depth is then dialled
by parameter 3 (§5), which no crossfeed plugin currently exposes.

---

## 1. Constants and notation

| symbol | value | meaning |
|---|---|---|
| `c` | 343 m/s | speed of sound |
| `a` | 0.0875 m | head radius (Brown & Duda's cited adult mean) |
| `ω0` | `c/a` = 3920 rad/s | **`f0 = 623.9 Hz`** — the head's acoustic corner |
| `θ` | 15…60° | virtual speaker azimuth from centre (param 1) |
| `g` | 10^(dB/20), −18…0 dB | crossfeed gain (param 2) |
| `κ` | 0…1 | Centre Focus (param 3) |
| `f_s` | 44.1…192 kHz | sample rate |

---

## 2. Layer 1 — spherical head diffraction (Brown & Duda 1998)

Brown & Duda, *A Structural Model for Binaural Sound Synthesis*, IEEE TSAP 6(5):476–488,
1998, eqs. (1)–(5). Rayleigh's rigid-sphere solution is approximated by a single
pole/zero pair whose zero moves with the angle of incidence `θ_inc`:

```
H(s, θ) = (ω0 + α(θ)·s) / (ω0 + s)                        … (BD-3)
```

with the pole pinned at `ω0 = c/a` and

```
α(θ) = 1.05 + 0.95·cos( min(θ_inc, 150°) · π / 150° )       … (BD-5)
```

so `α(0°) = 2.0` (**+6 dB high-frequency plateau** — the "bright spot"), `α(90°) = 1.05`,
`α(150°) = 0.1` (−20 dB). DC gain is exactly 1 for every `α`: this filter is a *shelf*, not
a low-pass. The plateau is the entire reason it does not sound muffled.

### 2.1 Re-parameterisation for a symmetric speaker pair (the one real deviation)

For a pair at ±θ, the near ear sees incidence `90°−θ` and the far ear `90°+θ`, so (BD-5)
gives two different `α` and — critically — the *ratio* of the two filters carries both the
ILD **and** a low-frequency group delay. Numerically applying (BD-5) directly to the pair
overshoots Kuhn's measured ITD ratio by ~35% (the paper's `α` fit was calibrated on-axis, not
on a ±30° pair).

Instead we derive the α pair from **one geometric scalar**, and it lands on Kuhn's ratio
exactly.

```
s(θ) = 0.5·(θ + sin θ)          [rad, θ in radians]
α_near = 1 + s/2                    α_far = 1 − s/2         … (PF-1)
```

Why: `α_near − α_far = s`, and the interaural filter's DC group delay is

```
τ_ITF(0) = (α_near − α_far)/ω0 = (a/c)·0.5·(θ + sin θ) = 0.5·τ_HF       … (PF-2)
```

i.e. **exactly half** of the Woodworth delay, which is Kuhn's 3:2 low-to-high ITD ratio.
Table 1 confirms 1.51 measured at every angle. So the frequency-dependent ITD falls out of
the shadow filter's own phase — no extra allpass network is needed for it.

| θ | s | α_near | α_far | τ_HF (Woodworth) | ILD @4 kHz | ILD @10 kHz |
|---|---|---|---|---|---|---|
| 15° | 0.259 | 1.1297 | 0.8703 | 132.8 µs | −2.2 dB | −2.3 dB |
| 30° | 0.512 | 1.2559 | 0.7441 | 261.1 µs | −4.4 dB | −4.5 dB |
| 45° | 0.746 | 1.3731 | 0.6269 | 380.7 µs | −6.5 dB | −6.8 dB |
| 60° | 0.957 | 1.4783 | 0.5217 | 488.1 µs | −8.4 dB | −9.0 dB |

(The `bs2b`-class ILD would be a straight low-pass; these values track published
KEMAR bare-head ILDs to within a couple of dB.)

### 2.2 The only filter primitive

Everything in the plugin is this one 5-multiply first-order shelf. Analog prototype
`H(s) = (ω0 + α_z·s)/(ω0 + α_p·s)`, bilinear with prewarping:

```
k  = 2·f_s ,      ω0' = k·tan(ω0/(2·f_s))
d  = ω0' + α_p·k
b0 = (ω0' + α_z·k)/d
b1 = (ω0' − α_z·k)/d
a1 = (ω0' − α_p·k)/d          … (PF-3)
```

`H(z) = (b0 + b1 z⁻¹)/(1 + a1 z⁻¹)`.

Properties (asserted in unit tests):
* DC gain `H(1) = 1` **exactly**, for every `(α_z, α_p)`, at every `f_s`;
* HF plateau `H(−1) = α_z/α_p`;
* pole at `z = (α_p·k − ω0')/(α_p·k + ω0')`, strictly inside the unit circle for all
  positive parameters ⇒ unconditionally stable;
* `|a1| < 1`, no branch switching, no coefficient blow-up at low `f_s`.

Used as: `ITF = shelf(α_far, α_near, ω0)` — DC 1, plateau `α_far/α_near`.

---

## 3. Layer 2 — ITD: Woodworth bulk delay + Kuhn's ratio

Woodworth's spherical-head path difference (Brown & Duda eq. 2):

```
τ_HF(θ) = (a/c)·(θ + sin θ)            … (W)
D = τ_HF · f_s  samples
```

Kuhn (JASA 62(1):157–167, 1977) measured that low-frequency ITD is 1.5× the high-frequency
ITD (LF asymptote `3a/c·sin θ`, HF `2a/c·sin θ`). That factor is delivered by (PF-2) — the
cross path is `ITF(z)·z^−D`, and the shelf's phase *is* the extra low-frequency delay.

**Measured (Table 1, `fs = 48 kHz`):**

| θ | ITD @50 Hz | ITD @4 kHz | ratio |
|---|---|---|---|
| 15° | 197.9 µs | 131.3 µs | **1.51** |
| 30° | 389.1 µs | 258.0 µs | **1.51** |
| 45° | 567.2 µs | 375.8 µs | **1.51** |
| 60° | 727.0 µs | 481.1 µs | **1.51** |

Measured from the *implemented* filter's impulse response (`examples/measure.rs`), not from
the closed form. `τ_HF` reads 1.2–1.4 % **below** Woodworth at 4 kHz because the shelf itself
contributes a small group *advance* at HF (−0.15 samples at 30°) — which is precisely the
frequency dependence being modelled. At the 50 Hz asymptote the bulk delay is exact.

Identical in microseconds at 44.1 / 48 / 96 / 192 kHz (spread 0.07 µs): the delay is
specified in seconds and sampled at `f_s`.

### 3.1 Fractional delay: 3rd-order Thiran

`D = n + d`, `n = floor(D) − 2`, `d = clamp(D − n, 2.001, 2.999)`:

```
a1 = −3(d−3)/(d+1)
a2 =  3(d−3)(d−2) / ((d+1)(d+2))
a3 = −(d−3)(d−2)(d−1) / ((d+1)(d+2)(d+3))
H(z) = (a3 + a2 z⁻¹ + a1 z⁻² + z⁻³)/(1 + a1 z⁻¹ + a2 z⁻² + a3 z⁻³)   … (T3)
```

Maximally flat group delay at DC and — crucially — **allpass**, so it cannot colour the comb.
`|C|` matches the shelf magnitude to <0.02 dB at 100 Hz…15 kHz (asserted in
`tests/dsp.rs::thiran_is_allpass`).

Third order, not second, for a specific reason the tests found: *every* Thiran order
degenerates into a marginal pole at integer delays (`d = 1, 2, 3`), and with a 2nd-order
network the fractional part lives in `[1,2]` so **both** ends are degenerate. Centring
`d ∈ [2,3]` keeps the poles away from `z = −1` for all but 0.1 % of the sweep, and the
clamp costs ≤0.001 samples (0.02 µs), so the *total* delay `n + d` stays continuous as the
knob moves — which is what stops the angle control from zippering
(`tests/dsp.rs::delay_tracks_angle_continuously` sweeps 15°→60° in 0.1° steps and asserts
the measured step matches the nominal step to 0.02 samples).

Residual: within 0.1 % of an integer `D` the network leaves a −60 dB Nyquist-frequency tail
decaying over ~800 samples. Below the noise floor of any master, but it is why the harness
uses 8192-sample impulse responses instead of 1024.

Zero look-ahead ⇒ **0 samples of latency**.

Reflection taps use linear interpolation instead: they are ≥18 dB down and wall absorption
is itself a HF loss, so the interpolator's own −2 dB at 10 kHz is physically the right
colour. Fractional delay quality matters only where the ITD is a cue.

---

## 4. Layer 3 — M/S kernel (comb budget instead of comb accident)

Any crossfeed is, in Mid/Side coordinates, a *gain asymmetry*:

```
L' = L + g·C·R          M' = M·(1 + g·C)       C = ITF(z)·z^−D
R' = R + g·C·L    ⟺     S' = S·(1 − g·C)
```

The comb lives equally in M and S — that is why it cannot be removed in M/S *per se*. What
M/S buys is **where** the comb is allowed to live. With Centre Focus `κ`:

```
M' = M·(1 + (1−κ)·g·C)
S' = S·(1 − g·C)
L' = M' + S'  ·  (1/(1+(1−κ)g))          R' = M' − S'  ·  (1/(1+(1−κ)g))
```

* `κ = 0` → exact physics: M and S see the same crossfeed.
* `κ = 1` → **the M path is untouched ⇒ dead-centre mono is bit-identical to bypass.**
  Measured: deepest null 0.00 dB (Table 2). Not "≈0" — structurally zero.
* `κ = 0.5` → the two extremes meet in the middle.

**Measured (Table 2/3, `g = −6 dB`), deepest null relative to that pan position's passband:**

| θ | κ=0.0 | κ=0.5 | κ=1.0 |
|---|---|---|---|
| 15° | −7.80 dB | −3.82 dB | 0.00 dB |
| 30° | −6.82 dB | −3.43 dB | 0.00 dB |
| 45° | −6.39 dB | −3.25 dB | 0.00 dB |
| 60° | −6.31 dB | −3.22 dB | 0.00 dB |

The worst case over *all* pan positions occurs at pan = centre (Table 3) — so protecting the
centre protects the worst case. Crossfeed-amount sweep at 30°/κ=0.5 (Table 6): −18 dB → −0.88
dB null, −6 dB → −3.43 dB, 0 dB → −6.80 dB. The control degrades gracefully.

### 4.1 Normalisation (exact, scalar — deliberately not a filter)

An earlier revision normalised with a matching shelf to flatten DC *and* HF. That was
**wrong and measurable**: two first-order shelves sharing a pole compose to a *second*-order
shape, producing a −2.6 dB hole at 1.7 kHz — a worse artifact than the one it removed.
Rejected. The shipped law is a scalar:

```
G_dc = 1 / (1 + (1−κ)·g)        ⇒ mono DC gain is exactly 0 dB (verified: ±0.01 dB at 20 Hz)
```

The HF asymptote is then allowed to sit where physics puts it: mono HF energy
−2.73 dB vs DC for κ=0 and +1.4 dB if the ipsilateral bright-spot shelf is folded in
(Table 4, row 2). Residual program-dependent level is the RMS matcher's job (§6).

`G_dc` is applied **before** the RMS detector, so the detector only measures what the matcher
actually has to fix. Getting this order wrong leaves a residual `G_dc²` error — measured as a
−1.94 dB level mismatch on every content type before it was fixed, and ±0.01 dB after.

---

## 5. Layer 4 — Early reflections (precedence, not reverb)

Image-source model of a reference control room **6.4 × 7.0 × 2.9 m**, head at (0,0,1.18),
speakers at `r = 1.7 m`, `h = 1.18 m`. Reflection coefficients: walls 0.75, ceiling 0.80,
floor 0.55. First- and second-order images computed; those inside the 6–20 ms precedence
window, ordered earliest-first, with gain `(1/r)·Π|ρ|` relative to the direct path:

| # | Δt (ms) | rel. level | azimuth | elev | surfaces |
|---|---|---|---|---|---|
| 1 | 6.23 | −9.0 dB | −30° | +64° | ceiling |
| 2 | 9.22 | −13.6 dB | −14° | +45° | front + ceiling |
| 3 | 11.78 | −13.1 dB | −75° | 0° | left wall |
| 4 | 13.14 | −18.9 dB | −75° | −22° | left wall + floor |
| 5 | 16.61 | −15.3 dB | +79° | 0° | right wall |
| 6 | 18.30 | −18.4 dB | +65° | 0° | right wall + front |

All tap levels are then offset so tap 1 sits at the Room knob value (default −22 dB) ⇒
individual taps land at −22 … −31 dB, matching the brief's −18…−26 dB window.

**Deliberate omission:** the true floor bounce (3.52 ms, −9.9 dB) and the front-wall bounce
(5.06 ms) are *inside* the precedence window (<6 ms). They sum perceptually with the direct
sound, so at −22 dB they add combing with no externalisation benefit. They belong to the
direct sound, not to this layer. The brief's 6–20 ms window is correct and we honour it.
Distinguishing feature vs. every crossfeed plugin: taps 5 and 6 arrive at **+79° / +65°**,
i.e. *behind* the interaural axis — that rear-quarter energy is what produces the "outside
the head" percept (and it is the only cue here that cannot be faked by a shelving filter).

### 5.1 Per-tap head-shadow, 2×2 symmetric per tap

Tap `i` at azimuth `az_i` (built for the left speaker; the right is its mirror). Ear
incidence is measured from the interaural axis (ears at ±90°), clamped to 150°:

```
α_same = γ · α( clamp(|az_i + 90°|) )        α_cross = γ · α( clamp(|az_i − 90°|) )
γ = 0.6   (broadband wall/air absorption, 1st-order HF loss folded into α)
```

then, with `z_i` = the tap's fractional delay,

```
L' += g_i·( shelf(α_same,1)·z_i·L + shelf(α_cross,1)·z_i·R )
R' += g_i·( shelf(α_cross,1)·z_i·L + shelf(α_same,1)·z_i·R )
```

12 shelves + 12 interpolated reads total. Symmetric by construction (the 2×2 is
`[[a,b],[b,a]]`), so it cannot shift the image.

---

## 6. Layer 5 — Air trim and RMS auto-gain

**Air** (param 5): one more `shelf(α_z, 1)` with `α_z = 10^(dB/20)`, anchored at 3.5 kHz,
±3 dB. Purely tonal; it is the only non-physical stage (and the honest place for a
"bright-spot" preference).

**Auto-gain** (param 6): two 300 ms EMAs (per-sample one-pole, `a = 1 − exp(−1/(τ·f_s))`) of
`Σ_ch x²` for the dry and wet signals,

```
G_dyn = clamp( sqrt(P_dry / P_wet), −9 dB, +9 dB ),  smoothed 30 ms
wet *= G_dc · G_dyn ;   bypass ⇒ output = dry, exactly, gain frozen
```

Rounding: the analytic `G_dc` already pins mono DC to 0 dB, so `G_dyn` only has to absorb
program-dependent residue (panned content, reflections, the leq of the comb). Measured
spread across content types (Table, §9) is −0.3…−1.1 dB for mono and −3…−5 dB for
hard-panned before matching — all inside the clamp. Gain is applied to the **wet path only**,
so Bypass stays a true, bit-identical bypass and A/B is level-matched by construction.
Freeze when `P_wet < 1e-12` so silence cannot ratchet the gain.

---

## 7. Signal flow

```
 in L ──┬─► shelf(α_far,α_near;ω0) ─► z^-D (3rd-order Thiran) ─► c_L ─┐
        │                                                           │
 in R ──┴─► shelf(α_far,α_near;ω0) ─► z^-D (3rd-order Thiran) ─► c_R ─┤
                                                                      │
 out L = L + ½(1−κ)g·(c_L − c_R) + g·c_R        (+ room)  (+ air)  ×G_dc  ×G_dyn
 out R = R + ½(1−κ)g·(c_R − c_L) + g·c_L        (+ room)  (+ air)  ×G_dc  ×G_dyn
                                    ▲
      6-tap bus ────────────────────┘  from L and R directly: 12 shelves + 12 interp reads
```

Written out that way — two filtered, delayed copies of the *input channels* rather than an
explicit M/S matrix — because it needs two `C` instances either way, and this form reduces to
`L + g·c_R` at κ = 0 (the physical answer) without any sum/difference arithmetic on the audio
before the delay line. The κ term is the M/S asymmetry of §4 in disguise: for L = R the two
coefficients collapse to `g·c`, which is `M(1 + g_m C)` exactly.

State per instance: two delay lines sized from `f_s` at `prepare()` (≤ 256 samples for
ITD at 60°/192 kHz, ≤ 8192 for the reflection taps), 14 first-order shelf states, one RMS
pair. **All buffers allocated in `activate()`**; `process()` does no allocation, takes no
locks, and grows no `Vec` — measured with a counting global allocator: **0 allocations
across 20 000 blocks with the angle and Centre Focus moving every block**. Parameter smoothing is `Smoother` over 30–50 ms; coefficients are
refreshed every 32 samples from the smoothed value (≈ 1/75 of a knob sweep per step), so a
fast drag cannot zipper. Denormals: nice-plug flushes them around `process()`, and no path
in this plugin relies on that — every feedback path is a single pole inside the unit circle.

---

## 8. Parameters

| # | name | range | default | smoothing |
|---|---|---|---|---|
| 1 | Speaker Angle | 15…60° | 30° | 50 ms |
| 2 | Crossfeed Amount | −18…0 dB | −6 dB | 30 ms |
| 3 | Centre Focus | 0…100 % | **50 %** | 30 ms |
| 4 | Room Externalisation | on/off + −30…−12 dB | on, −22 dB | 30 ms |
| 5 | Air / HF Trim | −3…+3 dB | 0 dB | 30 ms |
| 6 | Auto-Gain Match | on/off | on | — |
| 7 | Bypass | on/off | off | 5 ms crossfade |

Latency reported to host: **0**.

---

## 9. Verification — what was actually run

Everything below is reproducible from this repository. Results are from this build.

**`cargo test --release`** — 11 tests, all passing:

| test | claim it pins down |
|---|---|
| `head_shadow_alpha_anchors` | α(0) = 2.0, α(≥150) = 0.1, monotone |
| `alpha_pair_gives_kuhn_ratio` | `(α_n−α_f)/ω0 == ½·τ_Woodworth` to 1e-12, every angle |
| `shelf_dc_gain_is_exactly_unity` | `H(DC) = 1` to 1e-9 for 5 sample rates × 7 α_z × 5 α_p |
| `delay_tracks_angle_continuously` | 451 angles, 0.1° apart: delay error < 0.15 samples, step matches nominal to 0.02 |
| `thiran_is_allpass` | the delay network cannot colour the comb (< 0.02 dB deviation) |
| `centre_focus_one_is_bit_transparent_for_mono` | κ = 100 % ⇒ L = R passes through to 6e-8 |
| `bypass_is_bit_exact` | `to_bits()` equality on 1024 frames |
| `no_nans_across_the_parameter_space` | 36 parameter corners at 0.9 FS |
| `reset_clears_state`, `shelf_pole_is_stable`, `woodworth_reference_values` | … |

**`cargo run --release --example measure`** — 32 checks, all passing. Impulse-response
measurement of the real DSP core: shelf DC/plateau against α_far/α_near, ITD in samples
against Woodworth, the 3:2 ratio, sample-rate independence, the full pan × frequency null
scan, κ = 1 exactness, bit-exact bypass, RMS match error, NaN and subnormal audit. The
measured null depths equal the closed-form predictions to **0.00 dB** at every angle and every
κ — two independent code paths agreeing.

Auto-gain residuals after 4 s of noise: mono −0.00 dB, hard-panned −0.01 dB. Heap
activity during processing: **0 allocations / 20 000 blocks**.

**Steinberg VST 3.8.1 `validator`,** built from `github.com/steinbergmedia/vst3sdk` at
`3cdf9ca5`, run against `target/bundled/pairfield.vst3`:

```
-------------------------------------------------------------
Result: 47 tests passed, 0 tests failed
-------------------------------------------------------------
```

…including *Variable Block Size*, *Process Format*, *Parameter Bypass persistence*, and the
block- and sample-accurate parameter-change accuracy suites.

**Not verified by me:** how it sounds, and how the editor looks on your screen. Both need you.

## 10. Decisions taken at the Phase-1 review gate

| question | decision |
|---|---|
| Centre Focus default | **50 %** — minimax across pan positions (−3.4 dB worst null instead of −6.8 dB) |
| Ipsilateral bright spot | **out of the signal path**; it is common-mode, and the Air knob is where a bright-spot preference belongs |
| Room control shape | **two parameters** (on/off + level) rather than one 19-position stepped knob |
| Verification gate | unit tests + offline harness + Steinberg validator |
| Formats / name | **VST3 + CLAP**, `Pairfield` |

## 11. Deviations, risks, open questions

1. `α` is re-parameterised from `s(θ)` (PF-1) instead of (BD-5). Rationale §2.1; it buys the
   exact 3:2 ITD. If you want literal Brown–Duda, it is one function swap — but the ITD
   ratio then reads ~1.35 and the ILD jumps to −9 dB at 30°.
2. The bright spot is **out** of the signal path by default (it is common-mode, so it costs
   no interaural information; folding it in adds +1.4 dB HF energy, Table 4).
3. Reflection delays are fixed (not re-derived when θ changes). Recomputing them per
   parameter change needs coefficient smoothing on 12 filters and buys <0.5 dB of realism.
4. `κ` trades centre-comb for side-comb rather than removing comb. Removing comb entirely is
   impossible for any linear time-invariant "add a delayed copy" — the value here is that the
   budget is exposed and measured.
5. Not implemented: torso/shoulder (Algazi 2001) and pinna models. Both need per-subject
   parameters and would fight a mastering tool's neutrality.
6. Mono (1-channel) tracks are passed through untouched: the plugin declares 2-in/2-out.
   A mono track has no interaural cues to synthesise, but the room layer could still apply.
   Deliberately not done rather than half-done.
7. The editor has not been opened on a machine with a working audio-capable host during
   development. The drawing code is wired to the live parameters and compiles, but the first
   human look at it is yours.

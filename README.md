# Pairfield

Binaural crossfeed and virtual-monitor stage for open-back headphones. Stereo in, stereo
out, **zero latency**, VST3 + CLAP, Rust.

`ARCHITECTURE.md` is the design document and the reference for the maths. This file is how to
build it and how to check that it does what the document says.

## Why it is different

Most crossfeed is "add a low-passed, delayed copy of the other channel". That has two
structural faults: the cross path's magnitude is a guessed RC curve that keeps falling at
6 dB/oct (measured: −30.5 dB of crossfeed left at 10 kHz, which is where the top end
disconnects from the bottom), and the comb filter is whatever falls out.

Pairfield replaces the cross path with the Brown & Duda rigid-sphere head-shadow **shelf**,
whose high-frequency plateau is the point, and derives the α pair from one geometric scalar
so that the filter's own phase delivers Kuhn's 3:2 low/high ITD ratio exactly (measured
1.51). The remaining comb is a measured budget you control with **Centre Focus**, which trades
a deep null on the centre for a shallow one on the sides — and at 100 % makes dead-centre
mono bit-identical to bypass.

## Install from a release

Prebuilt bundles are attached to every tagged release (Linux x86_64, macOS universal,
Windows x86_64):

https://github.com/rodritimecoder/pairfield/releases

Unpack into `~/.vst3` and `~/.clap` (Linux), `~/Library/Audio/Plug-Ins/{VST3,CLAP}` (macOS),
or `%COMMONPROGRAMFILES%\VST3` and `%LOCALAPPDATA%\Programs\Common\CLAP` (Windows).

## Build

Requires Rust 1.88+, a C toolchain, and (for the editor) an OpenGL-capable desktop.

```sh
cargo nice-plug bundle pairfield --release
```

produces:

```
target/bundled/pairfield.vst3   →  ~/.vst3/
target/bundled/pairfield.clap   →  ~/.clap/
```

## Verify

```sh
cargo test --release                      # 11 DSP invariants
cargo run --release --example measure     # 33 impulse-response checks against the spec
```

`measure` exits non-zero if anything fails. It prints the measured ITD, head-shadow
plateau, comb depth per pan position, level-matching error, and allocation/NaN audits.

Plugin-format validation (Steinberg's validator, built separately):

```sh
git clone --depth 1 --shallow-submodules --recurse-submodules \
    https://github.com/steinbergmedia/vst3sdk
cmake -S vst3sdk -B vst3sdk-build -DCMAKE_BUILD_TYPE=Release \
      -DSMTG_ENABLE_VST3_PLUGIN_EXAMPLES=OFF -DSMTG_ENABLE_VSTGUI_SUPPORT=OFF
cmake --build vst3sdk-build --target validator -j
vst3sdk-build/bin/Release/validator target/bundled/pairfield.vst3
```

The editor target in that SDK needs X11 headers; drop
`public.sdk/samples/vst-hosting/editorhost` if you do not have them.

## Controls

| control | range | default | what it does |
|---|---|---|---|
| Speaker Angle | 15–60° | 30° | virtual speaker azimuth; drives the ITD and the head-shadow pole/zero pair |
| Crossfeed Amount | −18…0 dB | −6 dB | interaural level difference at low frequency |
| Centre Focus | 0–100 % | 50 % | 100 % leaves the centre untouched (bit-exact mono) |
| Room | on/off, −30…−12 dB | on, −22 dB | six early reflections, 6–19 ms, rear-quarter weighted |
| Air / HF Trim | −3…+3 dB | 0 dB | shelf at 3.5 kHz |
| Auto-Gain Match | on/off | on | RMS-matched wet path, so bypass is an honest A/B |
| Bypass | on/off | off | bit-exact |

## Limits

- 2-in / 2-out only; mono tracks pass through untouched.
- No individualised HRTF. The model is a sphere of radius 8.75 cm — deliberately, because
  a mastering tool has to be neutral for everyone.
- Torso, shoulder and pinna are not modelled (`ARCHITECTURE.md` §11).

## Licence

ISC — see [`LICENSE`](LICENSE). Use it, fork it, ship it in a commercial build; keep the
copyright notice. The release bundles are built from this source by
`.github/workflows/release.yml`.

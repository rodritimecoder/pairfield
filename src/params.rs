//! Exposed parameters. Ranges and defaults are the ones documented in ARCHITECTURE.md section 8.

use nice_plug::prelude::*;
use std::sync::Arc;

pub const MIN_ANGLE: f32 = 15.0;
pub const MAX_ANGLE: f32 = 60.0;

fn db_string(v: f32) -> String {
    format!("{v:+.1} dB")
}

#[derive(Params)]
pub struct PairfieldParams {
    /// Virtual speaker azimuth from centre. Drives both the Wooordworth ITD and the
    /// Brown & Duda head-shadow pole/zero pair.
    #[id = "angle"]
    pub speaker_angle: FloatParam,

    /// Crossfeed interaural level difference at DC.
    #[id = "amount"]
    pub crossfeed: FloatParam,

    /// 0 % = full physical M/S crossfeed, 100 % = the M path is left alone so dead-centre
    /// mono is bit-identical to bypass.
    #[id = "focus"]
    pub center_focus: FloatParam,

    #[id = "room_on"]
    pub room_on: BoolParam,

    #[id = "room_lvl"]
    pub room_level: FloatParam,

    #[id = "air"]
    pub air: FloatParam,

    #[id = "autogain"]
    pub auto_gain: BoolParam,

    #[id = "bypass"]
    pub bypass: BoolParam,
}

impl Default for PairfieldParams {
    fn default() -> Self {
        Self {
            speaker_angle: FloatParam::new(
                "Speaker Angle",
                30.0,
                FloatRange::Linear {
                    min: MIN_ANGLE,
                    max: MAX_ANGLE,
                },
            )
            .with_smoother(SmoothingStyle::Linear(50.0))
            .with_step_size(0.1)
            .with_unit(" °")
            .with_value_to_string(Arc::new(|v| format!("{v:.1} °"))),

            crossfeed: FloatParam::new(
                "Crossfeed Amount",
                -6.0,
                FloatRange::Linear {
                    min: -18.0,
                    max: 0.0,
                },
            )
            .with_smoother(SmoothingStyle::Linear(30.0))
            .with_step_size(0.01)
            .with_unit(" dB")
            .with_value_to_string(Arc::new(|v| format!("{v:+.2} dB"))),

            center_focus: FloatParam::new(
                "Centre Focus",
                50.0,
                FloatRange::Linear {
                    min: 0.0,
                    max: 100.0,
                },
            )
            .with_smoother(SmoothingStyle::Linear(30.0))
            .with_step_size(0.5)
            .with_unit(" %")
            .with_value_to_string(Arc::new(|v| format!("{v:.0} %"))),

            room_on: BoolParam::new("Room Externalisation", true),

            room_level: FloatParam::new(
                "Room Level",
                -22.0,
                FloatRange::Linear {
                    min: -30.0,
                    max: -12.0,
                },
            )
            .with_smoother(SmoothingStyle::Linear(30.0))
            .with_step_size(0.1)
            .with_unit(" dB")
            .with_value_to_string(Arc::new(db_string)),

            air: FloatParam::new(
                "Air / HF Trim",
                0.0,
                FloatRange::Linear {
                    min: -3.0,
                    max: 3.0,
                },
            )
            .with_smoother(SmoothingStyle::Linear(30.0))
            .with_step_size(0.01)
            .with_unit(" dB")
            .with_value_to_string(Arc::new(db_string)),

            auto_gain: BoolParam::new("Auto-Gain Match", true),
            bypass: BoolParam::new("Bypass", false),
        }
    }
}

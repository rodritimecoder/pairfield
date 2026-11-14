//! Pairfield — binaural crossfeed and virtual-monitor stage for open-back headphones.

pub mod dsp;
mod editor;
mod params;

use dsp::{DspParams, PairfieldDsp};
use editor::PairfieldEditor;
use nice_plug::editor::dpi::{LogicalSize, Size};
use nice_plug::prelude::*;
use nice_plug_egui::{
    create_egui_editor, EguiEditor, EguiEditorState, EguiNiceSettings, RepaintNotifier,
};
use params::PairfieldParams;
use std::sync::atomic::Ordering;
use std::sync::Arc;

/// Coefficient refresh cadence. 32 samples keeps the shelf/Thiran coefficients from stepping
/// audibly during a fast knob drag while costing only three transcendental calls per chunk.
const CHUNK: usize = 32;

const EDITOR_WIDTH: f64 = 560.0;
const EDITOR_HEIGHT: f64 = 430.0;

pub struct Pairfield {
    params: Arc<PairfieldParams>,
    dsp: PairfieldDsp,
    egui_state: Arc<EguiEditorState>,
    /// Published for the editor's gain readout. Written on the audio thread, read on the GUI
    /// thread; no locks.
    gain_meter: Arc<AtomicF32>,
}

impl Default for Pairfield {
    fn default() -> Self {
        Self {
            params: Arc::new(PairfieldParams::default()),
            dsp: PairfieldDsp::new(48000.0),
            egui_state: EguiEditorState::from_size(
                Size::Logical(LogicalSize::new(EDITOR_WIDTH, EDITOR_HEIGHT)),
                1.0,
            ),
            gain_meter: Arc::new(AtomicF32::new(1.0)),
        }
    }
}

impl Plugin for Pairfield {
    const NAME: &'static str = "Pairfield";
    const VENDOR: &'static str = "Pairfield Audio";
    const URL: &'static str = "https://example.com/pairfield";
    const EMAIL: &'static str = "hello@example.com";
    const VERSION: &'static str = env!("CARGO_PKG_VERSION");

    const AUDIO_IO_LAYOUTS: &'static [AudioIOLayout] = &[AudioIOLayout {
        main_input_channels: NonZeroU32::new(2),
        main_output_channels: NonZeroU32::new(2),
        ..AudioIOLayout::const_default()
    }];

    type Editor = EguiEditor<PairfieldEditor>;
    type SysExMessage = ();
    type BackgroundTask = ();

    fn params(&self) -> Arc<dyn Params> {
        self.params.clone()
    }

    fn editor(&mut self, _async_executor: AsyncExecutor<Self>) -> Option<Self::Editor> {
        let settings = EguiNiceSettings::new().with_tile(Self::NAME);
        create_egui_editor(
            self.egui_state.clone(),
            RepaintNotifier::new(),
            settings,
            PairfieldEditor::new(self.params.clone(), self.gain_meter.clone()),
        )
    }

    fn activate(
        &mut self,
        _layout: &AudioIOLayout,
        buffer_config: &BufferConfig,
        _context: &mut impl ActivateContext<Self>,
    ) -> bool {
        // All allocation happens here, never in `process`.
        self.dsp.prepare(buffer_config.sample_rate as f64);
        true
    }

    fn reset(&mut self) {
        self.dsp.reset();
    }

    fn process(
        &mut self,
        buffer: &mut Buffer,
        _aux: &mut AuxiliaryBuffers,
        _context: &mut impl ProcessContext<Self>,
    ) -> ProcessStatus {
        let frames = buffer.samples();
        if frames == 0 || buffer.channels() < 2 {
            return ProcessStatus::Normal;
        }

        let (left, right) = {
            let channels = buffer.as_slice();
            let (head, tail) = channels.split_at_mut(1);
            (&mut *head[0], &mut *tail[0])
        };

        let mut angle_buf = [0.0f32; CHUNK];
        let mut amount_buf = [0.0f32; CHUNK];
        let mut focus_buf = [0.0f32; CHUNK];
        let mut room_buf = [0.0f32; CHUNK];
        let mut air_buf = [0.0f32; CHUNK];

        let bypass = self.params.bypass.value();
        let room_on = self.params.room_on.value();
        let auto_gain = self.params.auto_gain.value();

        let mut pos = 0;
        while pos < frames {
            let len = CHUNK.min(frames - pos);

            self.params
                .speaker_angle
                .smoothed
                .next_block(&mut angle_buf[..len], len);
            self.params
                .crossfeed
                .smoothed
                .next_block(&mut amount_buf[..len], len);
            self.params
                .center_focus
                .smoothed
                .next_block(&mut focus_buf[..len], len);
            self.params
                .room_level
                .smoothed
                .next_block(&mut room_buf[..len], len);
            self.params.air.smoothed.next_block(&mut air_buf[..len], len);

            // Coefficients are refreshed from the value at the start of the chunk; the gains
            // inside the chunk therefore step by at most 1/75 of a full knob sweep.
            let p = DspParams {
                speaker_angle_deg: angle_buf[0],
                crossfeed_db: amount_buf[0],
                center_focus: focus_buf[0] / 100.0,
                room_active: room_on,
                room_db: room_buf[0],
                air_db: air_buf[0],
                auto_gain,
                bypass,
            };
            self.dsp.set_params(&p);
            self.dsp.process_block(&mut left[pos..pos + len], &mut right[pos..pos + len], &p);

            pos += len;
        }

        self.gain_meter
            .store(self.dsp.auto_gain().max(1e-6), Ordering::Relaxed);

        ProcessStatus::Normal
    }
}

impl Vst3Plugin for Pairfield {
    const VST3_CLASS_ID: [u8; 16] = *b"PairfieldVST3_01";
    const VST3_SUBCATEGORIES: &'static [Vst3SubCategory] =
        &[Vst3SubCategory::Fx, Vst3SubCategory::Tools];
}

impl ClapPlugin for Pairfield {
    const CLAP_ID: &'static str = "com.pairfield.binaural";
    const CLAP_DESCRIPTION: Option<&'static str> =
        Some("Binaural crossfeed and virtual-monitor stage for open-back headphones");
    const CLAP_MANUAL_URL: Option<&'static str> = Some(Self::URL);
    const CLAP_SUPPORT_URL: Option<&'static str> = None;
    const CLAP_FEATURES: &'static [ClapFeature] = &[
        ClapFeature::AudioEffect,
        ClapFeature::Stereo,
        ClapFeature::Utility,
    ];
}

nice_export_vst3!(Pairfield);
nice_export_clap!(Pairfield);

//! Minimal dark editor: a top-down head with draggable virtual speakers, the live transfer
//! functions, and the eight controls. Everything is painted by hand; no bitmaps.

use crate::dsp::{DspParams, PairfieldDsp};
use crate::params::PairfieldParams;
use egui::{Color32, FontId, Rect, Sense, Stroke, Vec2};
use nice_plug::context::gui::{GuiContext, ParamSetter};
use nice_plug::prelude::*;
use nice_plug_egui::widgets::ParamSlider;
use nice_plug_egui::{Frame, NiceEguiApp};
use std::sync::Arc;

const BG: Color32 = Color32::from_rgb(18, 19, 22);
const PANEL: Color32 = Color32::from_rgb(26, 28, 33);
const GRID: Color32 = Color32::from_rgb(44, 47, 54);
const DIM: Color32 = Color32::from_rgb(120, 126, 138);
const MONO: Color32 = Color32::from_rgb(236, 238, 242);
const CROSS: Color32 = Color32::from_rgb(72, 196, 178);
const ACCENT: Color32 = Color32::from_rgb(226, 168, 88);

pub struct PairfieldEditor {
    params: Arc<PairfieldParams>,
    gain: Arc<AtomicF32>,
    gui_ctx: Option<GuiContext>,
}

impl PairfieldEditor {
    pub fn new(params: Arc<PairfieldParams>, gain: Arc<AtomicF32>) -> Self {
        Self {
            params,
            gain,
            gui_ctx: None,
        }
    }

    fn snapshot(&self) -> DspParams {
        DspParams {
            speaker_angle_deg: self.params.speaker_angle.value(),
            crossfeed_db: self.params.crossfeed.value(),
            center_focus: self.params.center_focus.value() / 100.0,
            room_active: self.params.room_on.value(),
            room_db: self.params.room_level.value(),
            air_db: self.params.air.value(),
            auto_gain: self.params.auto_gain.value(),
            bypass: self.params.bypass.value(),
        }
    }

    fn draw(&mut self, ui: &mut egui::Ui) {
        let Some(gui_ctx) = self.gui_ctx.clone() else {
            ui.label("no gui context");
            return;
        };
        let setter = gui_ctx.param_setter();

        ui.painter().rect_filled(ui.max_rect(), 0.0, BG);

        // ---- header
        ui.horizontal(|ui| {
            ui.label(
                egui::RichText::new("PAIRFIELD")
                    .font(FontId::monospace(15.0))
                    .color(MONO)
                    .strong(),
            );
            ui.label(
                egui::RichText::new("binaural monitor · open-back")
                    .font(FontId::monospace(11.0))
                    .color(DIM),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let g = 20.0 * self.gain.load(std::sync::atomic::Ordering::Relaxed).log10();
                ui.label(
                    egui::RichText::new(format!("auto {g:+.1} dB"))
                        .font(FontId::monospace(11.0))
                        .color(if self.params.auto_gain.value() { CROSS } else { DIM }),
                );
            });
        });
        ui.add_space(4.0);

        // ---- head + curve side by side
        ui.horizontal(|ui| {
            let w = (ui.available_width() - 12.0) * 0.5;
            self.draw_head(ui, &setter, w, 196.0);
            self.draw_curve(ui, w, 196.0);
        });

        ui.add_space(6.0);

        // ---- controls
        let row = |ui: &mut egui::Ui, label: &str, param: &FloatParam, setter: &ParamSetter| {
            ui.horizontal(|ui| {
                ui.add_sized([118.0, 16.0], egui::Label::new(
                    egui::RichText::new(label).font(FontId::monospace(10.5)).color(DIM),
                ));
                ui.add_sized(
                    [ui.available_width() - 86.0, 16.0],
                    ParamSlider::for_param(param, setter),
                );
            });
        };

        row(ui, "speaker angle", &self.params.speaker_angle, &setter);
        row(ui, "crossfeed", &self.params.crossfeed, &setter);
        row(ui, "centre focus", &self.params.center_focus, &setter);
        row(ui, "air trim", &self.params.air, &setter);

        ui.add_space(2.0);
        ui.horizontal(|ui| {
            ui.add_sized([118.0, 16.0], egui::Label::new(
                egui::RichText::new("room").font(FontId::monospace(10.5)).color(DIM),
            ));
            let mut on = self.params.room_on.value();
            if ui
                .add(egui::Checkbox::without_text(&mut on))
                .on_hover_text("early reflections")
                .changed()
            {
                setter.set_parameter(&self.params.room_on, on);
            }
            ui.add_sized(
                [ui.available_width() - 120.0, 16.0],
                ParamSlider::for_param(&self.params.room_level, &setter),
            );
        });

        ui.add_space(4.0);
        ui.horizontal(|ui| {
            let mut ag = self.params.auto_gain.value();
            if ui
                .toggle_value(&mut ag, egui::RichText::new("auto-gain").size(11.0))
                .changed()
            {
                setter.set_parameter(&self.params.auto_gain, ag);
            }
            let mut by = self.params.bypass.value();
            if ui
                .toggle_value(&mut by, egui::RichText::new("bypass").size(11.0))
                .on_hover_text("bit-exact, level matched by the auto-gain")
                .changed()
            {
                setter.set_parameter(&self.params.bypass, by);
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let p = self.snapshot();
                let itd = crate::dsp::woodworth_delay_s(p.speaker_angle_deg as f64) * 1e6;
                ui.label(
                    egui::RichText::new(format!("ITD {:.0} µs", itd))
                        .font(FontId::monospace(10.5))
                        .color(DIM),
                );
            });
        });

        // ---- live indicator sweep so the header dot breathes
        ui.ctx().request_repaint();
    }

    /// Top-down head: two virtual speakers at ±angle, drag anywhere to set the angle.
    fn draw_head(&mut self, ui: &mut egui::Ui, setter: &ParamSetter, w: f32, h: f32) {
        let (resp, painter) = ui.allocate_painter(Vec2::new(w, h), Sense::click_and_drag());
        let rect = resp.rect;
        painter.rect_filled(rect, 4.0, PANEL);

        let center = egui::pos2(rect.center().x, rect.center().y + 6.0);
        let head_r = h * 0.20;
        let ring_r = h * 0.38;

        // floor grid
        for i in 0..=6 {
            let t = i as f32 / 6.0;
            painter.line_segment(
                [
                    egui::pos2(rect.left() + 8.0, rect.top() + 8.0 + t * (h - 16.0)),
                    egui::pos2(rect.right() - 8.0, rect.top() + 8.0 + t * (h - 16.0)),
                ],
                Stroke::new(1.0, Color32::from_rgb(31, 33, 38)),
            );
        }

        painter.circle_stroke(center, ring_r, Stroke::new(1.0, GRID));
        painter.circle_filled(center, head_r, Color32::from_rgb(35, 38, 44));
        painter.circle_stroke(center, head_r, Stroke::new(1.0, Color32::from_rgb(70, 75, 85)));
        // nose
        painter.add(egui::Shape::convex_polygon(
            vec![
                center + Vec2::new(-head_r * 0.28, -head_r * 0.94),
                center + Vec2::new(head_r * 0.28, -head_r * 0.94),
                center + Vec2::new(0.0, -head_r * 1.34),
            ],
            Color32::from_rgb(35, 38, 44),
            Stroke::new(1.0, Color32::from_rgb(70, 75, 85)),
        ));

        let angle = self.params.speaker_angle.value();
        let rad = angle.to_radians();

        // drag to set the angle
        if (resp.clicked() || resp.dragged()) && resp.interact_pointer_pos().is_some() {
            let pos = resp.interact_pointer_pos().unwrap();
            let dx = pos.x - center.x;
            let dy = center.y - pos.y; // up is "front"
            if dx.abs() > 1.0 || dy.abs() > 1.0 {
                let a = dx.atan2(dy).to_degrees().abs();
                let a = a.clamp(crate::params::MIN_ANGLE, crate::params::MAX_ANGLE);
                setter.set_parameter(&self.params.speaker_angle, a);
            }
        }

        for side in [-1.0f32, 1.0] {
            let dir = Vec2::new(rad.sin() * side, -rad.cos());
            let pos = center + dir * ring_r;
            let perp = Vec2::new(-dir.y, dir.x);
            let half = 9.0;
            let box_pts = vec![
                pos + perp * half + dir * 5.0,
                pos - perp * half + dir * 5.0,
                pos - perp * half - dir * 5.0,
                pos + perp * half - dir * 5.0,
            ];
            painter.add(egui::Shape::convex_polygon(
                box_pts,
                Color32::from_rgb(44, 48, 56),
                Stroke::new(1.0, ACCENT),
            ));
            // line from speaker to the ear
            let ear = center + Vec2::new(head_r * side, 0.0);
            painter.line_segment([pos, ear], Stroke::new(1.0, Color32::from_rgb(52, 58, 66)));
        }

        // crossfeed arc
        let arc_r = ring_r * 0.62;
        painter.line_segment(
            [
                center + Vec2::new(-arc_r, head_r * 0.0),
                center + Vec2::new(arc_r, 0.0),
            ],
            Stroke::new(1.0, Color32::from_rgb(60, 130, 122)),
        );

        painter.text(
            egui::pos2(rect.left() + 10.0, rect.top() + 8.0),
            egui::Align2::LEFT_TOP,
            format!(
                "{:.1}°  ·  ITD {:.0} µs",
                angle,
                crate::dsp::woodworth_delay_s(angle as f64) * 1e6
            ),
            FontId::monospace(10.0),
            DIM,
        );
    }

    /// Mono response (what a centred source does) and the cross-path ILD curve.
    fn draw_curve(&mut self, ui: &mut egui::Ui, w: f32, h: f32) {
        let (resp, painter) = ui.allocate_painter(Vec2::new(w, h), Sense::hover());
        let rect = resp.rect;
        painter.rect_filled(rect, 4.0, PANEL);

        let plot = Rect::from_min_max(
            egui::pos2(rect.left() + 30.0, rect.top() + 14.0),
            egui::pos2(rect.right() - 8.0, rect.bottom() - 16.0),
        );
        let db_min = -14.0f32;
        let db_max = 6.0f32;
        let fx = |f: f64| {
            let t = (f.log10() - 20f64.log10()) / (20000f64.log10() - 20f64.log10());
            plot.left() + t as f32 * plot.width()
        };
        let fy = |db: f32| {
            let t = (db - db_min) / (db_max - db_min);
            plot.bottom() - t * plot.height()
        };

        for f in [20.0, 100.0, 1000.0, 10000.0] {
            let x = fx(f);
            painter.line_segment(
                [egui::pos2(x, plot.top()), egui::pos2(x, plot.bottom())],
                Stroke::new(1.0, GRID),
            );
            painter.text(
                egui::pos2(x + 2.0, plot.bottom() + 2.0),
                egui::Align2::LEFT_TOP,
                if f >= 1000.0 {
                    format!("{:.0}k", f / 1000.0)
                } else {
                    format!("{f:.0}")
                },
                FontId::monospace(9.0),
                DIM,
            );
        }
        for db in [-12.0f32, -6.0, 0.0, 6.0] {
            let y = fy(db);
            painter.line_segment(
                [egui::pos2(plot.left(), y), egui::pos2(plot.right(), y)],
                Stroke::new(if db == 0.0 { 1.0 } else { 0.7 }, if db == 0.0 { Color32::from_rgb(70, 75, 85) } else { GRID }),
            );
            painter.text(
                egui::pos2(plot.left() - 4.0, y),
                egui::Align2::RIGHT_CENTER,
                format!("{db:+.0}"),
                FontId::monospace(9.0),
                DIM,
            );
        }

        let p = self.snapshot();
        let fs = 48000.0;
        let n = 160;
        let mut mono = Vec::with_capacity(n);
        let mut cross = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let f = 10f64.powf(20f64.log10() + t * (20000f64.log10() - 20f64.log10()));
            mono.push(egui::pos2(
                fx(f),
                fy(PairfieldDsp::mono_response_db(&p, fs, f) as f32),
            ));
            cross.push(egui::pos2(
                fx(f),
                fy(PairfieldDsp::cross_response_db(&p, fs, f) as f32),
            ));
        }
        painter.add(egui::Shape::line(cross, Stroke::new(1.2, CROSS)));
        painter.add(egui::Shape::line(mono, Stroke::new(1.6, MONO)));

        painter.text(
            egui::pos2(rect.left() + 10.0, rect.top() + 6.0),
            egui::Align2::LEFT_TOP,
            "mono sum",
            FontId::monospace(9.5),
            MONO,
        );
        painter.text(
            egui::pos2(rect.left() + 78.0, rect.top() + 6.0),
            egui::Align2::LEFT_TOP,
            "cross path",
            FontId::monospace(9.5),
            CROSS,
        );
    }
}

impl NiceEguiApp for PairfieldEditor {
    fn build(
        &mut self,
        _egui_ctx: egui::Context,
        nice_gui_ctx: GuiContext,
        _frame: &mut Frame,
    ) -> Result<(), nice_plug_egui::baseview::HandlerError> {
        self.gui_ctx = Some(nice_gui_ctx);
        Ok(())
    }

    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut Frame) {
        self.draw(ui);
    }

    fn editor_closed(&mut self) {
        self.gui_ctx = None;
    }
}

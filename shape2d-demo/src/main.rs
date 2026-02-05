use eframe::egui;
use egui_plot::{Line, Plot, Points, Polygon as PlotPolygon};
use shape2d::kernel::line::{BasicKernelF32WithCustomEpsilon, CapStyleF32, KernelF32};
use shape2d::triangle_kernel::TriangleKernelF32 as TriangleKernel;
use shape2d::{clean, clip, offset_raw, triangulate};

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1200.0, 800.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Shape2D Demo",
        options,
        Box::new(|_cc| Ok(Box::new(Demo::new()))),
    )
}

#[cfg(target_arch = "wasm32")]
fn main() {
    use eframe::wasm_bindgen::JsCast as _;

    console_error_panic_hook::set_once();

    wasm_bindgen_futures::spawn_local(async {
        let document = web_sys::window()
            .expect("No window")
            .document()
            .expect("No document");

        let canvas = document
            .get_element_by_id("canvas")
            .expect("Failed to find canvas")
            .dyn_into::<web_sys::HtmlCanvasElement>()
            .expect("Canvas element is not an HtmlCanvasElement");

        let web_options = eframe::WebOptions::default();

        eframe::WebRunner::new()
            .start(
                canvas,
                web_options,
                Box::new(|_cc| Ok(Box::new(Demo::new()))),
            )
            .await
            .expect("Failed to start eframe");
    });
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum WindingRule {
    Positive,
    NonZero,
    EvenOdd,
    GreaterOrEqual2,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum CapStyle {
    Arc,
    Bevel,
    Miter,
}

impl WindingRule {
    fn as_str(&self) -> &str {
        match self {
            WindingRule::Positive => "Positive",
            WindingRule::NonZero => "Non-Zero",
            WindingRule::EvenOdd => "Even-Odd",
            WindingRule::GreaterOrEqual2 => "2 or more (Intersection)",
        }
    }

    fn as_function(&self) -> Box<dyn Fn(i32) -> bool> {
        match self {
            WindingRule::Positive => Box::new(|w| w > 0),
            WindingRule::NonZero => Box::new(|w| w != 0),
            WindingRule::EvenOdd => Box::new(|w| w % 2 != 0),
            WindingRule::GreaterOrEqual2 => Box::new(|w| w >= 2),
        }
    }

    fn all() -> [WindingRule; 4] {
        [
            WindingRule::Positive,
            WindingRule::NonZero,
            WindingRule::EvenOdd,
            WindingRule::GreaterOrEqual2,
        ]
    }
}

impl CapStyle {
    fn as_str(&self) -> &str {
        match self {
            CapStyle::Arc => "Arc",
            CapStyle::Bevel => "Bevel",
            CapStyle::Miter => "Miter",
        }
    }

    fn to_cap_style_f32(&self) -> CapStyleF32 {
        match self {
            CapStyle::Arc => CapStyleF32::Arc {
                tolerance: ARC_TOLERANCE,
            },
            CapStyle::Bevel => CapStyleF32::Bevel,
            CapStyle::Miter => CapStyleF32::Miter { limit: MITER_LIMIT },
        }
    }

    fn all() -> [CapStyle; 3] {
        [CapStyle::Arc, CapStyle::Bevel, CapStyle::Miter]
    }
}

struct Demo {
    input_vertices: Vec<[f32; 2]>,
    input_edges: Vec<(u32, u32)>,
    processing_results: ProcessingResults,

    // UI state
    show_input: bool,
    show_cleaned: bool,
    show_clipped: bool,
    show_triangulation: bool,
    show_raw_offset: bool,
    show_cleaned_offset: bool,
    show_clipped_offset: bool,
    winding_rule: WindingRule,
    epsilon: f32,
    offset_amount: f32,
    cap_style: CapStyle,

    // Drag state
    dragging_vertex: Option<usize>,
    pointer_near_vertex: bool,

    // Edit popup state
    show_edit_popup: bool,
    edit_geometry_text: String,
}

const ARROW_SIZE_PIXELS: f64 = 12.0; // Arrow size in screen pixels
const ARC_TOLERANCE: f32 = 0.001;
const MITER_LIMIT: f32 = 4.;

struct ProcessingResults {
    kernel: BasicKernelF32WithCustomEpsilon,
    triangle_kernel: TriangleKernel,
    input_edges: Vec<(u32, u32)>,
    cleaned_edges: Vec<(u32, u32)>,
    clipped_edges: Vec<(u32, u32)>,
    raw_offset_edges: Vec<(u32, u32)>,
    cleaned_offset_edges: Vec<(u32, u32)>,
    clipped_offset_edges: Vec<(u32, u32)>,
    triangles: Vec<[u32; 3]>,
    error: Option<String>,
}

impl ProcessingResults {
    pub fn process(
        input_vertices: &[[f32; 2]],
        input_edges: &[(u32, u32)],
        winding_rule: WindingRule,
        epsilon: f32,
        offset_amount: f32,
        cap_style: CapStyleF32,
    ) -> Self {
        let input_edges = input_edges.to_vec();

        // Create kernel with vertices and epsilon
        let mut kernel = BasicKernelF32WithCustomEpsilon {
            points: input_vertices.to_vec(),
            epsilon,
        };

        // Step 1: Clean the edges (remove intersections)
        let cleaned_edges = clean(&mut kernel, input_edges.iter().copied());

        // Step 2: Clip with selected winding rule
        let winding_fn = winding_rule.as_function();
        let mut error = None;
        let clipped_edges = match clip(&mut kernel, cleaned_edges.iter().copied(), winding_fn) {
            Ok(edges) => edges,
            Err(e) => {
                error = Some(format!("Clipping error: {:?}", e));
                eprintln!("Clipping error! {:?}", e);
                dbg!(&kernel.points);
                dbg!(&cleaned_edges);
                Vec::new()
            }
        };

        // Step 3: Offset
        let raw_offset_edges = match offset_raw(
            &mut kernel,
            clipped_edges.iter().copied(),
            offset_amount,
            cap_style,
        ) {
            Ok(edges) => edges,
            Err(e) => {
                error = Some(format!("Offset error: {:?}", e));
                eprintln!("Offset error! {:?}", e);
                dbg!(&kernel.points);
                dbg!(&clipped_edges);
                Vec::new()
            }
        };

        // Step 4: Clean again
        let cleaned_offset_edges = clean(&mut kernel, raw_offset_edges.iter().copied());

        // Step 5: Clip again
        let clipped_offset_edges =
            match clip(&mut kernel, cleaned_offset_edges.iter().copied(), |f| f > 0) {
                Ok(edges) => edges,
                Err(e) => {
                    error = Some(format!("Offset clipping error: {:?}", e));
                    eprintln!("Offset clipping error! {:?}", e);
                    dbg!(&kernel.points);
                    dbg!(&cleaned_offset_edges);
                    Vec::new()
                }
            };

        // Step 6: Triangulate
        let mut triangle_kernel = TriangleKernel::new();
        let triangles = triangulate(
            &kernel,
            &mut triangle_kernel,
            clipped_offset_edges.iter().copied(),
        )
        .unwrap_or_else(|e| {
            error = Some(format!("Triangulation error: {:?}", e));
            eprintln!("Triangulation error! {:?}", e);
            dbg!(&kernel.points);
            dbg!(&cleaned_edges);
            dbg!(&clipped_edges);
            Vec::new()
        });

        Self {
            kernel,
            triangle_kernel,
            input_edges,
            cleaned_edges,
            clipped_edges,
            raw_offset_edges,
            cleaned_offset_edges,
            clipped_offset_edges,
            triangles,
            error,
        }
    }

    fn edges_to_arrows(&self, edges: &[(u32, u32)]) -> Vec<([f64; 2], [f64; 2])> {
        edges
            .iter()
            .map(|&(tail, head)| {
                let tail_pos = self.kernel.pt(tail);
                let head_pos = self.kernel.pt(head);
                (
                    [tail_pos[0] as f64, tail_pos[1] as f64],
                    [head_pos[0] as f64, head_pos[1] as f64],
                )
            })
            .collect()
    }

    fn edges_to_vertices(&self, edges: &[(u32, u32)]) -> Vec<[f64; 2]> {
        edges
            .iter()
            .flat_map(|&(tail, head)| {
                let tail = self.kernel.pt(tail);
                let head = self.kernel.pt(head);
                [
                    [tail[0] as f64, tail[1] as f64],
                    [head[0] as f64, head[1] as f64],
                ]
            })
            .collect()
    }

    fn triangles_to_polygons(&self) -> Vec<Vec<[f64; 2]>> {
        const MIN_AREA_THRESHOLD: f32 = 0.001; // Minimum triangle area to render

        self.triangles
            .iter()
            .filter_map(|&[v0, v1, v2]| {
                let p0 = self.triangle_kernel.pt(v0);
                let p1 = self.triangle_kernel.pt(v1);
                let p2 = self.triangle_kernel.pt(v2);

                // Filter out small sliver triangles
                let edge1_x = p1[0] - p0[0];
                let edge1_y = p1[1] - p0[1];
                let edge2_x = p2[0] - p0[0];
                let edge2_y = p2[1] - p0[1];
                let cross = edge1_x * edge2_y - edge1_y * edge2_x;
                let area = cross.abs() * 0.5;

                if area < MIN_AREA_THRESHOLD {
                    return None; // Skip small sliver triangles
                }

                Some(vec![
                    [p0[0] as f64, p0[1] as f64],
                    [p1[0] as f64, p1[1] as f64],
                    [p2[0] as f64, p2[1] as f64],
                ])
            })
            .collect()
    }

    fn triangles_to_arrows(&self) -> Vec<([f64; 2], [f64; 2])> {
        self.triangles
            .iter()
            .map(|&[v0, v1, v2]| {
                let p0 = self.triangle_kernel.pt(v0);
                let p1 = self.triangle_kernel.pt(v1);
                let p2 = self.triangle_kernel.pt(v2);

                [
                    ([p0[0] as f64, p0[1] as f64], [p1[0] as f64, p1[1] as f64]),
                    ([p1[0] as f64, p1[1] as f64], [p2[0] as f64, p2[1] as f64]),
                    ([p2[0] as f64, p2[1] as f64], [p0[0] as f64, p0[1] as f64]),
                ]
            })
            .flatten()
            .collect()
    }

    fn triangles_to_vertices(&self) -> Vec<[f64; 2]> {
        self.triangles
            .iter()
            .flat_map(|&[v0, v1, v2]| {
                let p0 = self.triangle_kernel.pt(v0);
                let p1 = self.triangle_kernel.pt(v1);
                let p2 = self.triangle_kernel.pt(v2);
                [
                    [p0[0] as f64, p0[1] as f64],
                    [p1[0] as f64, p1[1] as f64],
                    [p2[0] as f64, p2[1] as f64],
                ]
            })
            .collect()
    }
}

impl Demo {
    fn new() -> Self {
        // Create three shapes: a star, rectangle, and triangle
        // Star (5-pointed star, indices 0-9)
        let input_vertices = vec![
            [1., 4.0],
            [1.2351141, 3.3236067],
            [1.9510565, 3.309017],
            [1.3804226, 2.8763933],
            [1.5877852, 2.190983],
            [1.0, 2.6],
            [0.41221482, 2.190983],
            [0.61957735, 2.8763933],
            [0.04894346, 3.309017],
            [0.76488584, 3.3236067],
            [-0.2, 1.7],
            [2.3, 1.7],
            [2.3, 4.2],
            [-0.2, 4.2],
            [1.3, 1.],
            [3., 1.],
            [3., 4.],
        ];

        let input_edges = vec![
            // Star edges (0-9)
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
            // Rectangle edges (10-13)
            (10, 11),
            (11, 12),
            (12, 13),
            (13, 10),
            // Triangle edges (14-16)
            (14, 15),
            (15, 16),
            (16, 14),
        ];

        let winding_rule = WindingRule::Positive;
        let epsilon = 1e-5;
        let offset_amount = 0.1;
        let cap_style = CapStyle::Arc;
        let processing_results = ProcessingResults::process(
            &input_vertices,
            &input_edges,
            winding_rule,
            epsilon,
            offset_amount,
            cap_style.to_cap_style_f32(),
        );

        let edit_geometry_text = serde_json::to_string_pretty(&serde_json::json!({
            "vertices": input_vertices,
            "edges": input_edges,
        }))
        .unwrap();

        Self {
            input_vertices,
            input_edges,
            processing_results,
            show_input: true,
            show_cleaned: false,
            show_clipped: true,
            show_triangulation: false,
            show_raw_offset: false,
            show_cleaned_offset: false,
            show_clipped_offset: true,
            winding_rule,
            epsilon,
            offset_amount,
            cap_style,
            dragging_vertex: None,
            pointer_near_vertex: false,
            show_edit_popup: false,
            edit_geometry_text,
        }
    }

    fn update_from_text(&mut self) -> Result<(), String> {
        let json: serde_json::Value = serde_json::from_str(&self.edit_geometry_text)
            .map_err(|e| format!("Error parsing JSON: {}", e))?;

        let vertices: Vec<[f32; 2]> = serde_json::from_value(
            json.get("vertices")
                .ok_or("Missing 'vertices' key")?
                .clone(),
        )
        .map_err(|e| format!("Error parsing vertices: {}", e))?;

        let edges: Vec<(u32, u32)> =
            serde_json::from_value(json.get("edges").ok_or("Missing 'edges' key")?.clone())
                .map_err(|e| format!("Error parsing edges: {}", e))?;

        // Validate edges reference valid vertices
        let max_vertex = vertices.len() as u32;
        for &(a, b) in &edges {
            if a >= max_vertex || b >= max_vertex {
                return Err(format!(
                    "Edge ({}, {}) references invalid vertex (max index: {})",
                    a,
                    b,
                    max_vertex - 1
                ));
            }
        }

        self.input_vertices = vertices;
        self.input_edges = edges;
        self.processing_results = ProcessingResults::process(
            &self.input_vertices,
            &self.input_edges,
            self.winding_rule,
            self.epsilon,
            self.offset_amount,
            self.cap_style.to_cap_style_f32(),
        );

        Ok(())
    }

    fn sync_text_from_data(&mut self) {
        self.edit_geometry_text = serde_json::to_string_pretty(&serde_json::json!({
            "vertices": self.input_vertices,
            "edges": self.input_edges,
        }))
        .unwrap();
    }

    // Helper function to draw solid arrows
    fn draw_arrow(
        plot_ui: &mut egui_plot::PlotUi,
        origin: [f64; 2],
        tip: [f64; 2],
        color: egui::Color32,
        line_width: f32,
    ) {
        // Calculate arrow size in plot coordinates based on screen pixels
        let plot_bounds = plot_ui.plot_bounds();
        let screen_rect = plot_ui.response().rect;

        // Calculate pixels per plot unit (average of x and y scales)
        let pixels_per_unit_x = screen_rect.width() as f64 / plot_bounds.width();
        let pixels_per_unit_y = screen_rect.height() as f64 / plot_bounds.height();
        let pixels_per_unit = (pixels_per_unit_x + pixels_per_unit_y) / 2.0;

        // Convert arrow size from pixels to plot coordinates
        let arrow_size = ARROW_SIZE_PIXELS / pixels_per_unit;

        // Calculate direction vector
        let dx = tip[0] - origin[0];
        let dy = tip[1] - origin[1];
        let length = (dx * dx + dy * dy).sqrt();

        if length < 1e-6 {
            return; // Skip zero-length arrows
        }

        // Normalize direction
        let ux = dx / length;
        let uy = dy / length;

        // Calculate arrow base (where the head starts)
        let arrow_base = [tip[0] - ux * arrow_size, tip[1] - uy * arrow_size];

        // Draw line from origin to arrow base
        plot_ui.line(
            Line::new("", vec![origin, arrow_base])
                .color(color)
                .width(line_width),
        );

        // Calculate perpendicular vector for arrow head
        let perp_x = -uy;
        let perp_y = ux;
        let head_width = arrow_size * 0.5;

        // Create triangular arrow head
        let arrow_head = vec![
            tip,
            [
                arrow_base[0] + perp_x * head_width,
                arrow_base[1] + perp_y * head_width,
            ],
            [
                arrow_base[0] - perp_x * head_width,
                arrow_base[1] - perp_y * head_width,
            ],
        ];

        plot_ui.polygon(
            PlotPolygon::new("", arrow_head)
                .fill_color(color)
                .stroke(egui::Stroke::NONE),
        );
    }
}

impl Demo {
    fn colors(&self, ctx: &egui::Context) -> ColorScheme {
        let is_light = ctx.style().visuals.dark_mode == false;

        if is_light {
            ColorScheme {
                input: egui::Color32::from_rgb(140, 140, 140), // Darker gray
                cleaned: egui::Color32::from_rgb(200, 100, 0), // Darker orange
                clipped: egui::Color32::from_rgb(0, 180, 0),   // Darker green
                triangulation: egui::Color32::from_rgb(100, 120, 180), // Darker blue
                triangulation_fill: egui::Color32::from_rgba_unmultiplied(100, 120, 180, 128),
                raw_offset: egui::Color32::from_rgb(180, 0, 180), // Darker magenta
                cleaned_offset: egui::Color32::from_rgb(0, 150, 180), // Darker cyan
                clipped_offset: egui::Color32::from_rgb(180, 0, 0), // Darker red
            }
        } else {
            ColorScheme {
                input: egui::Color32::from_rgb(200, 200, 200), // Light gray
                cleaned: egui::Color32::from_rgb(255, 165, 0), // Orange
                clipped: egui::Color32::from_rgb(0, 255, 0),   // Bright green
                triangulation: egui::Color32::from_rgb(100, 149, 237), // Cornflower blue
                triangulation_fill: egui::Color32::from_rgba_unmultiplied(100, 149, 237, 128),
                raw_offset: egui::Color32::from_rgb(255, 0, 255), // Bright magenta
                cleaned_offset: egui::Color32::from_rgb(0, 255, 255), // Bright cyan
                clipped_offset: egui::Color32::from_rgb(255, 100, 100), // Light red
            }
        }
    }
}

struct ColorScheme {
    input: egui::Color32,
    cleaned: egui::Color32,
    clipped: egui::Color32,
    triangulation: egui::Color32,
    triangulation_fill: egui::Color32,
    raw_offset: egui::Color32,
    cleaned_offset: egui::Color32,
    clipped_offset: egui::Color32,
}

impl eframe::App for Demo {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let colors = self.colors(ctx);

        // Sidebar
        egui::SidePanel::left("controls").default_width(400.).show(ctx, |ui| {
            ui.heading("About");
            ui.separator();
            ui.label(
                "Shape2D is a geometry library. This demo showcases three of its algorithms:
 * Cleaning: Snapping together almost-equal geometry, and adding missing intersection points
 * Clipping: Removing and reversing edges that aren't on the boundary
 * Offsetting: Moving edges along their normal direction, and adding caps to the junctions
 * Triangulation: Tesselating the interior with triangles
 
It runs the input data through each of these steps sequentially, and shows you the output.
",
            );
            ui.heading("Display Options");
            ui.separator();

            let mut changed = false;

            ui.checkbox(&mut self.show_input, "Show Input Geometry");
            ui.checkbox(&mut self.show_cleaned, "Show Cleaned Edges");
            ui.checkbox(&mut self.show_clipped, "Show Clipped Result");
            ui.checkbox(&mut self.show_raw_offset, "Show Raw Offset Edges");
            ui.checkbox(&mut self.show_cleaned_offset, "Show Raw Offset Edges (cleaned)");
            ui.checkbox(&mut self.show_clipped_offset, "Show Final Offset Edges (clipped)");
            ui.checkbox(&mut self.show_triangulation, "Show Triangulation");

            ui.separator();
            ui.heading("Epsilon");
            ui.label(concat!("This defines how close geometry needs to get before it is considered coincident. ",
            "Generally it is kept at a very small constant value, just large enough to encompass floating-point rounding errors. ",
            "It is exposed here to make it easier to explore edge-cases that generally only happen on a microscopic level."));
            let epsilon_response =
                ui.add(egui::Slider::new(&mut self.epsilon, 1e-5..=1.0).logarithmic(true));

            if epsilon_response.changed() {
                self.processing_results = ProcessingResults::process(
                    &self.input_vertices,
                    &self.input_edges,
                    self.winding_rule,
                    self.epsilon,
                    self.offset_amount,
                    self.cap_style.to_cap_style_f32(),
                );
            }

            ui.separator();
            ui.heading("Winding Rule");
            ui.label("This defines how overlapping shapes are treated by the clipping step.");

            let current_rule = self.winding_rule;
            egui::ComboBox::from_id_salt("winding_rule")
                .selected_text(current_rule.as_str())
                .show_ui(ui, |ui| {
                    for rule in WindingRule::all() {
                        if ui
                            .selectable_value(&mut self.winding_rule, rule, rule.as_str())
                            .clicked()
                        {
                            changed = true;
                        }
                    }
                });

            if changed {
                self.processing_results = ProcessingResults::process(
                    &self.input_vertices,
                    &self.input_edges,
                    self.winding_rule,
                    self.epsilon,
                    self.offset_amount,
                    self.cap_style.to_cap_style_f32(),
                );
            }

            ui.separator();
            ui.heading("Offset Amount");

            let offset_response = ui.add(
                egui::Slider::new(&mut self.offset_amount, -1.0..=1.0)
                    .step_by(0.05)
            );

            if offset_response.changed() {
                self.processing_results = ProcessingResults::process(
                    &self.input_vertices,
                    &self.input_edges,
                    self.winding_rule,
                    self.epsilon,
                    self.offset_amount,
                    self.cap_style.to_cap_style_f32(),
                );
            }

            ui.separator();
            ui.heading("Offset Cap Style");
            let current_cap = self.cap_style;
            let mut cap_changed = false;
            egui::ComboBox::from_id_salt("cap_style")
                .selected_text(current_cap.as_str())
                .show_ui(ui, |ui| {
                    for style in CapStyle::all() {
                        if ui
                            .selectable_value(&mut self.cap_style, style, style.as_str())
                            .clicked()
                        {
                            cap_changed = true;
                        }
                    }
                });

            if cap_changed {
                self.processing_results = ProcessingResults::process(
                    &self.input_vertices,
                    &self.input_edges,
                    self.winding_rule,
                    self.epsilon,
                    self.offset_amount,
                    self.cap_style.to_cap_style_f32(),
                );
            }

            ui.separator();
            if ui.button("Edit Geometry as JSON...").clicked() {
                self.sync_text_from_data();
                self.show_edit_popup = true;
            }

            ui.separator();
            ui.heading("Legend");
            ui.horizontal(|ui| {
                let mut color = colors.input.to_array();
                ui.color_edit_button_srgba_unmultiplied(&mut color);
                ui.label("Input");
            });
            ui.horizontal(|ui| {
                let mut color = colors.cleaned.to_array();
                ui.color_edit_button_srgba_unmultiplied(&mut color);
                ui.label("Cleaned");
            });
            ui.horizontal(|ui| {
                let mut color = colors.clipped.to_array();
                ui.color_edit_button_srgba_unmultiplied(&mut color);
                ui.label("Clipped");
            });
            ui.horizontal(|ui| {
                let mut color = colors.raw_offset.to_array();
                ui.color_edit_button_srgba_unmultiplied(&mut color);
                ui.label("Raw Offset");
            });
            ui.horizontal(|ui| {
                let mut color = colors.cleaned_offset.to_array();
                ui.color_edit_button_srgba_unmultiplied(&mut color);
                ui.label("Raw Offset (cleaned)");
            });
            ui.horizontal(|ui| {
                let mut color = colors.clipped_offset.to_array();
                ui.color_edit_button_srgba_unmultiplied(&mut color);
                ui.label("Final Offset (clipped)");
            });
            ui.horizontal(|ui| {
                let mut color = colors.triangulation.to_array();
                ui.color_edit_button_srgba_unmultiplied(&mut color);
                ui.label("Triangulation");
            });
            ui.separator();
            ui.heading("Errors");
            if let Some(e) = &self.processing_results.error {
                ui.colored_label(egui::Color32::RED, e);
            } else {
                ui.label("No errors");
            }
        });

        // Main plot area
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Shape2D Demo");
            ui.label("Click and drag input vertices");

            // Disable plot dragging/scrolling when near a vertex or actively dragging
            let allow_plot_interaction =
                self.dragging_vertex.is_none() && !self.pointer_near_vertex;

            let plot_response = Plot::new("geometry_plot")
                .view_aspect(1.0)
                .data_aspect(1.0)
                .allow_drag(allow_plot_interaction)
                .allow_scroll(allow_plot_interaction)
                .show(ui, |plot_ui| {
                    // Handle vertex dragging
                    let pointer_pos = plot_ui.pointer_coordinate();
                    let pointer_down = plot_ui.ctx().input(|i| i.pointer.primary_down());
                    let pointer_released = plot_ui.ctx().input(|i| i.pointer.primary_released());

                    // Check if pointer is near a vertex
                    let mut near_vertex = false;
                    if let Some(pos) = pointer_pos {
                        let drag_threshold = 0.15; // Distance threshold in plot coordinates
                        for &vertex in self.input_vertices.iter() {
                            let dx = vertex[0] as f64 - pos.x;
                            let dy = vertex[1] as f64 - pos.y;
                            let dist = (dx * dx + dy * dy).sqrt();
                            if dist < drag_threshold {
                                near_vertex = true;
                                break;
                            }
                        }
                    }

                    // Start dragging if clicking near a vertex
                    if pointer_down && self.dragging_vertex.is_none() {
                        if let Some(pos) = pointer_pos {
                            let drag_threshold = 0.15; // Distance threshold in plot coordinates
                            for (i, &vertex) in self.input_vertices.iter().enumerate() {
                                let dx = vertex[0] as f64 - pos.x;
                                let dy = vertex[1] as f64 - pos.y;
                                let dist = (dx * dx + dy * dy).sqrt();
                                if dist < drag_threshold {
                                    self.dragging_vertex = Some(i);
                                    break;
                                }
                            }
                        }
                    }

                    // Update vertex position while dragging
                    if let Some(vertex_idx) = self.dragging_vertex {
                        if pointer_down {
                            if let Some(pos) = pointer_pos {
                                self.input_vertices[vertex_idx] = [pos.x as f32, pos.y as f32];
                                // Reprocess geometry in real-time
                                self.processing_results = ProcessingResults::process(
                                    &self.input_vertices,
                                    &self.input_edges,
                                    self.winding_rule,
                                    self.epsilon,
                                    self.offset_amount,
                                    self.cap_style.to_cap_style_f32(),
                                );
                            }
                        } else if pointer_released {
                            // Stop dragging
                            self.dragging_vertex = None;
                        }
                    }

                    // Filled triangles
                    if self.show_clipped {
                        // Draw filled triangles
                        for triangle_points in self.processing_results.triangles_to_polygons() {
                            plot_ui.polygon(
                                PlotPolygon::new("", triangle_points)
                                    .fill_color(colors.triangulation_fill)
                                    .stroke(egui::Stroke::NONE),
                            );
                        }
                    }

                    // Layer 1: Triangulation (filled, drawn first so it's underneath)
                    if self.show_triangulation {
                        // Draw triangle edges as arrows
                        let triangle_arrows = self.processing_results.triangles_to_arrows();
                        for (origin, tip) in triangle_arrows {
                            Self::draw_arrow(plot_ui, origin, tip, colors.triangulation, 2.0);
                        }

                        // Draw triangle vertices
                        let triangle_verts = self.processing_results.triangles_to_vertices();
                        if !triangle_verts.is_empty() {
                            plot_ui.points(
                                Points::new("", triangle_verts)
                                    .color(colors.triangulation)
                                    .radius(3.0),
                            );
                        }
                    }

                    // Layer 2: Input geometry (light gray)
                    if self.show_input {
                        let input_arrows = self
                            .processing_results
                            .edges_to_arrows(&self.processing_results.input_edges);
                        for (origin, tip) in input_arrows {
                            Self::draw_arrow(plot_ui, origin, tip, colors.input, 2.5);
                        }

                        // Draw input vertices (with highlighting for draggable/dragged)
                        for (i, &vertex) in self.input_vertices.iter().enumerate() {
                            let pos = [vertex[0] as f64, vertex[1] as f64];
                            let is_dragging = self.dragging_vertex == Some(i);

                            let (color, radius) = if is_dragging {
                                // Highlight the vertex being dragged
                                (egui::Color32::from_rgb(255, 255, 0), 6.0)
                            } else {
                                // Normal input vertex
                                (colors.input, 5.0)
                            };

                            plot_ui.points(
                                Points::new("", vec![pos])
                                    .color(color)
                                    .radius(radius)
                                    .shape(egui_plot::MarkerShape::Circle)
                                    .filled(true),
                            );
                        }
                    }

                    // Layer 3: Cleaned edges (orange)
                    if self.show_cleaned {
                        let cleaned_arrows = self
                            .processing_results
                            .edges_to_arrows(&self.processing_results.cleaned_edges);
                        for (origin, tip) in cleaned_arrows {
                            Self::draw_arrow(plot_ui, origin, tip, colors.cleaned, 2.0);
                        }

                        // Draw cleaned vertices
                        let cleaned_verts = self
                            .processing_results
                            .edges_to_vertices(&self.processing_results.cleaned_edges);
                        if !cleaned_verts.is_empty() {
                            plot_ui.points(
                                Points::new("", cleaned_verts)
                                    .color(colors.cleaned)
                                    .radius(3.5),
                            );
                        }
                    }

                    // Layer 4: Clipped result (green, topmost)
                    if self.show_clipped {
                        let clipped_arrows = self
                            .processing_results
                            .edges_to_arrows(&self.processing_results.clipped_edges);
                        for (origin, tip) in clipped_arrows {
                            Self::draw_arrow(plot_ui, origin, tip, colors.clipped, 2.5);
                        }

                        // Draw clipped vertices
                        let clipped_verts = self
                            .processing_results
                            .edges_to_vertices(&self.processing_results.clipped_edges);
                        if !clipped_verts.is_empty() {
                            plot_ui.points(
                                Points::new("", clipped_verts)
                                    .color(colors.clipped)
                                    .radius(4.0),
                            );
                        }
                    }

                    // Layer 5: Raw offset edges (magenta)
                    if self.show_raw_offset {
                        let raw_offset_arrows = self
                            .processing_results
                            .edges_to_arrows(&self.processing_results.raw_offset_edges);
                        for (origin, tip) in raw_offset_arrows {
                            Self::draw_arrow(plot_ui, origin, tip, colors.raw_offset, 2.0);
                        }

                        // Draw raw offset vertices
                        let raw_offset_verts = self
                            .processing_results
                            .edges_to_vertices(&self.processing_results.raw_offset_edges);
                        if !raw_offset_verts.is_empty() {
                            plot_ui.points(
                                Points::new("", raw_offset_verts)
                                    .color(colors.raw_offset)
                                    .radius(3.5),
                            );
                        }
                    }

                    // Layer 6: Cleaned offset edges (cyan)
                    if self.show_cleaned_offset {
                        let cleaned_offset_arrows = self
                            .processing_results
                            .edges_to_arrows(&self.processing_results.cleaned_offset_edges);
                        for (origin, tip) in cleaned_offset_arrows {
                            Self::draw_arrow(plot_ui, origin, tip, colors.cleaned_offset, 2.0);
                        }

                        // Draw cleaned offset vertices
                        let cleaned_offset_verts = self
                            .processing_results
                            .edges_to_vertices(&self.processing_results.cleaned_offset_edges);
                        if !cleaned_offset_verts.is_empty() {
                            plot_ui.points(
                                Points::new("", cleaned_offset_verts)
                                    .color(colors.cleaned_offset)
                                    .radius(3.5),
                            );
                        }
                    }

                    // Layer 7: Clipped offset edges (red)
                    if self.show_clipped_offset {
                        let clipped_offset_arrows = self
                            .processing_results
                            .edges_to_arrows(&self.processing_results.clipped_offset_edges);
                        for (origin, tip) in clipped_offset_arrows {
                            Self::draw_arrow(plot_ui, origin, tip, colors.clipped_offset, 2.5);
                        }

                        // Draw clipped offset vertices
                        let clipped_offset_verts = self
                            .processing_results
                            .edges_to_vertices(&self.processing_results.clipped_offset_edges);
                        if !clipped_offset_verts.is_empty() {
                            plot_ui.points(
                                Points::new("", clipped_offset_verts)
                                    .color(colors.clipped_offset)
                                    .radius(4.0),
                            );
                        }
                    }

                    // Return whether pointer is near a vertex for next frame
                    near_vertex
                });

            // Update state for next frame
            self.pointer_near_vertex = plot_response.inner;
        });

        // Edit geometry popup
        if self.show_edit_popup {
            let mut is_open = true;
            egui::Window::new("Edit Geometry")
                .open(&mut is_open)
                .resizable(true)
                .default_size([500.0, 600.0])
                .show(ctx, |ui| {
                    ui.label("Edit geometry as JSON:");
                    ui.label("Format: {\"vertices\": [[x1, y1], ...], \"edges\": [[v1, v2], ...]}");
                    ui.separator();

                    egui::ScrollArea::vertical()
                        .max_height(600.0)
                        .show(ui, |ui| {
                            ui.add(
                                egui::TextEdit::multiline(&mut self.edit_geometry_text)
                                    .font(egui::TextStyle::Monospace)
                                    .desired_width(f32::INFINITY)
                                    .desired_rows(30)
                                    .clip_text(true),
                            );
                        });
                    ui.separator();
                    ui.horizontal(|ui| {
                        if ui.button("Apply").clicked() {
                            match self.update_from_text() {
                                Ok(()) => {
                                    self.show_edit_popup = false;
                                }
                                Err(e) => {
                                    // Show error in the processing results
                                    self.processing_results.error = Some(e);
                                }
                            }
                        }
                        if ui.button("Cancel").clicked() {
                            self.show_edit_popup = false;
                        }
                    });

                    if let Some(err) = &self.processing_results.error {
                        ui.separator();
                        ui.colored_label(egui::Color32::RED, err);
                    }
                });

            if !is_open {
                self.show_edit_popup = false;
            }
        }
    }
}

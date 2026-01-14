use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints, Points, Polygon as PlotPolygon};
use polygon::kernel::polyline::F32 as Kernel;
use polygon::triangle_kernel::F32TriangleKernel as TriangleKernel;
use polygon::{clean, clip, triangulate};

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1200.0, 800.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Shape2D Interactive Demo",
        options,
        Box::new(|_cc| Ok(Box::new(InteractiveDemo::new()))),
    )
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum WindingRule {
    Positive,
    Negative,
    NonZero,
    EvenOdd,
    GreaterOrEqual2,
}

impl WindingRule {
    fn as_str(&self) -> &str {
        match self {
            WindingRule::Positive => "Positive",
            WindingRule::Negative => "Negative",
            WindingRule::NonZero => "Non-Zero",
            WindingRule::EvenOdd => "Even-Odd",
            WindingRule::GreaterOrEqual2 => "Intersection (>=2)",
        }
    }

    fn as_function(&self) -> Box<dyn Fn(i32) -> bool> {
        match self {
            WindingRule::Positive => Box::new(|w| w > 0),
            WindingRule::Negative => Box::new(|w| w < 0),
            WindingRule::NonZero => Box::new(|w| w != 0),
            WindingRule::EvenOdd => Box::new(|w| w % 2 != 0),
            WindingRule::GreaterOrEqual2 => Box::new(|w| w >= 2),
        }
    }

    fn all() -> Vec<WindingRule> {
        vec![
            WindingRule::Positive,
            WindingRule::Negative,
            WindingRule::NonZero,
            WindingRule::EvenOdd,
            WindingRule::GreaterOrEqual2,
        ]
    }
}

struct InteractiveDemo {
    input_vertices: Vec<[f32; 2]>,
    input_edges: Vec<(u32, u32)>,
    processing_results: ProcessingResults,

    // UI state
    show_input: bool,
    show_cleaned: bool,
    show_clipped: bool,
    show_triangulation: bool,
    winding_rule: WindingRule,

    // Drag state
    dragging_vertex: Option<usize>,
    pointer_near_vertex: bool,
}

struct ProcessingResults {
    kernel: Kernel,
    triangle_kernel: TriangleKernel,
    input_edges: Vec<(u32, u32)>,
    cleaned_edges: Vec<(u32, u32)>,
    clipped_edges: Vec<(u32, u32)>,
    triangles: Vec<(u32, u32, u32)>,
    error: Option<String>,
}

impl ProcessingResults {
    pub fn process(
        input_vertices: &[[f32; 2]],
        input_edges: &[(u32, u32)],
        winding_rule: WindingRule,
    ) -> Self {
        let input_edges = input_edges.to_vec();

        // Create kernel with vertices
        let mut kernel = Kernel::new(input_vertices.to_vec());

        // Step 1: Clean the edges (remove intersections)
        let cleaned_edges = clean(&mut kernel, input_edges.iter().copied());
        let cleaned_edges = clean(&mut kernel, cleaned_edges.iter().copied()); // XXX
        let cleaned_edges = clean(&mut kernel, cleaned_edges.iter().copied()); // XXX

        dbg!(&kernel);
        dbg!(&input_edges);
        dbg!(&cleaned_edges);

        // Step 2: Clip with selected winding rule
        let winding_fn = winding_rule.as_function();
        let clipped_edges = clip(&mut kernel, cleaned_edges.iter().copied(), winding_fn);

        // Step 3: Triangulate
        let mut error = None;
        let mut triangle_kernel = TriangleKernel::new();
        let triangles = triangulate(&kernel, &mut triangle_kernel, clipped_edges.iter().copied())
            .unwrap_or_else(|e| {
                error = Some(format!("Triangulation error: {:?}", e));
                Vec::new()
            });

        println!("Triangle vertices: {:?}", triangle_kernel.vertices);
        println!("Triangle indices: {:?}", triangles);

        Self {
            kernel,
            triangle_kernel,
            input_edges,
            cleaned_edges,
            clipped_edges,
            triangles,
            error,
        }
    }

    fn edges_to_arrows(&self, edges: &[(u32, u32)]) -> Vec<([f64; 2], [f64; 2])> {
        edges
            .iter()
            .map(|&(tail, head)| {
                let tail_pos = self.kernel.v(tail);
                let head_pos = self.kernel.v(head);
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
                let tail = self.kernel.v(tail);
                let head = self.kernel.v(head);
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
            .filter_map(|&(v0, v1, v2)| {
                let p0 = self.triangle_kernel.v(v0);
                let p1 = self.triangle_kernel.v(v1);
                let p2 = self.triangle_kernel.v(v2);

                // Filter out small sliver triangles by calculating area
                // Area = 0.5 * |cross product of two edge vectors|
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
        const MIN_AREA_THRESHOLD: f32 = 0.001;

        self.triangles
            .iter()
            .filter_map(|&(v0, v1, v2)| {
                let p0 = self.triangle_kernel.v(v0);
                let p1 = self.triangle_kernel.v(v1);
                let p2 = self.triangle_kernel.v(v2);

                // Filter out small sliver triangles
                let edge1_x = p1[0] - p0[0];
                let edge1_y = p1[1] - p0[1];
                let edge2_x = p2[0] - p0[0];
                let edge2_y = p2[1] - p0[1];
                let cross = edge1_x * edge2_y - edge1_y * edge2_x;
                let area = cross.abs() * 0.5;

                if area < MIN_AREA_THRESHOLD {
                    return None;
                }

                Some([
                    ([p0[0] as f64, p0[1] as f64], [p1[0] as f64, p1[1] as f64]),
                    ([p1[0] as f64, p1[1] as f64], [p2[0] as f64, p2[1] as f64]),
                    ([p2[0] as f64, p2[1] as f64], [p0[0] as f64, p0[1] as f64]),
                ])
            })
            .flatten()
            .collect()
    }

    fn triangles_to_vertices(&self) -> Vec<[f64; 2]> {
        self.triangles
            .iter()
            .flat_map(|&(v0, v1, v2)| {
                let p0 = self.triangle_kernel.v(v0);
                let p1 = self.triangle_kernel.v(v1);
                let p2 = self.triangle_kernel.v(v2);
                [
                    [p0[0] as f64, p0[1] as f64],
                    [p1[0] as f64, p1[1] as f64],
                    [p2[0] as f64, p2[1] as f64],
                ]
            })
            .collect()
    }
}

impl InteractiveDemo {
    fn new() -> Self {
        // Create two overlapping squares as initial geometry
        let input_vertices = vec![
            // First square (0-3)
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            //[0.0, 2.0],
            // Second square (4-7), offset to create overlap
            [1.0, 1.0],
            [3.0, 1.0],
            [3.0, 3.0],
            //[1.0, 3.0],
        ];

        let input_edges = vec![
            // First square
            (0, 1),
            (1, 2),
            (2, 0),
            //(3, 0),
            // Second square
            (3, 4),
            (4, 5),
            (5, 3),
            //(7, 4),
        ];
        let winding_rule = WindingRule::NonZero;
        let processing_results =
            ProcessingResults::process(&input_vertices, &input_edges, winding_rule);

        Self {
            input_vertices,
            input_edges,
            processing_results,
            show_input: true,
            show_cleaned: true,
            show_clipped: true,
            show_triangulation: true,
            winding_rule,
            dragging_vertex: None,
            pointer_near_vertex: false,
        }
    }

    // Helper function to draw solid arrows
    fn draw_arrow(
        plot_ui: &mut egui_plot::PlotUi,
        origin: [f64; 2],
        tip: [f64; 2],
        color: egui::Color32,
        line_width: f32,
        arrow_size: f64,
    ) {
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
            Line::new(vec![origin, arrow_base])
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
            PlotPolygon::new(arrow_head)
                .fill_color(color)
                .stroke(egui::Stroke::NONE),
        );
    }
}

impl eframe::App for InteractiveDemo {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Sidebar
        egui::SidePanel::left("controls").show(ctx, |ui| {
            ui.heading("Display Options");
            ui.separator();

            let mut changed = false;

            ui.checkbox(&mut self.show_input, "Show Input Geometry");
            ui.checkbox(&mut self.show_cleaned, "Show Cleaned Edges");
            ui.checkbox(&mut self.show_clipped, "Show Clipped Result");
            ui.checkbox(&mut self.show_triangulation, "Show Triangulation");

            ui.separator();
            ui.heading("Winding Rule");

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
                );
            }

            ui.separator();
            ui.heading("Legend");
            ui.horizontal(|ui| {
                ui.color_edit_button_srgba_unmultiplied(&mut [200, 200, 200, 255]);
                ui.label("Input");
            });
            ui.horizontal(|ui| {
                ui.color_edit_button_srgba_unmultiplied(&mut [255, 165, 0, 255]);
                ui.label("Cleaned");
            });
            ui.horizontal(|ui| {
                ui.color_edit_button_srgba_unmultiplied(&mut [0, 255, 0, 255]);
                ui.label("Clipped");
            });
            ui.horizontal(|ui| {
                ui.color_edit_button_srgba_unmultiplied(&mut [100, 149, 237, 255]);
                ui.label("Triangulation");
            });
            if let Some(e) = &self.processing_results.error {
                ui.separator();
                ui.label(e);
            }
        });

        // Main plot area
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Shape2D Geometry Processing Demo");
            ui.label("Click and drag input vertices to see real-time updates");

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
                                );
                            }
                        } else if pointer_released {
                            // Stop dragging
                            self.dragging_vertex = None;
                        }
                    }

                    // Layer 1: Triangulation (filled, drawn first so it's underneath)
                    if self.show_triangulation {
                        // Draw filled triangles
                        for triangle_points in self.processing_results.triangles_to_polygons() {
                            plot_ui.polygon(
                                PlotPolygon::new(PlotPoints::new(triangle_points))
                                    .fill_color(egui::Color32::from_rgba_unmultiplied(
                                        100, 149, 237, 128,
                                    ))
                                    .stroke(egui::Stroke::NONE),
                            );
                        }

                        // Draw triangle edges as arrows
                        let triangle_arrows = self.processing_results.triangles_to_arrows();
                        for (origin, tip) in triangle_arrows {
                            Self::draw_arrow(
                                plot_ui,
                                origin,
                                tip,
                                egui::Color32::from_rgb(100, 149, 237),
                                2.0,
                                0.08,
                            );
                        }

                        // Draw triangle vertices
                        let triangle_verts = self.processing_results.triangles_to_vertices();
                        if !triangle_verts.is_empty() {
                            plot_ui.points(
                                Points::new(triangle_verts)
                                    .color(egui::Color32::from_rgb(100, 149, 237))
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
                            Self::draw_arrow(
                                plot_ui,
                                origin,
                                tip,
                                egui::Color32::from_rgb(200, 200, 200),
                                2.5,
                                0.08,
                            );
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
                                (egui::Color32::from_rgb(200, 200, 200), 5.0)
                            };

                            plot_ui.points(
                                Points::new(vec![pos])
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
                            Self::draw_arrow(
                                plot_ui,
                                origin,
                                tip,
                                egui::Color32::from_rgb(255, 165, 0),
                                2.0,
                                0.08,
                            );
                        }

                        // Draw cleaned vertices
                        let cleaned_verts = self
                            .processing_results
                            .edges_to_vertices(&self.processing_results.cleaned_edges);
                        if !cleaned_verts.is_empty() {
                            plot_ui.points(
                                Points::new(cleaned_verts)
                                    .color(egui::Color32::from_rgb(255, 165, 0))
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
                            Self::draw_arrow(
                                plot_ui,
                                origin,
                                tip,
                                egui::Color32::from_rgb(0, 255, 0),
                                2.5,
                                0.08,
                            );
                        }

                        // Draw clipped vertices
                        let clipped_verts = self
                            .processing_results
                            .edges_to_vertices(&self.processing_results.clipped_edges);
                        if !clipped_verts.is_empty() {
                            plot_ui.points(
                                Points::new(clipped_verts)
                                    .color(egui::Color32::from_rgb(0, 255, 0))
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
    }
}

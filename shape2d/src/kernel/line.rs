use std::cmp::Ordering;

const EPSILON_MIN_F32: f32 = 1e-5;
const EPSILON_RATE_F32: f32 = 1e-5;

use crate::{
    kernel::{Edge, EdgeCoincidence, EdgeSide, Kernel, VertexEvent},
    rtree::Rect,
    sweep_line::{SweepLineChain, SweepLineEvent, SweepLineEventType, SweepLineSegment},
    triangle_kernel::{TriangleKernel, TriangleKernelF32},
};

pub trait KernelF32 {
    type Vertex: Copy + std::fmt::Debug + Ord;
    type Edge: Edge;

    fn epsilon(&self, fp_mag: f32) -> f32 {
        EPSILON_MIN_F32.max(EPSILON_RATE_F32 * fp_mag)
    }
    fn pt(&self, v: Self::Vertex) -> [f32; 2];
    fn vs(&self, e: Self::Edge) -> (Self::Vertex, Self::Vertex);

    fn merge_vertices(&mut self, a: Self::Vertex, b: Self::Vertex, pt: [f32; 2]) -> Self::Vertex;

    fn merge_edges(
        &mut self,
        a: Self::Edge,
        b: Self::Edge,
        coincidence: EdgeCoincidence,
    ) -> Self::Edge;

    fn new_edge_from_other(
        &mut self,
        old_edge: Self::Edge,
        start: Self::Vertex,
        end: Self::Vertex,
    ) -> Self::Edge;

    fn split_edge(&mut self, edge: Self::Edge, t: f32, pt: [f32; 2]) -> (Self::Edge, Self::Edge);

    fn new_vertex_from_other(&mut self, v: Self::Vertex, pt: [f32; 2]) -> Self::Vertex;

    fn new_cap_edge_from_vertex(
        &mut self,
        vertex: Self::Vertex,
        start: Self::Vertex,
        end: Self::Vertex,
        start_t: f32,
        end_t: f32,
    ) -> Self::Edge;
}

#[derive(Debug, Clone, Copy)]
pub enum CapStyleF32 {
    Arc { tolerance: f32 },
    Bevel,
    Miter { limit: f32 },
}

impl<K: KernelF32> Kernel for K
where
    (K::Vertex, K::Vertex): Edge,
{
    type Vertex = K::Vertex;
    type Edge = K::Edge;
    type Extents = ExtentsF32;
    type T = f32;
    type SweepLineEdgePortion = ();
    type SweepLineEventPoint = K::Vertex; // TODO rename this or see if it should be Point instead
    type TriangleKernel = TriangleKernelF32;
    type CapStyle = CapStyleF32;
    type OffsetAmount = f32;

    fn vertices_coincident(&self, a: Self::Vertex, b: Self::Vertex) -> bool {
        let pt_a = self.pt(a);
        let pt_b = self.pt(b);
        let fp_mag = fp_mag_pt_f32(pt_a).max(fp_mag_pt_f32(pt_b));
        points_coincident_f32(pt_a, pt_b, self.epsilon(fp_mag))
    }

    fn edges_coincident(&self, a: Self::Edge, b: Self::Edge) -> Option<EdgeCoincidence> {
        // Line segments are necessarily coincident if they share both vertices
        let (a_start, a_end) = self.vs(a);
        let (b_start, b_end) = self.vs(b);
        if a_start == b_start && a_end == b_end {
            Some(EdgeCoincidence::SameDirection)
        } else if a_start == b_end && a_end == b_start {
            Some(EdgeCoincidence::OppositeDirections)
        } else {
            None
        }
    }

    fn split(&self, vertex: Self::Vertex, edge: Self::Edge) -> Option<Self::T> {
        let vertex_pt = self.pt(vertex);
        let edge_vs = self.vs(edge);
        let edge_start_pt = self.pt(edge_vs.0);
        let edge_end_pt = self.pt(edge_vs.1);
        let fp_mag = fp_mag_pt_f32(vertex_pt)
            .max(fp_mag_pt_f32(edge_start_pt))
            .max(fp_mag_pt_f32(edge_end_pt));

        t_on_segment_f32(vertex_pt, edge_start_pt, edge_end_pt, self.epsilon(fp_mag))
    }

    fn intersection(&self, a: Self::Edge, b: Self::Edge) -> Option<(Self::T, Self::T)> {
        let a_vs = self.vs(a);
        let b_vs = self.vs(b);
        if a_vs.0 == b_vs.0 || a_vs.0 == b_vs.1 || a_vs.1 == b_vs.0 || a_vs.1 == b_vs.1 {
            // Segments that share an endpoint don't intersect
            return None;
        }

        intersect_segments_f32(
            self.pt(a_vs.0),
            self.pt(a_vs.1),
            self.pt(b_vs.0),
            self.pt(b_vs.1),
        )
    }

    fn merge_vertices(&mut self, a: Self::Vertex, b: Self::Vertex) -> Self::Vertex {
        let pt_a = self.pt(a);
        let pt_b = self.pt(b);
        let pt = merge_points_f32(pt_a, pt_b);
        <Self as KernelF32>::merge_vertices(self, a, b, pt)
    }

    fn merge_edges(
        &mut self,
        a: Self::Edge,
        b: Self::Edge,
        coincidence: EdgeCoincidence,
    ) -> Self::Edge {
        <Self as KernelF32>::merge_edges(self, a, b, coincidence)
    }

    fn split_edge(&mut self, edge: Self::Edge, t: Self::T) -> (Self::Edge, Self::Edge) {
        let (start_v, end_v) = self.vs(edge);
        let start_pt = self.pt(start_v);
        let end_pt = self.pt(end_v);
        let pt = pt_on_segment_f32(start_pt, end_pt, t);
        <Self as KernelF32>::split_edge(self, edge, t, pt)
    }

    fn replace_vertex_in_edge(
        &mut self,
        edge: Self::Edge,
        old_v: Self::Vertex,
        new_v: Self::Vertex,
    ) -> Option<Self::Edge> {
        let (mut start, mut end) = self.vs(edge);

        if start == old_v {
            if end == new_v {
                return None; // Self-edge
            }
            start = new_v;
        }
        if end == old_v {
            if start == new_v {
                return None; // Self-edge
            }
            end = new_v;
        }
        Some(self.new_edge_from_other(edge, start, end))
    }

    fn extents(&self, edges: impl Iterator<Item = Self::Edge>) -> Self::Extents {
        extents_f32(
            edges.flat_map(|e| {
                let (start, end) = self.vs(e);
                [self.pt(start), self.pt(end)]
            }),
            |fp_mag| self.epsilon(fp_mag),
        )
    }

    fn edge_bbox(&self, edge: Self::Edge, extents: &Self::Extents) -> Rect {
        let (edge_start, edge_end) = self.vs(edge);
        let edge_start_pt = self.pt(edge_start);
        let edge_end_pt = self.pt(edge_end);
        let fp_mag = fp_mag_pt_f32(edge_start_pt).max(fp_mag_pt_f32(edge_end_pt));
        segment_bbox_f32(edge_start_pt, edge_end_pt, *extents, self.epsilon(fp_mag))
    }

    fn sin_cmp(&self, common: Self::Vertex, a: Self::Vertex, b: Self::Vertex) -> Ordering {
        sin_cmp_f32(self.pt(common), self.pt(a), self.pt(b))
    }

    fn vertices_for_edge(&self, edge: Self::Edge) -> Option<(Self::Vertex, Self::Vertex)> {
        Some(self.vs(edge))
    }

    fn sweep_line_events_for_edge(
        &self,
        edge: Self::Edge,
    ) -> impl Iterator<Item = SweepLineEvent<Self>> {
        // Edges always have exactly 2 events.
        // We will use the `segment` data to store whether this edge is bottom (going  right) or top (going left)
        // based on how its endpoints sort in sweep-line order.
        let (start, end) = self.vs(edge);
        let chain = match sweep_line_cmp_f32(self.pt(start), self.pt(end)) {
            Ordering::Less => SweepLineChain::Bottom,
            Ordering::Equal => {
                panic!("Encountered a reflex edge (which are invalid)");
            }
            Ordering::Greater => SweepLineChain::Top,
        };

        let segment = SweepLineSegment::new(edge, (), chain);
        [
            SweepLineEvent {
                event_type: SweepLineEventType::Start,
                segment,
            },
            SweepLineEvent {
                event_type: SweepLineEventType::End,
                segment,
            },
        ]
        .into_iter()
    }

    fn sweep_line_event_cmp(&self, a: &SweepLineEvent<Self>, b: &SweepLineEvent<Self>) -> Ordering {
        let a_pt = self.pt(sweep_line_select_vertex(
            a.event_type,
            a.segment.chain,
            self.vs(a.segment.edge),
        ));
        let b_pt = self.pt(sweep_line_select_vertex(
            b.event_type,
            b.segment.chain,
            self.vs(b.segment.edge),
        ));

        // Compare first by event point (sweep-line order)
        sweep_line_cmp_f32(a_pt, b_pt)
            // Then by event type (End before Start)
            .then_with(|| a.event_type.cmp(&b.event_type))
            // Then by incidence angle, bottom-to-top
            .then_with(|| {
                let shared_event_type = a.event_type;
                let shared_pt = a_pt;
                let a_other_pt = self.pt(sweep_line_select_vertex(
                    shared_event_type.other(),
                    a.segment.chain,
                    self.vs(a.segment.edge),
                ));
                let b_other_pt = self.pt(sweep_line_select_vertex(
                    shared_event_type.other(),
                    b.segment.chain,
                    self.vs(b.segment.edge),
                ));
                match shared_event_type {
                    SweepLineEventType::End => sin_cmp_f32(shared_pt, b_other_pt, a_other_pt),
                    SweepLineEventType::Start => sin_cmp_f32(shared_pt, a_other_pt, b_other_pt),
                }
            })
    }

    fn sweep_line_event_point(&self, event: &SweepLineEvent<Self>) -> Self::SweepLineEventPoint {
        sweep_line_select_vertex(
            event.event_type,
            event.segment.chain,
            self.vs(event.segment.edge),
        )
    }

    fn sweep_line_segment_cmp(
        &self,
        segment: &SweepLineSegment<Self>,
        event_point: Self::SweepLineEventPoint,
    ) -> Ordering {
        let (start, end) = self.vs(segment.edge);
        let (left_i, right_i) = match segment.chain {
            SweepLineChain::Bottom => (start, end),
            SweepLineChain::Top => (end, start),
        };
        let left_pt = self.pt(left_i);
        let right_pt = self.pt(right_i);
        let common_pt = self.pt(event_point);
        sin_cmp_f32(left_pt, right_pt, common_pt)
    }

    fn sweep_line_event_point_to_triangle_vertex(
        &self,
        triangle_kernel: &mut Self::TriangleKernel,
        event_point: Self::SweepLineEventPoint,
    ) -> <Self::TriangleKernel as TriangleKernel>::Vertex {
        triangle_kernel.push_vertex(self.pt(event_point))
    }

    fn sweep_line_edge_segment_to_triangle_vertices(
        &self,
        _triangle_kernel: &mut Self::TriangleKernel,
        _segment: &SweepLineSegment<Self>,
    ) -> impl Iterator<Item = <Self::TriangleKernel as TriangleKernel>::Vertex> {
        // For polylines, we don't discretize edges into intermediate vertices
        // This could be extended for curves to discretize them
        None.into_iter()
    }

    fn vertex_event_cmp(&self, a: &VertexEvent<Self>, b: &VertexEvent<Self>) -> Ordering {
        let a_pt = self.pt(select_vertex(a.event_type, self.vs(a.edge)));
        let b_pt = self.pt(select_vertex(b.event_type, self.vs(b.edge)));

        // Compare first by event point (arbitrary; we'll pick sweep-line order)
        sweep_line_cmp_f32(a_pt, b_pt)
            // Then by incidence angle, we will have +X direction be the first,
            // and go CCW from there
            .then_with(|| {
                let shared_pt = a_pt;
                let a_other_pt = self.pt(select_vertex(a.event_type.other(), self.vs(a.edge)));
                let b_other_pt = self.pt(select_vertex(b.event_type.other(), self.vs(b.edge)));

                let a = [a_other_pt[0] - shared_pt[0], a_other_pt[1] - shared_pt[1]];
                let b = [b_other_pt[0] - shared_pt[0], b_other_pt[1] - shared_pt[1]];

                // Check coordinates to see if they're in different quadrants
                quadrant_f32(a).cmp(&quadrant_f32(b)).then_with(||
                    // Points are within 90 degrees of each other,
                    // so we can use sin_cmp to compare them
                    sin_cmp_f32(shared_pt, a_other_pt, b_other_pt))
            })
            .then_with(|| {
                // Then by event type--incoming before outgoing:
                //  ----->  sorted second
                // X
                //  <-----  sorted first
                // This sorting nudges two equal but opposite edges
                // in such a way that the area encountered between them is slightly negative,
                // encouraging the topology-finding algorithm
                // to avoid joining these two edges to each other if possible.
                a.event_type.cmp(&b.event_type).reverse()
            })
    }

    fn offset_edge_loops(
        &mut self,
        edge_loops: &[(u32, Self::Edge)],
        offset: Self::OffsetAmount,
        cap_style: Self::CapStyle,
    ) -> Vec<Self::Edge> {
        let mut result_edges = vec![];
        let mut raw_edges = vec![];

        fn offset_segment(a_pt: [f32; 2], b_pt: [f32; 2], offset: f32) -> ([f32; 2], [f32; 2]) {
            let offset_x = b_pt[1] - a_pt[1]; // Y
            let offset_y = a_pt[0] - b_pt[0]; // -X
            let norm_factor = offset / offset_x.hypot(offset_y);
            let offset_x = offset_x * norm_factor;
            let offset_y = offset_y * norm_factor;

            let new_a = [a_pt[0] + offset_x, a_pt[1] + offset_y];
            let new_b = [b_pt[0] + offset_x, b_pt[1] + offset_y];
            (new_a, new_b)
        }

        for edge_loop in edge_loops.chunk_by(|(i1, _), (i2, _)| i1 == i2) {
            // Offset each edge
            raw_edges.clear();
            for i in 0..edge_loop.len() {
                let (_, edge) = edge_loop[i];
                let (edge_start_v, edge_end_v) = self.vs(edge);
                let edge_start_pt = self.pt(edge_start_v);
                let edge_end_pt = self.pt(edge_end_v);
                // OPTIMIZATION:
                // Check for edge annihilation, which occurs if the segment is offset
                // beyond the intersection point of its two corner bisectors.
                // This works out to an offset limit of L / (tan(a/2) + tan(b/2))
                // where a and b are the corner angles and L is the segment length.
                // This can be rewritten as L * (cos(a) + 1.) * (cos(b) + 1.) /
                // ((cos(a) + 1.) * sin(b) + (cos(b) + 1.) * sin(a))
                // which can then be expressed using dot & cross products.
                {
                    let (_, prev_edge) = edge_loop[(i + edge_loop.len() - 1) % edge_loop.len()];
                    let (prev_edge_start_v, _) = self.vs(prev_edge);
                    let (_, next_edge) = edge_loop[(i + 1) % edge_loop.len()];
                    let (_, next_edge_end_v) = self.vs(next_edge);
                    let pt1 = self.pt(prev_edge_start_v);
                    let pt2 = edge_start_pt;
                    let pt3 = edge_end_pt;
                    let pt4 = self.pt(next_edge_end_v);
                    let l12 = [pt2[0] - pt1[0], pt2[1] - pt1[1]];
                    let l23 = [pt3[0] - pt2[0], pt3[1] - pt2[1]];
                    let l34 = [pt4[0] - pt3[0], pt4[1] - pt3[1]];
                    let a_dot = l12[0] * l23[0] + l12[1] * l23[1];
                    let mut a_cross = l12[0] * l23[1] - l12[1] * l23[0];
                    let b_dot = l23[0] * l34[0] + l23[1] * l34[1];
                    let mut b_cross = l23[0] * l34[1] - l23[1] * l34[0];

                    // Negate the angles based on the offset direction
                    if offset > 0. {
                        a_cross = -a_cross;
                        b_cross = -b_cross;
                    }

                    // Both corner angles must be concave
                    // in the direction of the offset
                    if a_cross > 0. && b_cross > 0. {
                        let len12 = l12[0].hypot(l12[1]);
                        let len23 = l23[0].hypot(l23[1]);
                        let len34 = l34[0].hypot(l34[1]);

                        // Compute sine and cosine using dot and cross product.
                        // To avoid division, both the numerator and denominator
                        // will have a common factor of len12 * len23 ^ 2 * len34
                        let a_cos_plus_1 = a_dot + len12 * len23;
                        let a_sin = a_cross;
                        let b_cos_plus_1 = b_dot + len23 * len34;
                        let b_sin = b_cross;

                        // Rearrange the inequality to avoid division
                        // (numerator and denominator are both always positive)
                        if len23 * a_cos_plus_1 * b_cos_plus_1
                            < offset.abs() * (a_cos_plus_1 * b_sin + b_cos_plus_1 * a_sin)
                        {
                            // Annihilate
                            continue;
                        }
                    }
                }

                // Offset the edge and push it to raw_edges
                let (offset_edge_start_pt, offset_edge_end_pt) =
                    offset_segment(edge_start_pt, edge_end_pt, offset);
                let offset_edge_start_v =
                    self.new_vertex_from_other(edge_start_v, offset_edge_start_pt);
                let offset_edge_end_v = self.new_vertex_from_other(edge_end_v, offset_edge_end_pt);
                let offset_edge =
                    self.new_edge_from_other(edge, offset_edge_start_v, offset_edge_end_v);
                raw_edges.push((edge, offset_edge));
            }

            // OPTIMIZATION:
            // Check each pair of offset edges for intersection
            // to save the cleaning / clipping steps some work.
            for prev_i in 0..raw_edges.len() {
                let next_i = (prev_i + 1) % raw_edges.len();
                let (prev_edge, prev_offset_edge) = raw_edges[prev_i];
                let (next_edge, next_offset_edge) = raw_edges[next_i];
                let (prev_offset_edge_start_v, prev_offset_edge_end_v) = self.vs(prev_offset_edge);
                let prev_offset_edge_start_pt = self.pt(prev_offset_edge_start_v);
                let prev_offset_edge_end_pt = self.pt(prev_offset_edge_end_v);
                let (next_offset_edge_start_v, next_offset_edge_end_v) = self.vs(next_offset_edge);
                let next_offset_edge_start_pt = self.pt(next_offset_edge_start_v);
                let next_offset_edge_end_pt = self.pt(next_offset_edge_end_v);

                let Some((t1, t2)) = intersect_segments_f32(
                    prev_offset_edge_start_pt,
                    prev_offset_edge_end_pt,
                    next_offset_edge_start_pt,
                    next_offset_edge_end_pt,
                ) else {
                    continue;
                };
                // Segments do intersect

                let (prev_offset_edge_shortened, _) =
                    <Self as Kernel>::split_edge(self, prev_offset_edge, t1);
                let split_vertex_a = self
                    .vertices_for_edge(prev_offset_edge_shortened)
                    .unwrap()
                    .1;

                let (_, next_offset_edge_shortened) =
                    <Self as Kernel>::split_edge(self, next_offset_edge, t2);
                let split_vertex_b = self
                    .vertices_for_edge(next_offset_edge_shortened)
                    .unwrap()
                    .0;

                let intersection_v =
                    <Self as Kernel>::merge_vertices(self, split_vertex_a, split_vertex_b);

                let prev_offset_edge_shortened = self
                    .replace_vertex_in_edge(
                        prev_offset_edge_shortened,
                        split_vertex_a,
                        intersection_v,
                    )
                    .unwrap();
                let next_offset_edge_shortened = self
                    .replace_vertex_in_edge(
                        next_offset_edge_shortened,
                        split_vertex_b,
                        intersection_v,
                    )
                    .unwrap();

                raw_edges[prev_i] = (prev_edge, prev_offset_edge_shortened);
                raw_edges[next_i] = (next_edge, next_offset_edge_shortened);
            }

            // Restore topology by iterating over each corner
            // and adding a cap if necessary
            for prev_i in 0..raw_edges.len() {
                let next_i = (prev_i + 1) % raw_edges.len();
                let (prev_edge, prev_offset_edge) = raw_edges[prev_i];
                let (_next_edge, next_offset_edge) = raw_edges[next_i];

                result_edges.push(prev_offset_edge);

                let (prev_offset_edge_start_v, prev_offset_edge_end_v) = self.vs(prev_offset_edge);
                let (next_offset_edge_start_v, next_offset_edge_end_v) = self.vs(next_offset_edge);

                // If this corner is already connected, there is nothing to do
                if prev_offset_edge_end_v == next_offset_edge_start_v {
                    continue;
                }

                let prev_offset_edge_start_pt = self.pt(prev_offset_edge_start_v);
                let prev_offset_edge_end_pt = self.pt(prev_offset_edge_end_v);
                let next_offset_edge_start_pt = self.pt(next_offset_edge_start_v);
                let next_offset_edge_end_pt = self.pt(next_offset_edge_end_v);

                let (_, original_v) = self.vs(prev_edge);
                // Could also have done (original_v, _) = self.vs(next_edge);
                // These should be equal, unless an intermediate edge was annihilated,
                // in which case we are dealing with a concave corner
                // and expect this cap to be clipped away.

                let original_pt = self.pt(original_v);

                let fp_mag = fp_mag_pt_f32(prev_offset_edge_end_pt)
                    .max(fp_mag_pt_f32(next_offset_edge_start_pt))
                    .max(fp_mag_pt_f32(original_pt));

                // If the two vertices are coincident,
                // or the corner is concave,
                // cap the corner with a straight line
                // (simple & topology-preserving.)
                // Given the annihilation step, I believe this is sufficient to guarantee a correct result.
                if points_coincident_f32(
                    prev_offset_edge_end_pt,
                    next_offset_edge_start_pt,
                    self.epsilon(fp_mag),
                ) || matches!(
                    (
                        sin_cmp_f32(
                            original_pt,
                            prev_offset_edge_end_pt,
                            next_offset_edge_start_pt
                        ),
                        offset >= 0.
                    ),
                    (Ordering::Greater, true) | (Ordering::Less, false)
                ) {
                    result_edges.push(self.new_cap_edge_from_vertex(
                        original_v,
                        prev_offset_edge_end_v,
                        next_offset_edge_start_v,
                        0.,
                        1.,
                    ));

                    continue;
                }

                // Now we know we are dealing with a convex corner
                // that needs a cap.
                match cap_style {
                    CapStyleF32::Arc { tolerance } => {
                        // Calculate vectors from center to both points
                        let dx_a = prev_offset_edge_end_pt[0] - original_pt[0];
                        let dy_a = prev_offset_edge_end_pt[1] - original_pt[1];
                        let dx_b = next_offset_edge_start_pt[0] - original_pt[0];
                        let dy_b = next_offset_edge_start_pt[1] - original_pt[1];

                        // Calculate the angle between the two vectors
                        let cross = dx_a * dy_b - dy_a * dx_b; // 2D cross product (z-component)
                        let dot = dx_a * dx_b + dy_a * dy_b; // dot product
                        let delta_angle = cross.atan2(dot);

                        // Calculate maximum angle per segment based on tolerance
                        // The sagitta (deviation) of a chord from an arc is: s = r * (1 - cos(θ/2))
                        // Solving for θ when s = tolerance: θ = 2 * arccos(1 - tolerance/r)
                        let cos_arg = (1.0 - tolerance.max(self.epsilon(fp_mag)) / offset.abs())
                            .clamp(-1.0, 1.0);
                        let max_angle_per_segment = 2.0 * cos_arg.acos();
                        let num_segments =
                            (delta_angle.abs() / max_angle_per_segment).ceil().max(1.) as u32;

                        // Interpolate the arc
                        let mut prev_vertex = prev_offset_edge_end_v;
                        let mut prev_t = 0.;
                        let inv_num_segments = 1. / num_segments as f32;
                        for i in 1..num_segments {
                            let t = i as f32 * inv_num_segments;
                            let theta = t * delta_angle;
                            let c = theta.cos();
                            let s = theta.sin();
                            // Rotate [dx_a, dy_a] by theta
                            let dx = c * dx_a - s * dy_a;
                            let dy = s * dx_a + c * dy_a;
                            let v = self.new_vertex_from_other(
                                original_v,
                                [original_pt[0] + dx, original_pt[1] + dy],
                            );
                            result_edges.push(self.new_cap_edge_from_vertex(
                                original_v,
                                prev_vertex,
                                v,
                                prev_t,
                                t,
                            ));
                            prev_vertex = v;
                            prev_t = t;
                        }
                        result_edges.push(self.new_cap_edge_from_vertex(
                            original_v,
                            prev_vertex,
                            next_offset_edge_start_v,
                            prev_t,
                            1.,
                        ));
                    }
                    CapStyleF32::Bevel => {
                        result_edges.push(self.new_cap_edge_from_vertex(
                            original_v,
                            prev_offset_edge_end_v,
                            next_offset_edge_start_v,
                            0.,
                            1.,
                        ));
                    }
                    CapStyleF32::Miter { limit } => {
                        let limit = limit * offset.abs();

                        let (miter_t_a, miter_t_b) = intersect_lines_f32(
                            prev_offset_edge_start_pt,
                            prev_offset_edge_end_pt,
                            next_offset_edge_start_pt,
                            next_offset_edge_end_pt,
                        );
                        let miter_pt_a = pt_on_segment_f32(
                            prev_offset_edge_start_pt,
                            prev_offset_edge_end_pt,
                            miter_t_a,
                        );
                        let miter_pt_b = pt_on_segment_f32(
                            next_offset_edge_start_pt,
                            next_offset_edge_end_pt,
                            miter_t_b,
                        );
                        let miter_pt = merge_points_f32(miter_pt_a, miter_pt_b);
                        let d_a = [
                            miter_pt[0] - prev_offset_edge_end_pt[0],
                            miter_pt[1] - prev_offset_edge_end_pt[1],
                        ];
                        let d_a_sq = d_a[0] * d_a[0] + d_a[1] * d_a[1];
                        if d_a_sq < limit * limit {
                            // Emit a true miter
                            let miter_v = self.new_vertex_from_other(original_v, miter_pt);
                            result_edges.push(self.new_cap_edge_from_vertex(
                                original_v,
                                prev_offset_edge_end_v,
                                miter_v,
                                0.,
                                0.5,
                            ));
                            result_edges.push(self.new_cap_edge_from_vertex(
                                original_v,
                                miter_v,
                                next_offset_edge_start_v,
                                0.5,
                                1.,
                            ));
                        } else {
                            // Emit a clipped miter
                            // Extend edges `a` and `b` by `limit`

                            let a_vec = [
                                prev_offset_edge_end_pt[0] - prev_offset_edge_start_pt[0],
                                prev_offset_edge_end_pt[1] - prev_offset_edge_start_pt[1],
                            ];
                            let a_len = a_vec[0].hypot(a_vec[1]);
                            let a_factor = limit / a_len;
                            let extended_a_pt = [
                                prev_offset_edge_end_pt[0] + a_vec[0] * a_factor,
                                prev_offset_edge_end_pt[1] + a_vec[1] * a_factor,
                            ];
                            let extended_a_v =
                                self.new_vertex_from_other(original_v, extended_a_pt);

                            let b_vec = [
                                next_offset_edge_start_pt[0] - next_offset_edge_end_pt[0],
                                next_offset_edge_start_pt[1] - next_offset_edge_end_pt[1],
                            ];
                            let b_len = b_vec[0].hypot(b_vec[1]);
                            let b_factor = limit / b_len;
                            let extended_b_pt = [
                                next_offset_edge_start_pt[0] + b_vec[0] * b_factor,
                                next_offset_edge_start_pt[1] + b_vec[1] * b_factor,
                            ];
                            let extended_b_v =
                                self.new_vertex_from_other(original_v, extended_b_pt);

                            let clipped_segment_length = (extended_b_pt[0] - extended_a_pt[0])
                                .hypot(extended_b_pt[1] - extended_a_pt[1]);
                            let miter_fraction = limit / (2. * limit + clipped_segment_length);

                            result_edges.push(self.new_cap_edge_from_vertex(
                                original_v,
                                prev_offset_edge_end_v,
                                extended_a_v,
                                0.,
                                miter_fraction,
                            ));
                            result_edges.push(self.new_cap_edge_from_vertex(
                                original_v,
                                extended_a_v,
                                extended_b_v,
                                miter_fraction,
                                1. - miter_fraction,
                            ));
                            result_edges.push(self.new_cap_edge_from_vertex(
                                original_v,
                                extended_b_v,
                                next_offset_edge_start_v,
                                1. - miter_fraction,
                                1.,
                            ));
                        }
                    }
                }
            }
        }
        result_edges
    }
}

#[repr(transparent)]
#[derive(Clone, Copy, Debug)]
pub struct BitOrdPointF32(pub [f32; 2]);

impl PartialEq for BitOrdPointF32 {
    fn eq(&self, other: &Self) -> bool {
        self.0[0].to_bits() == other.0[0].to_bits() && self.0[1].to_bits() == other.0[1].to_bits()
    }
}

impl Eq for BitOrdPointF32 {}

impl PartialOrd for BitOrdPointF32 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for BitOrdPointF32 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0[0]
            .to_bits()
            .cmp(&other.0[0].to_bits())
            .then_with(|| self.0[1].to_bits().cmp(&other.0[1].to_bits()))
    }
}

impl BitOrdPointF32 {
    const MIN: Self = BitOrdPointF32([f32::from_bits(u32::MIN), f32::from_bits(u32::MIN)]);
    const MAX: Self = BitOrdPointF32([f32::from_bits(u32::MAX), f32::from_bits(u32::MAX)]);
}

impl From<[f32; 2]> for BitOrdPointF32 {
    fn from(value: [f32; 2]) -> Self {
        Self(value)
    }
}

impl From<BitOrdPointF32> for [f32; 2] {
    fn from(value: BitOrdPointF32) -> Self {
        value.0
    }
}

pub struct DirectKernelF32 {}

impl Edge for (BitOrdPointF32, BitOrdPointF32) {
    const MIN: Self = (BitOrdPointF32::MIN, BitOrdPointF32::MIN);
    const MAX: Self = (BitOrdPointF32::MAX, BitOrdPointF32::MAX);

    fn reversed(self) -> Self {
        (self.1, self.0)
    }
}

impl KernelF32 for DirectKernelF32 {
    type Vertex = BitOrdPointF32;
    type Edge = (Self::Vertex, Self::Vertex);

    fn pt(&self, v: Self::Vertex) -> [f32; 2] {
        v.into()
    }

    fn vs(&self, e: Self::Edge) -> (Self::Vertex, Self::Vertex) {
        e
    }

    fn merge_vertices(&mut self, _a: Self::Vertex, _b: Self::Vertex, pt: [f32; 2]) -> Self::Vertex {
        pt.into()
    }

    fn merge_edges(
        &mut self,
        _a: Self::Edge,
        _b: Self::Edge,
        _coincidence: EdgeCoincidence,
    ) -> Self::Edge {
        panic!("Not possible for edges to share endpoints without being equal");
    }

    fn new_edge_from_other(
        &mut self,
        _old_edge: Self::Edge,
        start: Self::Vertex,
        end: Self::Vertex,
    ) -> Self::Edge {
        (start, end)
    }

    fn split_edge(&mut self, edge: Self::Edge, _t: f32, pt: [f32; 2]) -> (Self::Edge, Self::Edge) {
        ((edge.0, pt.into()), (pt.into(), edge.1))
    }

    fn new_vertex_from_other(&mut self, _v: Self::Vertex, pt: [f32; 2]) -> Self::Vertex {
        pt.into()
    }

    fn new_cap_edge_from_vertex(
        &mut self,
        _vertex: Self::Vertex,
        start: Self::Vertex,
        end: Self::Vertex,
        _start_t: f32,
        _end_t: f32,
    ) -> Self::Edge {
        (start, end)
    }
}

pub struct BasicKernelF32 {
    pub points: Vec<[f32; 2]>,
}

impl BasicKernelF32 {
    fn insert_vertex(&mut self, pt: [f32; 2]) -> u32 {
        let i = self.points.len() as u32;
        self.points.push(pt);
        i
    }
}

impl KernelF32 for BasicKernelF32 {
    type Vertex = u32;
    type Edge = (Self::Vertex, Self::Vertex);

    fn pt(&self, v: Self::Vertex) -> [f32; 2] {
        self.points[v as usize]
    }

    fn vs(&self, e: Self::Edge) -> (Self::Vertex, Self::Vertex) {
        e
    }

    fn merge_vertices(&mut self, _a: Self::Vertex, _b: Self::Vertex, pt: [f32; 2]) -> Self::Vertex {
        self.insert_vertex(pt)
    }

    fn merge_edges(
        &mut self,
        _a: Self::Edge,
        _b: Self::Edge,
        _coincidence: EdgeCoincidence,
    ) -> Self::Edge {
        panic!("Not possible for edges to share endpoints without being equal");
    }

    fn new_edge_from_other(
        &mut self,
        _old_edge: Self::Edge,
        start: Self::Vertex,
        end: Self::Vertex,
    ) -> Self::Edge {
        (start, end)
    }

    fn split_edge(&mut self, edge: Self::Edge, _t: f32, pt: [f32; 2]) -> (Self::Edge, Self::Edge) {
        let v = self.insert_vertex(pt);
        ((edge.0, v), (v, edge.1))
    }

    fn new_vertex_from_other(&mut self, _v: Self::Vertex, pt: [f32; 2]) -> Self::Vertex {
        self.insert_vertex(pt)
    }

    fn new_cap_edge_from_vertex(
        &mut self,
        _vertex: Self::Vertex,
        start: Self::Vertex,
        end: Self::Vertex,
        _start_t: f32,
        _end_t: f32,
    ) -> Self::Edge {
        (start, end)
    }
}

pub struct BasicKernelF32WithCustomEpsilon {
    pub points: Vec<[f32; 2]>,
    pub epsilon: f32,
}

impl BasicKernelF32WithCustomEpsilon {
    fn insert_vertex(&mut self, pt: [f32; 2]) -> u32 {
        let i = self.points.len() as u32;
        self.points.push(pt);
        i
    }
}

impl KernelF32 for BasicKernelF32WithCustomEpsilon {
    type Vertex = u32;
    type Edge = (Self::Vertex, Self::Vertex);

    fn epsilon(&self, fp_mag: f32) -> f32 {
        self.epsilon.max(EPSILON_RATE_F32 * fp_mag)
    }

    fn pt(&self, v: Self::Vertex) -> [f32; 2] {
        self.points[v as usize]
    }

    fn vs(&self, e: Self::Edge) -> (Self::Vertex, Self::Vertex) {
        e
    }

    fn merge_vertices(&mut self, _a: Self::Vertex, _b: Self::Vertex, pt: [f32; 2]) -> Self::Vertex {
        self.insert_vertex(pt)
    }

    fn merge_edges(
        &mut self,
        _a: Self::Edge,
        _b: Self::Edge,
        _coincidence: EdgeCoincidence,
    ) -> Self::Edge {
        panic!("Not possible for edges to share endpoints without being equal");
    }

    fn new_edge_from_other(
        &mut self,
        _old_edge: Self::Edge,
        start: Self::Vertex,
        end: Self::Vertex,
    ) -> Self::Edge {
        (start, end)
    }

    fn split_edge(&mut self, edge: Self::Edge, _t: f32, pt: [f32; 2]) -> (Self::Edge, Self::Edge) {
        let v = self.insert_vertex(pt);
        ((edge.0, v), (v, edge.1))
    }

    fn new_vertex_from_other(&mut self, _v: Self::Vertex, pt: [f32; 2]) -> Self::Vertex {
        self.insert_vertex(pt)
    }

    fn new_cap_edge_from_vertex(
        &mut self,
        _vertex: Self::Vertex,
        start: Self::Vertex,
        end: Self::Vertex,
        _start_t: f32,
        _end_t: f32,
    ) -> Self::Edge {
        (start, end)
    }
}

pub struct ColorKernelF32 {
    pub points: Vec<[f32; 2]>,
    pub colors: Vec<[f32; 4]>,
}

impl ColorKernelF32 {
    fn insert_vertex(&mut self, pt: [f32; 2], color: [f32; 4]) -> u32 {
        let i = self.points.len() as u32;
        self.points.push(pt);
        self.colors.push(color);
        i
    }

    fn color(&self, v: u32) -> [f32; 4] {
        self.colors[v as usize]
    }
}

fn average_colors(a: [f32; 4], b: [f32; 4]) -> [f32; 4] {
    [
        0.5 * (a[0] + b[0]),
        0.5 * (a[1] + b[1]),
        0.5 * (a[2] + b[2]),
        0.5 * (a[3] + b[3]),
    ]
}

fn weighted_average_colors(a: [f32; 4], b: [f32; 4], t: f32) -> [f32; 4] {
    [
        (1. - t) * a[0] + t * b[0],
        (1. - t) * a[1] + t * b[1],
        (1. - t) * a[2] + t * b[2],
        (1. - t) * a[3] + t * b[3],
    ]
}

impl KernelF32 for ColorKernelF32 {
    type Vertex = u32;
    type Edge = (Self::Vertex, Self::Vertex);

    fn pt(&self, v: Self::Vertex) -> [f32; 2] {
        self.points[v as usize]
    }

    fn vs(&self, e: Self::Edge) -> (Self::Vertex, Self::Vertex) {
        e
    }

    fn merge_vertices(&mut self, a: Self::Vertex, b: Self::Vertex, pt: [f32; 2]) -> Self::Vertex {
        self.insert_vertex(pt, average_colors(self.color(a), self.color(b)))
    }

    fn merge_edges(
        &mut self,
        _a: Self::Edge,
        _b: Self::Edge,
        _coincidence: EdgeCoincidence,
    ) -> Self::Edge {
        panic!("Not possible for edges to share endpoints without being equal");
    }

    fn new_edge_from_other(
        &mut self,
        _old_edge: Self::Edge,
        start: Self::Vertex,
        end: Self::Vertex,
    ) -> Self::Edge {
        (start, end)
    }

    fn split_edge(&mut self, edge: Self::Edge, t: f32, pt: [f32; 2]) -> (Self::Edge, Self::Edge) {
        let v = self.insert_vertex(
            pt,
            weighted_average_colors(self.color(edge.0), self.color(edge.1), t),
        );
        ((edge.0, v), (v, edge.1))
    }

    fn new_vertex_from_other(&mut self, v: Self::Vertex, pt: [f32; 2]) -> Self::Vertex {
        self.insert_vertex(pt, self.color(v))
    }

    fn new_cap_edge_from_vertex(
        &mut self,
        _vertex: Self::Vertex,
        start: Self::Vertex,
        end: Self::Vertex,
        _start_t: f32,
        _end_t: f32,
    ) -> Self::Edge {
        (start, end)
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct ExtentsF32 {
    scale: [f32; 2],
    offset: [f32; 2],
}

#[inline]
pub fn quadrant_f32([x, y]: [f32; 2]) -> u8 {
    if y > 0. && x <= 0. {
        1
    } else if y <= 0. && x < 0. {
        2
    } else if y < 0. && x >= 0. {
        3
    } else {
        0 // Also includes [0, 0]
    }
}

#[inline]
pub fn sweep_line_cmp_f32(a: [f32; 2], b: [f32; 2]) -> Ordering {
    a[0].partial_cmp(&b[0])
        .unwrap_or(Ordering::Equal)
        .then_with(|| a[1].partial_cmp(&b[1]).unwrap_or(Ordering::Equal))
}

#[inline]
pub fn sin_cmp_f32(common: [f32; 2], a: [f32; 2], b: [f32; 2]) -> Ordering {
    let ax = a[0] - common[0];
    let ay = a[1] - common[1];
    let bx = b[0] - common[0];
    let by = b[1] - common[1];

    // Check sign of cross product
    (ay * bx).partial_cmp(&(ax * by)).unwrap_or(Ordering::Equal)
}

#[inline]
fn points_coincident_f32(a: [f32; 2], b: [f32; 2], epsilon: f32) -> bool {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    dx * dx + dy * dy < epsilon * epsilon
}

#[inline]
fn pt_on_segment_f32(segment_start: [f32; 2], segment_end: [f32; 2], t: f32) -> [f32; 2] {
    let segment_dx = segment_end[0] - segment_start[0];
    let segment_dy = segment_end[1] - segment_start[1];

    [
        segment_start[0] + t * segment_dx,
        segment_start[1] + t * segment_dy,
    ]
}

#[inline]
fn t_on_segment_f32(
    p: [f32; 2],
    segment_start: [f32; 2],
    segment_end: [f32; 2],
    epsilon: f32,
) -> Option<f32> {
    let segment_dx = segment_end[0] - segment_start[0];
    let segment_dy = segment_end[1] - segment_start[1];
    let segment_len_sq = segment_dx * segment_dx + segment_dy * segment_dy;

    // Project point onto line
    // The dot product is the parallel distance along the line times the segment length
    let to_p_x = p[0] - segment_start[0];
    let to_p_y = p[1] - segment_start[1];
    let dot = segment_dx * to_p_x + segment_dy * to_p_y;

    // Check if projection lies outside the segment
    if dot <= 0. || dot >= segment_len_sq {
        return None;
    }

    // Get orthogonal distance to line using cross product
    let cross = segment_dy * to_p_x - segment_dx * to_p_y;
    // The cross product is the orthogonal distance to the line times the segment length
    // Since we only care about absolute distance,
    // we can do our comparison with the square of this value
    // and avoid a square root.

    if cross * cross >= epsilon * epsilon * segment_len_sq {
        return None;
    }

    let t = dot / segment_len_sq;
    Some(t)
}

#[inline]
fn intersect_segments_f32(
    a_start: [f32; 2],
    a_end: [f32; 2],
    b_start: [f32; 2],
    b_end: [f32; 2],
) -> Option<(f32, f32)> {
    // Perform separating axis test
    // Note: We will also return "no intersection" for any colinear segments.
    if matches!(
        (
            sin_cmp_f32(a_start, b_start, a_end),
            sin_cmp_f32(a_start, b_end, a_end)
        ),
        (Ordering::Equal, _)
            | (_, Ordering::Equal)
            | (Ordering::Greater, Ordering::Greater)
            | (Ordering::Less, Ordering::Less)
    ) {
        return None;
    }
    if matches!(
        (
            sin_cmp_f32(b_start, a_start, b_end),
            sin_cmp_f32(b_start, a_end, b_end)
        ),
        (Ordering::Equal, _)
            | (_, Ordering::Equal)
            | (Ordering::Greater, Ordering::Greater)
            | (Ordering::Less, Ordering::Less)
    ) {
        return None;
    }

    // We know the segments intersect, so fall back to line-line intersection algorithm
    Some(intersect_lines_f32(a_start, a_end, b_start, b_end))
}

#[inline]
fn merge_points_f32(a: [f32; 2], b: [f32; 2]) -> [f32; 2] {
    [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])]
}

fn extents_f32(
    mut points: impl Iterator<Item = [f32; 2]>,
    epsilon_for_mag: impl Fn(f32) -> f32,
) -> ExtentsF32 {
    let Some(first) = points.next() else {
        return Default::default();
    };
    let (min, max) = points.fold((first, first), |(min, max), point| {
        (
            [min[0].min(point[0]), min[1].min(point[1])],
            [max[0].max(point[0]), max[1].max(point[1])],
        )
    });

    let fp_mag = fp_mag_pt_f32(min).max(fp_mag_pt_f32(max));
    let epsilon = epsilon_for_mag(fp_mag);

    let scale = [
        u16::MAX as f32 / (max[0] - min[0] + 2. * epsilon),
        u16::MAX as f32 / (max[1] - min[1] + 2. * epsilon),
    ];
    let offset = [
        (-min[0] - epsilon) * scale[0],
        (-min[1] - epsilon) * scale[1],
    ];

    ExtentsF32 { scale, offset }
}

#[inline]
fn segment_bbox_f32(
    segment_start: [f32; 2],
    segment_end: [f32; 2],
    extents: ExtentsF32,
    epsilon: f32,
) -> Rect {
    let min = [
        segment_start[0].min(segment_end[0]) - epsilon,
        segment_start[1].min(segment_end[1]) - epsilon,
    ];
    let max = [
        segment_start[0].max(segment_end[0]) + epsilon,
        segment_start[1].max(segment_end[1]) + epsilon,
    ];

    let min = [
        f32_to_u16(extents.scale[0], extents.offset[0], min[0]),
        f32_to_u16(extents.scale[1], extents.offset[1], min[1]),
    ];
    let max = [
        f32_to_u16(extents.scale[0], extents.offset[0], max[0]),
        f32_to_u16(extents.scale[1], extents.offset[1], max[1]),
    ];

    Rect { min, max }
}

#[inline]
fn f32_to_u16(scale: f32, offset: f32, value: f32) -> u16 {
    (value * scale + offset).max(0.).min(u16::MAX as f32) as u16
}

#[inline]
fn sweep_line_select_vertex<T>(
    event_type: SweepLineEventType,
    edge_type: SweepLineChain,
    edge: (T, T),
) -> T {
    match (event_type, edge_type) {
        (SweepLineEventType::Start, SweepLineChain::Bottom)
        | (SweepLineEventType::End, SweepLineChain::Top) => edge.0,
        (SweepLineEventType::Start, SweepLineChain::Top)
        | (SweepLineEventType::End, SweepLineChain::Bottom) => edge.1,
    }
}

#[inline]
fn select_vertex<T>(event_type: EdgeSide, edge: (T, T)) -> T {
    match event_type {
        EdgeSide::Tail => edge.0,
        EdgeSide::Head => edge.1,
    }
}

#[inline]
fn intersect_lines_f32(
    a_start: [f32; 2],
    a_end: [f32; 2],
    b_start: [f32; 2],
    b_end: [f32; 2],
) -> (f32, f32) {
    let da = [a_end[0] - a_start[0], a_end[1] - a_start[1]];
    let db = [b_end[0] - b_start[0], b_end[1] - b_start[1]];

    // Inverse of the determinant
    let inv_det = 1. / (da[0] * db[1] - da[1] * db[0]);

    // Vector from a_start to b_start
    let diff = [b_start[0] - a_start[0], b_start[1] - a_start[1]];

    // Solve for t and u using Cramer's rule
    let t = inv_det * (diff[0] * db[1] - diff[1] * db[0]);
    let u = inv_det * (diff[0] * da[1] - diff[1] * da[0]);

    (t, u)
}

#[inline]
fn fp_mag_pt_f32(pt: [f32; 2]) -> f32 {
    pt[0].abs().max(pt[1].abs())
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use super::*;

    #[test]
    fn test_points_coincident_same_point() {
        let a = [0.5_f32, 0.5];
        let b = [0.5, 0.5];
        assert!(points_coincident_f32(a, b, EPSILON_MIN_F32));
    }

    #[test]
    fn test_points_coincident_within_epsilon() {
        let a = [0.5_f32, 0.5];
        let b = [0.500001, 0.500001];
        assert!(points_coincident_f32(a, b, EPSILON_MIN_F32));
    }

    #[test]
    fn test_points_coincident_outside_epsilon() {
        let a = [0.5_f32, 0.5];
        let b = [0.51, 0.5];
        assert!(!points_coincident_f32(a, b, EPSILON_MIN_F32));
    }

    #[test]
    fn test_points_coincident_far_apart() {
        let a = [0.0_f32, 0.0];
        let b = [1.0, 1.0];
        assert!(!points_coincident_f32(a, b, EPSILON_MIN_F32));
    }

    #[test]
    fn test_point_on_segment_midpoint() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let mid = [0.5, 0.0];
        let t = t_on_segment_f32(mid, start, end, EPSILON_MIN_F32).unwrap();
        let reconstructed = pt_on_segment_f32(start, end, t);
        assert!(points_coincident_f32(reconstructed, mid, EPSILON_MIN_F32));
    }

    #[test]
    fn test_point_on_segment_at_start() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.001, 0.0];
        let t = t_on_segment_f32(p, start, end, EPSILON_MIN_F32).unwrap();
        let reconstructed = pt_on_segment_f32(start, end, t);
        assert!(points_coincident_f32(reconstructed, p, EPSILON_MIN_F32));
    }

    #[test]
    fn test_point_on_segment_at_end() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.999, 0.0];
        let t = t_on_segment_f32(p, start, end, EPSILON_MIN_F32).unwrap();
        let reconstructed = pt_on_segment_f32(start, end, t);
        assert!(points_coincident_f32(reconstructed, p, EPSILON_MIN_F32));
    }

    #[test]
    fn test_point_on_segment_not_on_line() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.5, 0.1];
        assert!(t_on_segment_f32(p, start, end, EPSILON_MIN_F32).is_none());
    }

    #[test]
    fn test_point_on_segment_beyond_end() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [1.5, 0.0];
        assert!(t_on_segment_f32(p, start, end, EPSILON_MIN_F32).is_none());
    }

    #[test]
    fn test_point_on_segment_before_start() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [-0.5, 0.0];
        assert!(t_on_segment_f32(p, start, end, EPSILON_MIN_F32).is_none());
    }

    #[test]
    fn test_point_on_segment_diagonal() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 1.0];
        let mid = [0.5, 0.5];
        let t = t_on_segment_f32(mid, start, end, EPSILON_MIN_F32).unwrap();
        let reconstructed = pt_on_segment_f32(start, end, t);
        assert!(points_coincident_f32(reconstructed, mid, EPSILON_MIN_F32));
    }

    #[test]
    fn test_point_on_segment_within_epsilon_margin() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.5, 0.000001]; // Very close to the line
        let t = t_on_segment_f32(p, start, end, EPSILON_MIN_F32).unwrap();
        let reconstructed = pt_on_segment_f32(start, end, t);
        assert!(points_coincident_f32(reconstructed, p, EPSILON_MIN_F32));
    }

    #[test]
    fn test_intersect_segments_crossing_lines() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 1.0];
        let b_start = [0.0, 1.0];
        let b_end = [1.0, 0.0];

        let (t, u) = intersect_segments_f32(a_start, a_end, b_start, b_end).unwrap();

        // Verify the intersection point is [0.5, 0.5] from both segments
        let intersection_from_a = pt_on_segment_f32(a_start, a_end, t);
        let intersection_from_b = pt_on_segment_f32(b_start, b_end, u);
        assert!(points_coincident_f32(
            intersection_from_a,
            [0.5, 0.5],
            EPSILON_MIN_F32
        ));
        assert!(points_coincident_f32(
            intersection_from_b,
            [0.5, 0.5],
            EPSILON_MIN_F32
        ));
    }

    #[test]
    fn test_intersect_segments_t_junction() {
        let a_start = [0.0_f32, 0.5];
        let a_end = [1.0, 0.5];
        let b_start = [0.5, 0.0];
        let b_end = [0.5, 1.0];

        let (t, u) = intersect_segments_f32(a_start, a_end, b_start, b_end).unwrap();

        // Verify the intersection point is [0.5, 0.5] from both segments
        let intersection_from_a = pt_on_segment_f32(a_start, a_end, t);
        let intersection_from_b = pt_on_segment_f32(b_start, b_end, u);
        assert!(points_coincident_f32(
            intersection_from_a,
            [0.5, 0.5],
            EPSILON_MIN_F32
        ));
        assert!(points_coincident_f32(
            intersection_from_b,
            [0.5, 0.5],
            EPSILON_MIN_F32
        ));
    }

    #[test]
    fn test_intersect_segments_parallel() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 0.0];
        let b_start = [0.0, 1.0];
        let b_end = [1.0, 1.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end);
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_no_overlap() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.6, 0.4];
        let b_end = [1.0, 0.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end);
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_collinear() {
        // Two collinear segments that share an endpoint
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.5, 0.5];
        let b_end = [1.0, 1.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end);
        // Collinear segments return None (det is too small)
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_almost_parallel() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 0.0];
        let b_start = [0.0, 0.0];
        let b_end = [1.0, 0.000001]; // Almost parallel

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end);
        // Should be None because det is too small
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_close() {
        let a0 = [2.177941, 2.2992249];
        let a1 = [2.1857128, 2.2375135];
        let b0 = [2.1857483, 2.2371936];
        let b1 = [2.1862082, 2.2374542];

        let intersection = intersect_segments_f32(a0, a1, b0, b1);
        assert!(intersection.is_none());
    }

    #[test]
    fn test_merge_points_simple() {
        let a = [0.0_f32, 0.0];
        let b = [1.0, 1.0];
        let merged = merge_points_f32(a, b);
        assert!(points_coincident_f32(merged, [0.5, 0.5], EPSILON_MIN_F32));
    }

    #[test]
    fn test_merge_points_same_point() {
        let a = [0.5_f32, 0.5];
        let b = [0.5, 0.5];
        let merged = merge_points_f32(a, b);
        assert!(points_coincident_f32(merged, [0.5, 0.5], EPSILON_MIN_F32));
    }

    #[test]
    fn test_merge_points_negative_coords() {
        let a = [-1.0_f32, -1.0];
        let b = [1.0, 1.0];
        let merged = merge_points_f32(a, b);
        assert!(points_coincident_f32(merged, [0.0, 0.0], EPSILON_MIN_F32));
    }

    #[test]
    fn test_segment_bbox_horizontal() {
        let start = [0.2_f32, 0.5];
        let end = [0.8, 0.5];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), |_| EPSILON_MIN_F32);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON_MIN_F32);

        assert!(bbox.min[0] < bbox.max[0]);
        assert!(bbox.min[1] <= bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bbox_vertical() {
        let start = [0.5_f32, 0.2];
        let end = [0.5, 0.8];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), |_| EPSILON_MIN_F32);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON_MIN_F32);

        assert!(bbox.min[0] <= bbox.max[0]);
        assert!(bbox.min[1] < bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bbox_diagonal() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 1.0];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), |_| EPSILON_MIN_F32);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON_MIN_F32);

        assert!(bbox.min[0] < bbox.max[0]);
        assert!(bbox.min[1] < bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bbox_point() {
        let start = [0.5_f32, 0.5];
        let end = [0.5, 0.5];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), |_| EPSILON_MIN_F32);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON_MIN_F32);

        assert!(bbox.min[0] <= bbox.max[0]);
        assert!(bbox.min[1] <= bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bboxes_overlap() {
        // Segment 1 is vertical
        let start1 = [0.5_f32, 0.];
        let end1 = [0.5, 1.];

        // Segment 2 is horizontal
        let start2 = [0.5_f32, 0.5];
        let end2 = [1., 0.5];

        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), |_| EPSILON_MIN_F32);
        let bbox1 = segment_bbox_f32(start1, end1, extents, EPSILON_MIN_F32);
        let bbox2 = segment_bbox_f32(start2, end2, extents, EPSILON_MIN_F32);

        assert!(bbox1.overlaps(&bbox2));
    }

    #[test]
    fn test_segment_bboxes_dont_overlap() {
        let start1 = [0.1_f32, 0.1];
        let end1 = [0.2, 0.2];

        let start2 = [0.2_f32, 0.4];
        let end2 = [0.1, 0.5];

        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), |_| EPSILON_MIN_F32);
        let bbox1 = segment_bbox_f32(start1, end1, extents, EPSILON_MIN_F32);
        let bbox2 = segment_bbox_f32(start2, end2, extents, EPSILON_MIN_F32);

        assert!(!bbox1.overlaps(&bbox2));
    }

    #[test]
    fn test_sweep_line_cmp() {
        // X takes precedence over Y
        let a = [0.0_f32, 1.0];
        let b = [1.0, 0.0];
        assert_eq!(sweep_line_cmp_f32(a, b), Ordering::Less);

        // Y is compared when X is equal
        let c = [0.5_f32, 0.0];
        let d = [0.5, 1.0];
        assert_eq!(sweep_line_cmp_f32(c, d), Ordering::Less);
    }

    #[test]
    fn test_sin_cmp() {
        let origin = [0.0_f32, 0.0];
        let right = [1.0, 0.0];
        let up = [0.0, 1.0];
        let diagonal = [1.0, 1.0];
        let diagonal2 = [2.0, 2.0];

        // Counterclockwise ordering
        assert_eq!(sin_cmp_f32(origin, right, up), Ordering::Less);
        assert_eq!(sin_cmp_f32(origin, up, right), Ordering::Greater);

        // Collinear points
        assert_eq!(sin_cmp_f32(origin, diagonal, diagonal2), Ordering::Equal);
    }

    #[test]
    fn test_quadrant_positive_x() {
        assert_eq!(quadrant_f32([1.0, 0.0]), 0);
        assert_eq!(quadrant_f32([5.0, 0.0]), 0);
    }

    #[test]
    fn test_quadrant_positive_y_negative_x() {
        assert_eq!(quadrant_f32([-1.0, 1.0]), 1);
        assert_eq!(quadrant_f32([-5.0, 5.0]), 1);
        assert_eq!(quadrant_f32([0.0, 1.0]), 1);
    }

    #[test]
    fn test_quadrant_negative_x_negative_y() {
        assert_eq!(quadrant_f32([-1.0, -1.0]), 2);
        assert_eq!(quadrant_f32([-5.0, -5.0]), 2);
        assert_eq!(quadrant_f32([-1.0, 0.0]), 2);
    }

    #[test]
    fn test_quadrant_negative_y_positive_x() {
        assert_eq!(quadrant_f32([1.0, -1.0]), 3);
        assert_eq!(quadrant_f32([5.0, -5.0]), 3);
    }

    #[test]
    fn test_quadrant_origin() {
        assert_eq!(quadrant_f32([0.0, 0.0]), 0);
    }

    #[test]
    fn test_quadrant_positive_x_axis() {
        assert_eq!(quadrant_f32([1.0, 0.0]), 0);
    }

    #[test]
    fn test_quadrant_negative_x_axis() {
        assert_eq!(quadrant_f32([-1.0, 0.0]), 2);
    }

    #[test]
    fn test_quadrant_positive_y_axis() {
        assert_eq!(quadrant_f32([0.0, 1.0]), 1);
    }

    #[test]
    fn test_quadrant_negative_y_axis() {
        assert_eq!(quadrant_f32([0.0, -1.0]), 3);
    }

    #[test]
    fn test_intersect_lines_perpendicular() {
        let a1 = [0.0, 0.0];
        let a2 = [1.0, 0.0]; // Horizontal line through origin
        let b1 = [0.5, -1.0];
        let b2 = [0.5, 1.0]; // Vertical line through x=0.5

        let (t, u) = intersect_lines_f32(a1, a2, b1, b2);
        let intersection_from_a = pt_on_segment_f32(a1, a2, t);
        let intersection_from_b = pt_on_segment_f32(b1, b2, u);
        assert!(points_coincident_f32(
            intersection_from_a,
            [0.5, 0.0],
            EPSILON_MIN_F32
        ));
        assert!(points_coincident_f32(
            intersection_from_b,
            [0.5, 0.0],
            EPSILON_MIN_F32
        ));
    }

    #[test]
    fn test_intersect_lines_diagonal() {
        let a1 = [0.0, 0.0];
        let a2 = [1.0, 1.0]; // Line y = x
        let b1 = [0.0, 1.0];
        let b2 = [1.0, 0.0]; // Line y = -x + 1

        let (t, u) = intersect_lines_f32(a1, a2, b1, b2);
        let intersection_from_a = pt_on_segment_f32(a1, a2, t);
        let intersection_from_b = pt_on_segment_f32(b1, b2, u);
        assert!(points_coincident_f32(
            intersection_from_a,
            [0.5, 0.5],
            EPSILON_MIN_F32
        ));
        assert!(points_coincident_f32(
            intersection_from_b,
            [0.5, 0.5],
            EPSILON_MIN_F32
        ));
    }

    #[test]
    fn test_intersect_lines_at_origin() {
        let a1 = [-1.0, 0.0];
        let a2 = [1.0, 0.0]; // Horizontal line through origin
        let b1 = [0.0, -1.0];
        let b2 = [0.0, 1.0]; // Vertical line through origin

        let (t, u) = intersect_lines_f32(a1, a2, b1, b2);
        let intersection_from_a = pt_on_segment_f32(a1, a2, t);
        let intersection_from_b = pt_on_segment_f32(b1, b2, u);
        assert!(points_coincident_f32(
            intersection_from_a,
            [0.0, 0.0],
            EPSILON_MIN_F32
        ));
        assert!(points_coincident_f32(
            intersection_from_b,
            [0.0, 0.0],
            EPSILON_MIN_F32
        ));
    }

    #[test]
    fn test_intersect_lines_various_angles() {
        // Line 1: from (0,0) to (2,1) → y = x/2
        let a1 = [0.0, 0.0];
        let a2 = [2.0, 1.0];
        // Line 2: from (0,1) to (2,0) → y = 1 - x/2
        let b1 = [0.0, 1.0];
        let b2 = [2.0, 0.0];

        let (t, u) = intersect_lines_f32(a1, a2, b1, b2);
        // These lines should intersect at (1, 0.5)
        let intersection_from_a = pt_on_segment_f32(a1, a2, t);
        let intersection_from_b = pt_on_segment_f32(b1, b2, u);
        assert!(points_coincident_f32(
            intersection_from_a,
            [1.0, 0.5],
            EPSILON_MIN_F32
        ));
        assert!(points_coincident_f32(
            intersection_from_b,
            [1.0, 0.5],
            EPSILON_MIN_F32
        ));
    }

    #[test]
    fn test_intersect_lines_negative_coords() {
        let a1 = [-2.0, -2.0];
        let a2 = [2.0, 2.0]; // Line through origin, slope 1
        let b1 = [-2.0, 2.0];
        let b2 = [2.0, -2.0]; // Line through origin, slope -1

        let (t, u) = intersect_lines_f32(a1, a2, b1, b2);
        let intersection_from_a = pt_on_segment_f32(a1, a2, t);
        let intersection_from_b = pt_on_segment_f32(b1, b2, u);
        assert!(points_coincident_f32(
            intersection_from_a,
            [0.0, 0.0],
            EPSILON_MIN_F32
        ));
        assert!(points_coincident_f32(
            intersection_from_b,
            [0.0, 0.0],
            EPSILON_MIN_F32
        ));
    }
}

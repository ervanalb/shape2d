use std::cmp::Ordering;

const EPSILON_MIN_F32: f32 = 1e-5;
const EPSILON_RATE_F32: f32 = 1e-5;

use crate::{
    kernel::{EdgeSide, Kernel, VertexEvent},
    rtree::Rect,
    sweep_line::{SweepLineChain, SweepLineEvent, SweepLineEventType, SweepLineSegment},
    triangle_kernel::{TriangleKernel, TriangleKernelF32},
};

/// Trait for providing epsilon values to the F32 kernel
pub trait EpsilonProviderF32 {
    /// Returns the epsilon value to use for geometric comparisons
    fn value(&self, fp_mag: f32) -> f32;
}

/// Implementation for f32
impl EpsilonProviderF32 for f32 {
    #[inline]
    fn value(&self, fp_mag: f32) -> f32 {
        self.max(EPSILON_RATE_F32 * fp_mag)
    }
}

/// Constant epsilon provider - returns a compile-time constant
#[derive(Clone, Debug)]
pub struct DefaultEpsilon;

impl EpsilonProviderF32 for DefaultEpsilon {
    #[inline]
    fn value(&self, fp_mag: f32) -> f32 {
        EPSILON_MIN_F32.max(EPSILON_RATE_F32 * fp_mag)
    }
}

#[derive(Debug, Clone)]
pub struct F32<E: EpsilonProviderF32 = DefaultEpsilon> {
    pub vertices: Vec<[f32; 2]>,
    pub epsilon: E,
}

impl<E: EpsilonProviderF32> F32<E> {
    pub fn v(&self, i: u32) -> [f32; 2] {
        self.vertices[i as usize]
    }

    fn push_vertex(&mut self, pt: [f32; 2]) -> u32 {
        let i = self.vertices.len();
        self.vertices.push(pt);
        i as u32
    }
}

impl F32<DefaultEpsilon> {
    pub fn new() -> Self {
        Self {
            vertices: vec![],
            epsilon: DefaultEpsilon,
        }
    }

    pub fn new_with_vertices(vertices: Vec<[f32; 2]>) -> Self {
        Self {
            vertices,
            epsilon: DefaultEpsilon,
        }
    }
}

impl<E: EpsilonProviderF32> F32<E> {
    pub fn new_with_vertices_and_epsilon(vertices: Vec<[f32; 2]>, epsilon: E) -> Self {
        Self { vertices, epsilon }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum CapStyleF32 {
    Arc { tolerance: f32 },
    Bevel,
    Miter { limit: f32 },
}

impl<E: EpsilonProviderF32> Kernel for F32<E> {
    type Vertex = u32;
    type Edge = (u32, u32);
    type Extents = ExtentsF32;
    type Point = [f32; 2];
    type SweepLineEdgePortion = ();
    type SweepLineEventPoint = u32;
    type TriangleKernel = TriangleKernelF32;
    type CapStyle = CapStyleF32;
    type OffsetAmount = f32;

    fn vertices_coincident(&self, a: Self::Vertex, b: Self::Vertex) -> bool {
        let pt_a = self.v(a);
        let pt_b = self.v(b);
        let fp_mag = fp_mag_pt_f32(pt_a).max(fp_mag_pt_f32(pt_b));
        points_coincident_f32(self.v(a), self.v(b), self.epsilon.value(fp_mag))
    }

    fn edges_coincident(&self, _a: Self::Edge, _b: Self::Edge) -> bool {
        // Line segments will never be coincident unless they share endpoints,
        // in which case they will simply be equal
        false
    }

    fn vertex_on_edge(&self, vertex: Self::Vertex, edge: Self::Edge) -> Option<Self::Point> {
        let vertex_pt = self.v(vertex);
        let edge_start_pt = self.v(edge.0);
        let edge_end_pt = self.v(edge.1);
        let fp_mag = fp_mag_pt_f32(vertex_pt)
            .max(fp_mag_pt_f32(edge_start_pt))
            .max(fp_mag_pt_f32(edge_end_pt));

        point_on_segment_f32(
            vertex_pt,
            edge_start_pt,
            edge_end_pt,
            self.epsilon.value(fp_mag),
        )
    }

    fn intersection(&self, a: Self::Edge, b: Self::Edge) -> Option<Self::Point> {
        if a.0 == b.0 || a.0 == b.1 || a.1 == b.0 || a.1 == b.1 {
            // Segments that share an endpoint don't intersect
            return None;
        }

        let pt_a0 = self.v(a.0);
        let pt_a1 = self.v(a.1);
        let pt_b0 = self.v(b.0);
        let pt_b1 = self.v(b.1);
        let fp_mag = fp_mag_pt_f32(pt_a0)
            .max(fp_mag_pt_f32(pt_a1))
            .max(fp_mag_pt_f32(pt_b0))
            .max(fp_mag_pt_f32(pt_b1));

        Some(intersect_segments_f32(
            self.v(a.0),
            self.v(a.1),
            self.v(b.0),
            self.v(b.1),
            self.epsilon.value(fp_mag),
        )?)
    }

    fn push_vertex(&mut self, intersection: Self::Point) -> Self::Vertex {
        self.push_vertex(intersection)
    }

    fn merged_vertex(&mut self, a: Self::Vertex, b: Self::Vertex) -> Self::Vertex {
        self.push_vertex(merge_points_f32(self.v(a), self.v(b)))
    }

    fn merged_edges(&mut self, _a: Self::Edge, _b: Self::Edge) -> (Self::Edge, Self::Edge) {
        panic!("Not possible to merge line segments");
    }

    fn extents(&self, edges: impl Iterator<Item = Self::Edge>) -> Self::Extents {
        extents_f32(
            edges.flat_map(|(a, b)| [self.v(a), self.v(b)]),
            &self.epsilon,
        )
    }

    fn edge_bbox(&self, edge: Self::Edge, extents: &Self::Extents) -> Rect {
        let edge_start_pt = self.v(edge.0);
        let edge_end_pt = self.v(edge.1);
        let fp_mag = fp_mag_pt_f32(edge_start_pt).max(fp_mag_pt_f32(edge_end_pt));
        segment_bbox_f32(
            edge_start_pt,
            edge_end_pt,
            *extents,
            self.epsilon.value(fp_mag),
        )
    }

    fn sin_cmp(&self, common: Self::Vertex, a: Self::Vertex, b: Self::Vertex) -> Ordering {
        sin_cmp_f32(self.v(common), self.v(a), self.v(b))
    }

    fn vertices_for_edge(&self, edge: Self::Edge) -> Option<(Self::Vertex, Self::Vertex)> {
        Some((edge.0, edge.1))
    }

    fn replace_vertex_in_edge(
        &mut self,
        mut edge: Self::Edge,
        old_v: Self::Vertex,
        new_v: Self::Vertex,
    ) -> Option<Self::Edge> {
        if edge.0 == old_v {
            if edge.1 == new_v {
                return None; // Reflex edge
            }
            edge.0 = new_v;
        }
        if edge.1 == old_v {
            if edge.0 == new_v {
                return None; // Reflex edge
            }
            edge.1 = new_v;
        }
        Some(edge)
    }

    fn split_edge(&self, edge: Self::Edge, vertex: Self::Vertex) -> (Self::Edge, Self::Edge) {
        ((edge.0, vertex), (vertex, edge.1))
    }

    fn sweep_line_events_for_edge(
        &self,
        edge: Self::Edge,
    ) -> impl Iterator<Item = SweepLineEvent<Self>> {
        // Edges always have exactly 2 events.
        // We will use the `segment` data to store whether this edge is bottom (going  right) or top (going left)
        // based on how its endpoints sort in sweep-line order.
        let chain = match sweep_line_cmp_f32(self.v(edge.0), self.v(edge.1)) {
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
        let a_pt = self.v(sweep_line_select_vertex(
            a.event_type,
            a.segment.chain,
            a.segment.edge,
        ));
        let b_pt = self.v(sweep_line_select_vertex(
            b.event_type,
            b.segment.chain,
            b.segment.edge,
        ));

        // Compare first by event point (sweep-line order)
        sweep_line_cmp_f32(a_pt, b_pt)
            // Then by event type (End before Start)
            .then_with(|| a.event_type.cmp(&b.event_type))
            // Then by incidence angle, bottom-to-top
            .then_with(|| {
                let shared_event_type = a.event_type;
                let shared_pt = a_pt;
                let a_other_pt = self.v(sweep_line_select_vertex(
                    shared_event_type.other(),
                    a.segment.chain,
                    a.segment.edge,
                ));
                let b_other_pt = self.v(sweep_line_select_vertex(
                    shared_event_type.other(),
                    b.segment.chain,
                    b.segment.edge,
                ));
                match shared_event_type {
                    SweepLineEventType::End => sin_cmp_f32(shared_pt, b_other_pt, a_other_pt),
                    SweepLineEventType::Start => sin_cmp_f32(shared_pt, a_other_pt, b_other_pt),
                }
            })
    }

    fn sweep_line_event_point(&self, event: &SweepLineEvent<Self>) -> Self::SweepLineEventPoint {
        sweep_line_select_vertex(event.event_type, event.segment.chain, event.segment.edge)
    }

    fn sweep_line_segment_cmp(
        &self,
        segment: &SweepLineSegment<Self>,
        event_point: Self::SweepLineEventPoint,
    ) -> Ordering {
        let (left_i, right_i) = match segment.chain {
            SweepLineChain::Bottom => (segment.edge.0, segment.edge.1),
            SweepLineChain::Top => (segment.edge.1, segment.edge.0),
        };
        let left_pt = self.v(left_i);
        let right_pt = self.v(right_i);
        let common_pt = self.v(event_point);
        sin_cmp_f32(left_pt, right_pt, common_pt)
    }

    fn sweep_line_event_point_to_triangle_vertex(
        &self,
        triangle_kernel: &mut Self::TriangleKernel,
        event_point: Self::SweepLineEventPoint,
    ) -> <Self::TriangleKernel as TriangleKernel>::Vertex {
        triangle_kernel.push_vertex(self.v(event_point))
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
        let a_pt = self.v(select_vertex(a.event_type, a.edge));
        let b_pt = self.v(select_vertex(b.event_type, b.edge));

        // Compare first by event point (arbitrary; we'll pick sweep-line order)
        sweep_line_cmp_f32(a_pt, b_pt)
            // Then by incidence angle, we will have +X direction be the first,
            // and go CCW from there
            .then_with(|| {
                let shared_pt = a_pt;
                let a_other_pt = self.v(select_vertex(a.event_type.other(), a.edge));
                let b_other_pt = self.v(select_vertex(b.event_type.other(), b.edge));

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

    fn offset_edge(&mut self, edge: Self::Edge, offset: Self::OffsetAmount) -> Self::Edge {
        let a_pt = self.v(edge.0);
        let b_pt = self.v(edge.1);

        let offset_x = b_pt[1] - a_pt[1]; // Y
        let offset_y = a_pt[0] - b_pt[0]; // -X
        let norm_factor = offset / offset_x.hypot(offset_y);
        let offset_x = offset_x * norm_factor;
        let offset_y = offset_y * norm_factor;

        let new_a = self.push_vertex([a_pt[0] + offset_x, a_pt[1] + offset_y]);
        let new_b = self.push_vertex([b_pt[0] + offset_x, b_pt[1] + offset_y]);

        (new_a, new_b)
    }

    fn cap_corner(
        &mut self,
        incoming_edge: Self::Edge,
        outgoing_edge: Self::Edge,
        original_vertex: Self::Vertex,
        offset_amount: Self::OffsetAmount,
        cap_style: Self::CapStyle,
        mut emit_edge: impl FnMut(Self::Edge),
    ) {
        let ix_a = incoming_edge.1;
        let ix_b = outgoing_edge.0;
        let pt_a = self.v(ix_a);
        let pt_b = self.v(ix_b);
        let original_pt = self.v(original_vertex);

        let fp_mag = fp_mag_pt_f32(pt_a)
            .max(fp_mag_pt_f32(pt_b))
            .max(fp_mag_pt_f32(original_pt));

        // For an interior corner, draw lines to the original vertex
        // as per https://mcmains.me.berkeley.edu/pubs/DAC05OffsetPolygon.pdf
        if matches!(
            (sin_cmp_f32(original_pt, pt_a, pt_b), offset_amount >= 0.),
            (Ordering::Greater, true) | (Ordering::Less, false)
        ) {
            emit_edge((ix_a, original_vertex));
            emit_edge((original_vertex, ix_b));
            return;
        }

        // If the two vertices are coincident, connect them with a bevel
        if points_coincident_f32(pt_a, pt_b, self.epsilon.value(fp_mag)) {
            emit_edge((ix_a, ix_b));
            return;
        }

        match cap_style {
            CapStyleF32::Arc { tolerance } => {
                // Calculate vectors from center to both points
                let dx_a = pt_a[0] - original_pt[0];
                let dy_a = pt_a[1] - original_pt[1];
                let dx_b = pt_b[0] - original_pt[0];
                let dy_b = pt_b[1] - original_pt[1];

                // Calculate the angle between the two vectors
                let cross = dx_a * dy_b - dy_a * dx_b; // 2D cross product (z-component)
                let dot = dx_a * dx_b + dy_a * dy_b; // dot product
                let delta_angle = cross.atan2(dot);

                // Calculate maximum angle per segment based on tolerance
                // The sagitta (deviation) of a chord from an arc is: s = r * (1 - cos(θ/2))
                // Solving for θ when s = tolerance: θ = 2 * arccos(1 - tolerance/r)
                let cos_arg = (1.0
                    - tolerance.max(self.epsilon.value(fp_mag)) / offset_amount.abs())
                .clamp(-1.0, 1.0);
                let max_angle_per_segment = 2.0 * cos_arg.acos();
                let num_segments =
                    (delta_angle.abs() / max_angle_per_segment).ceil().max(1.) as u32;

                // Interpolate the arc
                let mut prev_vertex = ix_a;
                let factor = delta_angle / num_segments as f32;
                for i in 1..num_segments {
                    let theta = i as f32 * factor;
                    let c = theta.cos();
                    let s = theta.sin();
                    // Rotate [dx_a, dy_a] by theta
                    let dx = c * dx_a - s * dy_a;
                    let dy = s * dx_a + c * dy_a;
                    let v = self.push_vertex([original_pt[0] + dx, original_pt[1] + dy]);
                    emit_edge((prev_vertex, v));
                    prev_vertex = v;
                }
                emit_edge((prev_vertex, ix_b));
            }
            CapStyleF32::Bevel => {
                emit_edge((ix_a, ix_b));
            }
            CapStyleF32::Miter { limit } => {
                let limit = limit * offset_amount.abs();
                let before_pt_a = self.v(incoming_edge.0);
                let after_pt_b = self.v(outgoing_edge.1);

                let miter_pt = intersect_lines_f32(before_pt_a, pt_a, pt_b, after_pt_b);
                let d_a = [miter_pt[0] - pt_a[0], miter_pt[1] - pt_a[1]];
                let d_a_sq = d_a[0] * d_a[0] + d_a[1] * d_a[1];
                if d_a_sq < limit * limit {
                    // Emit a true miter
                    let miter_v = self.push_vertex(miter_pt);
                    emit_edge((ix_a, miter_v));
                    emit_edge((miter_v, ix_b));
                } else {
                    // Emit a clipped miter
                    // Extended edges `a` and `b` by `limit`

                    let a_vec = [pt_a[0] - before_pt_a[0], pt_a[1] - before_pt_a[1]];
                    let a_len = a_vec[0].hypot(a_vec[1]);
                    let a_factor = limit / a_len;
                    let extended_a_pt =
                        [pt_a[0] + a_vec[0] * a_factor, pt_a[1] + a_vec[1] * a_factor];
                    let extended_a_vertex = self.push_vertex(extended_a_pt);

                    let b_vec = [pt_b[0] - after_pt_b[0], pt_b[1] - after_pt_b[1]];
                    let b_len = b_vec[0].hypot(b_vec[1]);
                    let b_factor = limit / b_len;
                    let extended_b_pt =
                        [pt_b[0] + b_vec[0] * b_factor, pt_b[1] + b_vec[1] * b_factor];
                    let extended_b_vertex = self.push_vertex(extended_b_pt);

                    emit_edge((ix_a, extended_a_vertex));
                    emit_edge((extended_a_vertex, extended_b_vertex));
                    emit_edge((extended_b_vertex, ix_b));
                }
            }
        }
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
fn point_on_segment_f32(
    p: [f32; 2],
    segment_start: [f32; 2],
    segment_end: [f32; 2],
    epsilon: f32,
) -> Option<[f32; 2]> {
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
    Some([
        segment_start[0] + segment_dx * t,
        segment_start[1] + segment_dy * t,
    ])
}

#[inline]
fn intersect_segments_f32(
    a_start: [f32; 2],
    a_end: [f32; 2],
    b_start: [f32; 2],
    b_end: [f32; 2],
    _epsilon: f32, // XXX
) -> Option<[f32; 2]> {
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

fn extents_f32<E: EpsilonProviderF32>(
    mut points: impl Iterator<Item = [f32; 2]>,
    epsilon_provider: &E,
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
    let epsilon = epsilon_provider.value(fp_mag);

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
) -> [f32; 2] {
    // Promote to homogeneous coordinates
    let a_start_h = [a_start[0], a_start[1], 1.0];
    let a_end_h = [a_end[0], a_end[1], 1.0];
    let b_start_h = [b_start[0], b_start[1], 1.0];
    let b_end_h = [b_end[0], b_end[1], 1.0];

    // Lines using Plucker coordinates
    let a = cross_f32(a_start_h, a_end_h);
    let b = cross_f32(b_start_h, b_end_h);

    // Calculate intersection
    let intersection_h = cross_f32(a, b);

    // Divide out perspective coordinate
    let z_inv = 1. / intersection_h[2];
    let first_guess = [intersection_h[0] * z_inv, intersection_h[1] * z_inv];

    // Subtract the first guess from our coordinates and do it all again
    // to gain accuracy

    let a_start_h = [
        a_start[0] - first_guess[0],
        a_start[1] - first_guess[1],
        1.0,
    ];
    let a_end_h = [a_end[0] - first_guess[0], a_end[1] - first_guess[1], 1.0];
    let b_start_h = [
        b_start[0] - first_guess[0],
        b_start[1] - first_guess[1],
        1.0,
    ];
    let b_end_h = [b_end[0] - first_guess[0], b_end[1] - first_guess[1], 1.0];
    let a = cross_f32(a_start_h, a_end_h);
    let b = cross_f32(b_start_h, b_end_h);
    let intersection_h = cross_f32(a, b);
    let z_inv = 1. / intersection_h[2];
    let intersection = [intersection_h[0] * z_inv, intersection_h[1] * z_inv];

    [
        first_guess[0] + intersection[0],
        first_guess[1] + intersection[1],
    ]
}

#[inline]
fn cross_f32(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
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
        assert!(points_coincident_f32(
            point_on_segment_f32(mid, start, end, EPSILON_MIN_F32).unwrap(),
            mid,
            EPSILON_MIN_F32
        ));
    }

    #[test]
    fn test_point_on_segment_at_start() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.001, 0.0];
        assert!(points_coincident_f32(
            point_on_segment_f32(p, start, end, EPSILON_MIN_F32).unwrap(),
            p,
            EPSILON_MIN_F32,
        ));
    }

    #[test]
    fn test_point_on_segment_at_end() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.999, 0.0];
        assert!(points_coincident_f32(
            point_on_segment_f32(p, start, end, EPSILON_MIN_F32).unwrap(),
            p,
            EPSILON_MIN_F32,
        ));
    }

    #[test]
    fn test_point_on_segment_not_on_line() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.5, 0.1];
        assert!(point_on_segment_f32(p, start, end, EPSILON_MIN_F32).is_none());
    }

    #[test]
    fn test_point_on_segment_beyond_end() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [1.5, 0.0];
        assert!(point_on_segment_f32(p, start, end, EPSILON_MIN_F32).is_none());
    }

    #[test]
    fn test_point_on_segment_before_start() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [-0.5, 0.0];
        assert!(point_on_segment_f32(p, start, end, EPSILON_MIN_F32).is_none());
    }

    #[test]
    fn test_point_on_segment_diagonal() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 1.0];
        let mid = [0.5, 0.5];
        assert!(points_coincident_f32(
            point_on_segment_f32(mid, start, end, EPSILON_MIN_F32).unwrap(),
            mid,
            EPSILON_MIN_F32,
        ));
    }

    #[test]
    fn test_point_on_segment_within_epsilon_margin() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.5, 0.000001]; // Very close to the line
        assert!(points_coincident_f32(
            point_on_segment_f32(p, start, end, EPSILON_MIN_F32).unwrap(),
            p,
            EPSILON_MIN_F32,
        ));
    }

    #[test]
    fn test_intersect_segments_crossing_lines() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 1.0];
        let b_start = [0.0, 1.0];
        let b_end = [1.0, 0.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON_MIN_F32);
        assert!(points_coincident_f32(
            result.unwrap(),
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

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON_MIN_F32);
        assert!(points_coincident_f32(
            result.unwrap(),
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

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON_MIN_F32);
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_no_overlap() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.6, 0.4];
        let b_end = [1.0, 0.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON_MIN_F32);
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_collinear() {
        // Two collinear segments that share an endpoint
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.5, 0.5];
        let b_end = [1.0, 1.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON_MIN_F32);
        // Collinear segments return None (det is too small)
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_almost_parallel() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 0.0];
        let b_start = [0.0, 0.0];
        let b_end = [1.0, 0.000001]; // Almost parallel

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON_MIN_F32);
        // Should be None because det is too small
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_close() {
        let a0 = [2.177941, 2.2992249];
        let a1 = [2.1857128, 2.2375135];
        let b0 = [2.1857483, 2.2371936];
        let b1 = [2.1862082, 2.2374542];

        let intersection = intersect_segments_f32(a0, a1, b0, b1, EPSILON_MIN_F32);
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
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), &DefaultEpsilon);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON_MIN_F32);

        assert!(bbox.min[0] < bbox.max[0]);
        assert!(bbox.min[1] <= bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bbox_vertical() {
        let start = [0.5_f32, 0.2];
        let end = [0.5, 0.8];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), &DefaultEpsilon);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON_MIN_F32);

        assert!(bbox.min[0] <= bbox.max[0]);
        assert!(bbox.min[1] < bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bbox_diagonal() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 1.0];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), &DefaultEpsilon);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON_MIN_F32);

        assert!(bbox.min[0] < bbox.max[0]);
        assert!(bbox.min[1] < bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bbox_point() {
        let start = [0.5_f32, 0.5];
        let end = [0.5, 0.5];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), &DefaultEpsilon);
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

        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), &DefaultEpsilon);
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

        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), &DefaultEpsilon);
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

        let result = intersect_lines_f32(a1, a2, b1, b2);
        assert!(points_coincident_f32(result, [0.5, 0.0], EPSILON_MIN_F32));
    }

    #[test]
    fn test_intersect_lines_diagonal() {
        let a1 = [0.0, 0.0];
        let a2 = [1.0, 1.0]; // Line y = x
        let b1 = [0.0, 1.0];
        let b2 = [1.0, 0.0]; // Line y = -x + 1

        let result = intersect_lines_f32(a1, a2, b1, b2);
        assert!(points_coincident_f32(result, [0.5, 0.5], EPSILON_MIN_F32));
    }

    #[test]
    fn test_intersect_lines_at_origin() {
        let a1 = [-1.0, 0.0];
        let a2 = [1.0, 0.0]; // Horizontal line through origin
        let b1 = [0.0, -1.0];
        let b2 = [0.0, 1.0]; // Vertical line through origin

        let result = intersect_lines_f32(a1, a2, b1, b2);
        assert!(points_coincident_f32(result, [0.0, 0.0], EPSILON_MIN_F32));
    }

    #[test]
    fn test_intersect_lines_various_angles() {
        // Line 1: from (0,0) to (2,1) → y = x/2
        let a1 = [0.0, 0.0];
        let a2 = [2.0, 1.0];
        // Line 2: from (0,1) to (2,0) → y = 1 - x/2
        let b1 = [0.0, 1.0];
        let b2 = [2.0, 0.0];

        let result = intersect_lines_f32(a1, a2, b1, b2);
        // These lines should intersect at (1, 0.5)
        assert!(points_coincident_f32(result, [1.0, 0.5], EPSILON_MIN_F32));
    }

    #[test]
    fn test_intersect_lines_negative_coords() {
        let a1 = [-2.0, -2.0];
        let a2 = [2.0, 2.0]; // Line through origin, slope 1
        let b1 = [-2.0, 2.0];
        let b2 = [2.0, -2.0]; // Line through origin, slope -1

        let result = intersect_lines_f32(a1, a2, b1, b2);
        assert!(points_coincident_f32(result, [0.0, 0.0], EPSILON_MIN_F32));
    }
}

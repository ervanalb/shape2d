use std::cmp::Ordering;

const DEFAULT_EPSILON_F32: f32 = 0.1; // XXX

use crate::{
    kernel::{Few, Kernel, SweepLineChain, SweepLineEvent, SweepLineEventType, SweepLineSegment},
    rtree::Rect,
    triangle_kernel::{F32TriangleKernel, TriangleKernel},
};

/// Trait for providing epsilon values to the F32 kernel
pub trait EpsilonProviderF32 {
    /// Returns the epsilon value to use for geometric comparisons
    fn value(&self) -> f32;
}

/// Implementation for f32 - epsilon() returns the value itself
impl EpsilonProviderF32 for f32 {
    #[inline]
    fn value(&self) -> f32 {
        *self
    }
}

/// Constant epsilon provider - returns a compile-time constant
#[derive(Debug)]
pub struct ConstantEpsilon;

impl EpsilonProviderF32 for ConstantEpsilon {
    #[inline]
    fn value(&self) -> f32 {
        DEFAULT_EPSILON_F32
    }
}

#[derive(Debug, Clone)]
pub struct F32<E: EpsilonProviderF32 = ConstantEpsilon> {
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

impl F32<ConstantEpsilon> {
    pub fn new(vertices: Vec<[f32; 2]>) -> Self {
        Self {
            vertices,
            epsilon: ConstantEpsilon,
        }
    }
}

impl<E: EpsilonProviderF32> F32<E> {
    pub fn new_with_epsilon(vertices: Vec<[f32; 2]>, epsilon: E) -> Self {
        Self { vertices, epsilon }
    }
}

impl<E: EpsilonProviderF32> Kernel for F32<E> {
    type Vertex = u32;
    type Edge = (u32, u32);
    type Extents = ExtentsF32;
    type Intersection = [f32; 2];
    type SweepLineEdgePortion = ();
    type SweepLineEventPoint = u32;
    type TriangleKernel = F32TriangleKernel;

    fn vertices_coincident(&self, a: Self::Vertex, b: Self::Vertex) -> bool {
        points_coincident_f32(self.v(a), self.v(b), self.epsilon.value())
    }

    fn edges_coincident(&self, _a: Self::Edge, _b: Self::Edge) -> bool {
        // Line segments will never be coincident unless they share endpoints,
        // in which case they will simply be equal
        false
    }

    fn vertex_on_edge(&self, vertex: Self::Vertex, edge: Self::Edge) -> bool {
        point_on_segment_f32(
            self.v(vertex),
            self.v(edge.0),
            self.v(edge.1),
            self.epsilon.value(),
        )
    }

    fn intersection(&self, a: Self::Edge, b: Self::Edge) -> Option<Self::Intersection> {
        if a.0 == b.0 || a.0 == b.1 || a.1 == b.0 || a.1 == b.1 {
            // Segments that share an endpoint don't intersect
            return None;
        }

        Some(intersect_segments_f32(
            self.v(a.0),
            self.v(a.1),
            self.v(b.0),
            self.v(b.1),
            self.epsilon.value(),
        )?)
    }

    fn intersection_vertex(&mut self, intersection: Self::Intersection) -> Self::Vertex {
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
            self.epsilon.value(),
        )
    }

    fn edge_bbox(&self, edge: Self::Edge, extents: &Self::Extents) -> Rect {
        segment_bbox_f32(
            self.v(edge.0),
            self.v(edge.1),
            *extents,
            self.epsilon.value(),
        )
    }

    fn sin_cmp(&self, common: Self::Vertex, a: Self::Vertex, b: Self::Vertex) -> Ordering {
        sin_cmp_f32(self.v(common), self.v(a), self.v(b))
    }

    fn vertices_for_edge(&self, edge: Self::Edge) -> Few<Self::Vertex> {
        if edge.0 == edge.1 {
            Few::One(edge.0)
        } else {
            Few::Two(edge.0, edge.1)
        }
    }

    fn replace_vertex_in_edge(
        &self,
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

    fn sweep_line_event_cmp_bottom_up(
        &self,
        a: &SweepLineEvent<Self>,
        b: &SweepLineEvent<Self>,
    ) -> Ordering {
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
                    SweepLineEventType::End => sin_cmp_f32(shared_pt, a_other_pt, b_other_pt),
                    SweepLineEventType::Start => sin_cmp_f32(shared_pt, b_other_pt, a_other_pt),
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
        sin_cmp_f32(left_pt, common_pt, right_pt)
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

    fn sweep_line_event_cmp_clockwise(
        &self,
        a: &SweepLineEvent<Self>,
        b: &SweepLineEvent<Self>,
    ) -> Ordering {
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
            // Then by incidence angle, clockwise
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
                // Clockwise sorting always uses the same order
                sin_cmp_f32(shared_pt, a_other_pt, b_other_pt)
            })
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct ExtentsF32 {
    scale: [f32; 2],
    offset: [f32; 2],
}

#[inline]
pub fn sweep_line_cmp_f32(a: [f32; 2], b: [f32; 2]) -> Ordering {
    a[0].partial_cmp(&b[0])
        .unwrap_or(Ordering::Equal)
        .then_with(|| a[1].partial_cmp(&b[1]).unwrap_or(Ordering::Equal))
}

#[inline]
pub fn sin_cmp_f32(common: [f32; 2], a: [f32; 2], b: [f32; 2]) -> Ordering {
    // Compute cross product: (a - origin) x (b - origin)
    let ax = a[0] - common[0];
    let ay = a[1] - common[1];
    let bx = b[0] - common[0];
    let by = b[1] - common[1];
    let cross = ax * by - ay * bx;

    cross.partial_cmp(&0.).unwrap_or(Ordering::Equal)
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
) -> bool {
    let segment_dx = segment_end[0] - segment_start[0];
    let segment_dy = segment_end[1] - segment_start[1];
    let segment_len_sq = segment_dx * segment_dx + segment_dy * segment_dy;

    // Check if point is collinear with segment using cross product
    let to_p_x = p[0] - segment_start[0];
    let to_p_y = p[1] - segment_start[1];
    let cross = segment_dy * to_p_x - segment_dx * to_p_y;

    // Use squared epsilon for consistent margin
    if cross * cross >= epsilon * epsilon * segment_len_sq {
        return false;
    }

    // Check if point is within segment bounds using dot product
    let dot = to_p_x * segment_dx + to_p_y * segment_dy;
    dot >= 0. && dot <= segment_len_sq
}

#[inline]
fn intersect_segments_f32(
    a_start: [f32; 2],
    a_end: [f32; 2],
    b_start: [f32; 2],
    b_end: [f32; 2],
    epsilon: f32,
) -> Option<[f32; 2]> {
    let da_x = a_end[0] - a_start[0];
    let da_y = a_end[1] - a_start[1];
    let db_x = b_end[0] - b_start[0];
    let db_y = b_end[1] - b_start[1];

    let det = da_x * db_y - da_y * db_x;

    let dx3 = b_start[0] - a_start[0];
    let dy3 = b_start[1] - a_start[1];

    let ta_det = dx3 * db_y - dy3 * db_x;
    let tb_det = dx3 * da_y - dy3 * da_x;

    // Check if intersection is within both edge segments
    let s = det.signum(); // Ensure inequalities compare the right way
    if s * ta_det > 0. && s * ta_det < s * det && s * tb_det > 0. && s * tb_det < s * det {
        let inv_det = 1. / det;

        // Compute the intersection point two different ways and make sure they agree
        // (they might disagree if inv_det is huge)
        let pt_a = [
            a_start[0] + ta_det * da_x * inv_det,
            a_start[1] + ta_det * da_y * inv_det,
        ];
        let pt_b = [
            b_start[0] + tb_det * db_x * inv_det,
            b_start[1] + tb_det * db_y * inv_det,
        ];

        if points_coincident_f32(pt_a, pt_b, epsilon) {
            Some(pt_a)
        } else {
            None
        }
    } else {
        None
    }
}

#[inline]
fn merge_points_f32(a: [f32; 2], b: [f32; 2]) -> [f32; 2] {
    [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])]
}

#[inline]
fn extents_f32(mut points: impl Iterator<Item = [f32; 2]>, epsilon: f32) -> ExtentsF32 {
    let Some(first) = points.next() else {
        return Default::default();
    };
    let (min, max) = points.fold((first, first), |(min, max), point| {
        (
            [
                min[0].min(point[0] - epsilon),
                min[1].min(point[1] - epsilon),
            ],
            [
                max[0].max(point[0] + epsilon),
                max[1].max(point[1] + epsilon),
            ],
        )
    });

    let scale = [
        u16::MAX as f32 / (max[0] - min[0]),
        u16::MAX as f32 / (max[1] - min[1]),
    ];
    let offset = [-min[0] * scale[0], -min[1] * scale[1]];

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

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use super::*;

    const EPSILON: f32 = 1e-5;

    #[test]
    fn test_points_coincident_same_point() {
        let a = [0.5_f32, 0.5];
        let b = [0.5, 0.5];
        assert!(points_coincident_f32(a, b, EPSILON));
    }

    #[test]
    fn test_points_coincident_within_epsilon() {
        let a = [0.5_f32, 0.5];
        let b = [0.500001, 0.500001];
        assert!(points_coincident_f32(a, b, EPSILON));
    }

    #[test]
    fn test_points_coincident_outside_epsilon() {
        let a = [0.5_f32, 0.5];
        let b = [0.51, 0.5];
        assert!(!points_coincident_f32(a, b, EPSILON));
    }

    #[test]
    fn test_points_coincident_far_apart() {
        let a = [0.0_f32, 0.0];
        let b = [1.0, 1.0];
        assert!(!points_coincident_f32(a, b, EPSILON));
    }

    #[test]
    fn test_point_on_segment_midpoint() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let mid = [0.5, 0.0];
        assert!(point_on_segment_f32(mid, start, end, EPSILON));
    }

    #[test]
    fn test_point_on_segment_at_start() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.0, 0.0];
        assert!(point_on_segment_f32(p, start, end, EPSILON));
    }

    #[test]
    fn test_point_on_segment_at_end() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [1.0, 0.0];
        assert!(point_on_segment_f32(p, start, end, EPSILON));
    }

    #[test]
    fn test_point_on_segment_not_on_line() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.5, 0.1];
        assert!(!point_on_segment_f32(p, start, end, EPSILON));
    }

    #[test]
    fn test_point_on_segment_beyond_end() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [1.5, 0.0];
        assert!(!point_on_segment_f32(p, start, end, EPSILON));
    }

    #[test]
    fn test_point_on_segment_before_start() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [-0.5, 0.0];
        assert!(!point_on_segment_f32(p, start, end, EPSILON));
    }

    #[test]
    fn test_point_on_segment_diagonal() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 1.0];
        let mid = [0.5, 0.5];
        assert!(point_on_segment_f32(mid, start, end, EPSILON));
    }

    #[test]
    fn test_point_on_segment_within_epsilon_margin() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.5, 0.000001]; // Very close to the line
        assert!(point_on_segment_f32(p, start, end, EPSILON));
    }

    #[test]
    fn test_intersect_segments_crossing_lines() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 1.0];
        let b_start = [0.0, 1.0];
        let b_end = [1.0, 0.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON);
        assert!(points_coincident_f32(result.unwrap(), [0.5, 0.5], EPSILON));
    }

    #[test]
    fn test_intersect_segments_t_junction() {
        let a_start = [0.0_f32, 0.5];
        let a_end = [1.0, 0.5];
        let b_start = [0.5, 0.0];
        let b_end = [0.5, 1.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON);
        assert!(points_coincident_f32(result.unwrap(), [0.5, 0.5], EPSILON));
    }

    #[test]
    fn test_intersect_segments_parallel() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 0.0];
        let b_start = [0.0, 1.0];
        let b_end = [1.0, 1.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON);
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_no_overlap() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.6, 0.4];
        let b_end = [1.0, 0.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON);
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_collinear() {
        // Two collinear segments that share an endpoint
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.5, 0.5];
        let b_end = [1.0, 1.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON);
        // Collinear segments return None (det is too small)
        assert!(result.is_none());
    }

    #[test]
    fn test_intersect_segments_almost_parallel() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 0.0];
        let b_start = [0.0, 0.0];
        let b_end = [1.0, 0.000001]; // Almost parallel

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON);
        // Should be None because det is too small
        assert!(result.is_none());
    }

    #[test]
    fn test_merge_points_simple() {
        let a = [0.0_f32, 0.0];
        let b = [1.0, 1.0];
        let merged = merge_points_f32(a, b);
        assert!(points_coincident_f32(merged, [0.5, 0.5], EPSILON));
    }

    #[test]
    fn test_merge_points_same_point() {
        let a = [0.5_f32, 0.5];
        let b = [0.5, 0.5];
        let merged = merge_points_f32(a, b);
        assert!(points_coincident_f32(merged, [0.5, 0.5], EPSILON));
    }

    #[test]
    fn test_merge_points_negative_coords() {
        let a = [-1.0_f32, -1.0];
        let b = [1.0, 1.0];
        let merged = merge_points_f32(a, b);
        assert!(points_coincident_f32(merged, [0.0, 0.0], EPSILON));
    }

    #[test]
    fn test_segment_bbox_horizontal() {
        let start = [0.2_f32, 0.5];
        let end = [0.8, 0.5];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), EPSILON);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON);

        assert!(bbox.min[0] < bbox.max[0]);
        assert!(bbox.min[1] <= bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bbox_vertical() {
        let start = [0.5_f32, 0.2];
        let end = [0.5, 0.8];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), EPSILON);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON);

        assert!(bbox.min[0] <= bbox.max[0]);
        assert!(bbox.min[1] < bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bbox_diagonal() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 1.0];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), EPSILON);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON);

        assert!(bbox.min[0] < bbox.max[0]);
        assert!(bbox.min[1] < bbox.max[1]);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_segment_bbox_point() {
        let start = [0.5_f32, 0.5];
        let end = [0.5, 0.5];
        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), EPSILON);
        let bbox = segment_bbox_f32(start, end, extents, EPSILON);

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

        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), EPSILON);
        let bbox1 = segment_bbox_f32(start1, end1, extents, EPSILON);
        let bbox2 = segment_bbox_f32(start2, end2, extents, EPSILON);

        assert!(bbox1.overlaps(&bbox2));
    }

    #[test]
    fn test_segment_bboxes_dont_overlap() {
        let start1 = [0.1_f32, 0.1];
        let end1 = [0.2, 0.2];

        let start2 = [0.2_f32, 0.4];
        let end2 = [0.1, 0.5];

        let extents = extents_f32([[0., 0.], [1., 1.]].into_iter(), EPSILON);
        let bbox1 = segment_bbox_f32(start1, end1, extents, EPSILON);
        let bbox2 = segment_bbox_f32(start2, end2, extents, EPSILON);

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
        assert_eq!(sin_cmp_f32(origin, right, up), Ordering::Greater);
        assert_eq!(sin_cmp_f32(origin, up, right), Ordering::Less);

        // Collinear points
        assert_eq!(sin_cmp_f32(origin, diagonal, diagonal2), Ordering::Equal);
    }
}

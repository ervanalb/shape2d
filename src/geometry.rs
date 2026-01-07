use crate::Rect;
use std::cmp::Ordering;

pub trait Geometry {
    type Vertex: Copy + Ord;
    type Edge: Copy + Ord;
    type Extents;
    type Intersection;

    /// Check if two vertices are coincident (at the same location)
    fn vertices_coincident(&self, a: Self::Vertex, b: Self::Vertex) -> bool;

    /// Check if this vertex lies on the edge from edge_start to edge_end
    fn vertex_on_edge(&self, vertex: Self::Vertex, edge: Self::Edge) -> bool;

    /// See if two edges intersect
    fn intersection(&self, a: Self::Edge, b: Self::Edge) -> Option<Self::Intersection>;

    /// Merge two vertices, returning the merged result
    fn merged_vertex(&mut self, a: Self::Vertex, b: Self::Vertex) -> Self::Vertex;

    /// Creates and returns the vertex for an intersection
    fn intersection_vertex(&mut self, intersection: Self::Intersection) -> Self::Vertex;

    /// Generate an "extents" object from a list of edges,
    /// which is used when calculating the edge_bbox
    fn extents(&self, edges: impl Iterator<Item = Self::Edge>) -> Self::Extents;

    /// Compute the axis-aligned bounding box of an edge
    fn edge_bbox(&self, edge: Self::Edge, extents: &Self::Extents) -> Rect;

    /// Compare vertices in sweep line order (first by X, then by Y)
    fn sweep_line_cmp(&self, a: Self::Vertex, b: Self::Vertex) -> Ordering;

    /// Compare the angular order of two vectors from self to a and self to b.
    /// Returns the sign of the cross product: (a - self) x (b - self).
    /// Greater = b is counterclockwise from a (positive cross product)
    /// Equal = collinear (zero cross product)
    /// Less = b is clockwise from a (negative cross product)
    fn sin_cmp(&self, common: Self::Vertex, a: Self::Vertex, b: Self::Vertex) -> Ordering;

    // Returns the vertices which are endpoints for this edge
    fn vertices_for_edge(&self, edge: Self::Edge) -> FewVertices<Self::Vertex>;

    // Replaces instances of old_v with new_v in edge
    fn replace_vertex_in_edge(
        &self,
        edge: &mut Self::Edge,
        old_v: Self::Vertex,
        new_v: Self::Vertex,
    );

    // Splits an edge at the given vertex, returning two new edges
    fn split_edge(&self, edge: Self::Edge, vertex: Self::Vertex) -> (Self::Edge, Self::Edge);
}

pub enum FewVertices<V> {
    Zero,
    One(V),
    Two(V, V),
}

impl<V: Copy> Iterator for FewVertices<V> {
    type Item = V;

    fn next(&mut self) -> Option<Self::Item> {
        let (result, next) = match self {
            Self::Zero => (None, Self::Zero),
            Self::One(a) => (Some(*a), Self::Zero),
            Self::Two(a, b) => (Some(*a), Self::One(*b)),
        };
        *self = next;
        result
    }
}

struct MyGeometry {
    vertices: Vec<[f32; 2]>,
}

impl MyGeometry {
    const EPSILON: f32 = 1e-5;

    fn v(&self, i: u32) -> [f32; 2] {
        self.vertices[i as usize]
    }

    fn push_vertex(&mut self, pt: [f32; 2]) -> u32 {
        let i = self.vertices.len();
        self.vertices.push(pt);
        i as u32
    }
}

impl Geometry for MyGeometry {
    type Vertex = u32;
    type Edge = (u32, u32);
    type Extents = ExtentsF32;
    type Intersection = [f32; 2];

    fn vertices_coincident(&self, a: Self::Vertex, b: Self::Vertex) -> bool {
        points_coincident_f32(self.v(a), self.v(b), Self::EPSILON)
    }

    fn vertex_on_edge(&self, vertex: Self::Vertex, edge: Self::Edge) -> bool {
        point_on_segment_f32(
            self.v(vertex),
            self.v(edge.0),
            self.v(edge.1),
            Self::EPSILON,
        )
    }

    fn intersection(&self, a: Self::Edge, b: Self::Edge) -> Option<Self::Intersection> {
        Some(intersect_segments_f32(
            self.v(a.0),
            self.v(a.1),
            self.v(b.0),
            self.v(b.1),
            Self::EPSILON,
        )?)
    }

    fn intersection_vertex(&mut self, intersection: Self::Intersection) -> Self::Vertex {
        self.push_vertex(intersection)
    }

    fn merged_vertex(&mut self, a: Self::Vertex, b: Self::Vertex) -> Self::Vertex {
        self.push_vertex(merge_points_f32(self.v(a), self.v(b)))
    }

    fn extents(&self, edges: impl Iterator<Item = Self::Edge>) -> Self::Extents {
        extents_f32(
            edges.flat_map(|(a, b)| [self.v(a), self.v(b)]),
            Self::EPSILON,
        )
    }

    fn edge_bbox(&self, edge: Self::Edge, extents: &Self::Extents) -> Rect {
        segment_bbox_f32(self.v(edge.0), self.v(edge.1), *extents, Self::EPSILON)
    }

    fn sweep_line_cmp(&self, a: Self::Vertex, b: Self::Vertex) -> Ordering {
        sweep_line_cmp_f32(self.v(a), self.v(b))
    }

    fn sin_cmp(&self, common: Self::Vertex, a: Self::Vertex, b: Self::Vertex) -> Ordering {
        sin_cmp_f32(self.v(common), self.v(a), self.v(b))
    }

    fn vertices_for_edge(&self, edge: Self::Edge) -> FewVertices<Self::Vertex> {
        if edge.0 == edge.1 {
            FewVertices::One(edge.0)
        } else {
            FewVertices::Two(edge.0, edge.1)
        }
    }

    fn replace_vertex_in_edge(
        &self,
        edge: &mut Self::Edge,
        old_v: Self::Vertex,
        new_v: Self::Vertex,
    ) {
        if edge.0 == old_v {
            edge.0 = new_v;
        }
        if edge.1 == old_v {
            edge.1 = new_v;
        }
    }

    fn split_edge(&self, edge: Self::Edge, vertex: Self::Vertex) -> (Self::Edge, Self::Edge) {
        ((edge.0, vertex), (vertex, edge.1))
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct ExtentsF32 {
    scale: [f32; 2],
    offset: [f32; 2],
}

#[inline]
fn sweep_line_cmp_f32(a: [f32; 2], b: [f32; 2]) -> Ordering {
    a[0].partial_cmp(&b[0])
        .unwrap_or(Ordering::Equal)
        .then_with(|| a[1].partial_cmp(&b[1]).unwrap_or(Ordering::Equal))
}

#[inline]
fn sin_cmp_f32(common: [f32; 2], a: [f32; 2], b: [f32; 2]) -> Ordering {
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
    let dx1 = a_end[0] - a_start[0];
    let dy1 = a_end[1] - a_start[1];
    let dx2 = b_end[0] - b_start[0];
    let dy2 = b_end[1] - b_start[1];

    let det = dx1 * dy2 - dy1 * dx2;

    // Check if parallel/coincident (need abs() which requires conversion)
    // TODO(Eric): Find a way to remove this epsilon? or use a different epsilon?
    if det.abs() < epsilon {
        return None;
    }

    let dx3 = b_start[0] - a_start[0];
    let dy3 = b_start[1] - a_start[1];

    // TODO(Eric): Use multiplication instead of division where possible
    let det_inv = 1. / det;
    let t1 = (dx3 * dy2 - dy3 * dx2) * det_inv;
    let t2 = (dx3 * dy1 - dy3 * dx1) * det_inv;

    // Check if intersection is within both edge segments
    if t1 >= 0. && t1 <= 1. && t2 >= 0. && t2 <= 1. {
        Some([a_start[0] + t1 * dx1, a_start[1] + t1 * dy1])
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

#[cfg(test)]
mod tests {
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
    fn test_intersect_segments_at_endpoint() {
        // Two segments that share an endpoint but cross
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.5, 0.5];
        let b_end = [1.0, 0.0];

        let result = intersect_segments_f32(a_start, a_end, b_start, b_end, EPSILON);
        assert!(points_coincident_f32(result.unwrap(), [0.5, 0.5], EPSILON));
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

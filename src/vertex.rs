use crate::Rect;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Trait that vertex types must implement to be used with the clean algorithm
pub trait Vertex: Clone {
    type Extents: Default;

    /// Check if two vertices are coincident (at the same location)
    fn is_coincident(&self, other: &Self) -> bool;

    /// Check if this vertex lies on the edge from edge_start to edge_end
    fn is_on_edge(&self, edge_start: &Self, edge_end: &Self) -> bool;

    /// Compute the intersection point of two edges, if they intersect
    fn from_intersection(
        a_start: &Self,
        a_end: &Self,
        b_start: &Self,
        b_end: &Self,
    ) -> Option<Self>;

    /// Merge this vertex with another, returning the merged result
    fn merged_with(&self, other: &Self) -> Self;

    /// Generate an "extents" object from a list of edges,
    /// which is used when calculating the edge_bbox
    fn extents<'a>(edges: impl Iterator<Item = (&'a Self, &'a Self)>) -> Self::Extents
    where
        Self: 'a;

    /// Compute the axis-aligned bounding box of an edge
    fn edge_bbox(edge_start: &Self, edge_end: &Self, extents: &Self::Extents) -> Rect;
}

/// Trait for epsilon-based geometric operations on [T; 2] vertices
/// The implementing type (Self) represents the epsilon tolerance value
pub trait Epsilon: Copy {
    type Extents;

    /// Check if two [Self; 2] vertices are coincident (at the same location)
    fn is_coincident(self, a: &[Self; 2], b: &[Self; 2]) -> bool;

    /// Check if vertex p lies on the edge from edge_start to edge_end
    fn is_on_edge(self, p: &[Self; 2], edge_start: &[Self; 2], edge_end: &[Self; 2]) -> bool;

    /// Compute the intersection point of two edges, if they intersect
    fn intersect(
        self,
        a_start: &[Self; 2],
        a_end: &[Self; 2],
        b_start: &[Self; 2],
        b_end: &[Self; 2],
    ) -> Option<[Self; 2]>;

    /// Merge two vertices by averaging their coordinates, returning the merged result
    fn merged_with(self, a: &[Self; 2], b: &[Self; 2]) -> [Self; 2];

    /// Generate an "extents" object from a list of edges,
    /// which is used when calculating the edge_bbox
    fn extents(self, vertices: impl Iterator<Item = ([Self; 2], [Self; 2])>) -> Self::Extents;

    /// Compute the axis-aligned bounding box of an edge
    fn edge_bbox(
        self,
        edge_start: &[Self; 2],
        edge_end: &[Self; 2],
        extents: &Self::Extents,
    ) -> Rect;
}

trait MapToU16<T> {
    fn map(&self, pos: [T; 2]) -> [u16; 2];
}

/// Helper trait to provide float-specific specialization
trait Float:
    Copy
    + Default
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn half() -> Self;
    fn min(self, other: Self) -> Self;
    fn max(self, other: Self) -> Self;
    fn saturating_into_u16(self) -> u16;
    fn u16_max() -> Self;
}

impl Float for f32 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn half() -> Self {
        0.5
    }
    fn min(self, other: Self) -> Self {
        f32::min(self, other)
    }
    fn max(self, other: Self) -> Self {
        f32::max(self, other)
    }
    fn saturating_into_u16(self) -> u16 {
        self.max(0.).min(u16::MAX as f32) as u16
    }
    fn u16_max() -> Self {
        u16::MAX as f32
    }
}

impl Float for f64 {
    fn zero() -> Self {
        0.0
    }
    fn one() -> Self {
        1.0
    }
    fn half() -> Self {
        0.5
    }
    fn min(self, other: Self) -> Self {
        f64::min(self, other)
    }
    fn max(self, other: Self) -> Self {
        f64::max(self, other)
    }
    fn saturating_into_u16(self) -> u16 {
        self.max(0.).min(u16::MAX as f64) as u16
    }
    fn u16_max() -> Self {
        u16::MAX as f64
    }
}

#[derive(Debug, Default)]
pub struct FloatExtents<T> {
    scale: [T; 2],
    offset: [T; 2],
}

impl<T: Float> MapToU16<T> for FloatExtents<T> {
    fn map(&self, pos: [T; 2]) -> [u16; 2] {
        [
            (pos[0] * self.scale[0] + self.offset[0]).saturating_into_u16(),
            (pos[1] * self.scale[1] + self.offset[1]).saturating_into_u16(),
        ]
    }
}

impl<T: Float> Epsilon for T {
    type Extents = FloatExtents<T>;

    fn is_coincident(self, a: &[T; 2], b: &[T; 2]) -> bool {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        // Use squared epsilon to avoid sqrt and maintain consistent margin
        let epsilon_sq = self * self;
        dx * dx + dy * dy < epsilon_sq
    }

    fn is_on_edge(self, p: &[T; 2], edge_start: &[T; 2], edge_end: &[T; 2]) -> bool {
        let edge_dx = edge_end[0] - edge_start[0];
        let edge_dy = edge_end[1] - edge_start[1];
        let edge_len_sq = edge_dx * edge_dx + edge_dy * edge_dy;

        // Check if point is collinear with edge using cross product
        let to_p_x = p[0] - edge_start[0];
        let to_p_y = p[1] - edge_start[1];
        let cross = edge_dy * to_p_x - edge_dx * to_p_y;

        // Use squared epsilon for consistent margin
        let epsilon_sq = self * self;
        if cross * cross >= epsilon_sq * edge_len_sq {
            return false;
        }

        // Check if point is within edge bounds using dot product
        let dot = to_p_x * edge_dx + to_p_y * edge_dy;

        // Allow small epsilon margin on both ends
        // Use squared epsilon for consistency
        dot >= -epsilon_sq && dot <= edge_len_sq + epsilon_sq
    }

    fn intersect(
        self,
        a_start: &[T; 2],
        a_end: &[T; 2],
        b_start: &[T; 2],
        b_end: &[T; 2],
    ) -> Option<[T; 2]> {
        let dx1 = a_end[0] - a_start[0];
        let dy1 = a_end[1] - a_start[1];
        let dx2 = b_end[0] - b_start[0];
        let dy2 = b_end[1] - b_start[1];

        let det = dx1 * dy2 - dy1 * dx2;

        // Check if parallel/coincident (need abs() which requires conversion)
        let det_abs = if det >= T::zero() { det } else { -det };
        if det_abs < self {
            return None;
        }

        let dx3 = b_start[0] - a_start[0];
        let dy3 = b_start[1] - a_start[1];

        // TODO(Eric): Use multiplication instead of division where possible
        let det_inv = T::one() / det;
        let t1 = (dx3 * dy2 - dy3 * dx2) * det_inv;
        let t2 = (dx3 * dy1 - dy3 * dx1) * det_inv;

        let neg_epsilon = -self;
        let one_plus_epsilon = T::one() + self;

        // Check if intersection is within both edge segments (with epsilon margin)
        if t1 >= neg_epsilon
            && t1 <= one_plus_epsilon
            && t2 >= neg_epsilon
            && t2 <= one_plus_epsilon
        {
            Some([a_start[0] + t1 * dx1, a_start[1] + t1 * dy1])
        } else {
            None
        }
    }

    fn merged_with(self, a: &[T; 2], b: &[T; 2]) -> [T; 2] {
        let half = T::half();
        [half * (a[0] + b[0]), half * (a[1] + b[1])]
    }

    fn extents(self, mut edges: impl Iterator<Item = ([Self; 2], [Self; 2])>) -> Self::Extents {
        let Some((min, max)) = edges.next() else {
            return Default::default();
        };
        let (min, max) = edges.fold((min, max), |(min, max), (start, end)| {
            (
                [
                    min[0].min(start[0].min(end[0]) - self),
                    min[1].min(start[1].min(end[1]) - self),
                ],
                [
                    max[0].max(start[0].max(end[0]) + self),
                    max[1].max(start[1].max(end[1]) + self),
                ],
            )
        });

        let scale = [
            T::u16_max() / (max[0] - min[0]),
            T::u16_max() / (max[1] - min[1]),
        ];
        let offset = [-min[0] * scale[0], -min[1] * scale[1]];

        FloatExtents { scale, offset }
    }

    fn edge_bbox(self, edge_start: &[T; 2], edge_end: &[T; 2], extents: &Self::Extents) -> Rect {
        let min = [
            edge_start[0].min(edge_end[0]) - self,
            edge_start[1].min(edge_end[1]) - self,
        ];
        let max = [
            edge_start[0].max(edge_end[0]) + self,
            edge_start[1].max(edge_end[1]) + self,
        ];

        let min = extents.map(min);
        let max = extents.map(max);

        Rect {
            min_x: min[0],
            max_x: max[0],
            min_y: min[1],
            max_y: max[1],
        }
    }
}

/// Type alias for [f32; 2] vertices
pub type VertexF32 = [f32; 2];

/// Type alias for [f64; 2] vertices
pub type VertexF64 = [f64; 2];

/// Default epsilon value for f32 vertex operations
pub const DEFAULT_EPSILON_F32: f32 = 1e-5;

/// Default epsilon value for f64 vertex operations
pub const DEFAULT_EPSILON_F64: f64 = 1e-10;

impl Vertex for VertexF32 {
    type Extents = <f32 as Epsilon>::Extents;

    fn is_coincident(&self, other: &Self) -> bool {
        DEFAULT_EPSILON_F32.is_coincident(self, other)
    }

    fn is_on_edge(&self, edge_start: &Self, edge_end: &Self) -> bool {
        DEFAULT_EPSILON_F32.is_on_edge(self, edge_start, edge_end)
    }

    fn from_intersection(
        a_start: &Self,
        a_end: &Self,
        b_start: &Self,
        b_end: &Self,
    ) -> Option<Self> {
        DEFAULT_EPSILON_F32.intersect(a_start, a_end, b_start, b_end)
    }

    fn merged_with(&self, other: &Self) -> Self {
        DEFAULT_EPSILON_F32.merged_with(self, other)
    }

    fn extents<'a>(edges: impl Iterator<Item = (&'a Self, &'a Self)>) -> Self::Extents
    where
        Self: 'a,
    {
        DEFAULT_EPSILON_F32.extents(edges.map(|(&a, &b)| (a, b)))
    }

    fn edge_bbox(edge_start: &Self, edge_end: &Self, extents: &Self::Extents) -> Rect {
        DEFAULT_EPSILON_F32.edge_bbox(edge_start, edge_end, extents)
    }
}

impl Vertex for VertexF64 {
    type Extents = <f64 as Epsilon>::Extents;

    fn is_coincident(&self, other: &Self) -> bool {
        DEFAULT_EPSILON_F64.is_coincident(self, other)
    }

    fn is_on_edge(&self, edge_start: &Self, edge_end: &Self) -> bool {
        DEFAULT_EPSILON_F64.is_on_edge(self, edge_start, edge_end)
    }

    fn from_intersection(
        a_start: &Self,
        a_end: &Self,
        b_start: &Self,
        b_end: &Self,
    ) -> Option<Self> {
        DEFAULT_EPSILON_F64.intersect(a_start, a_end, b_start, b_end)
    }

    fn merged_with(&self, other: &Self) -> Self {
        DEFAULT_EPSILON_F64.merged_with(self, other)
    }

    fn extents<'a>(edges: impl Iterator<Item = (&'a Self, &'a Self)>) -> Self::Extents
    where
        Self: 'a,
    {
        DEFAULT_EPSILON_F64.extents(edges.map(|(&a, &b)| (a, b)))
    }

    fn edge_bbox(edge_start: &Self, edge_end: &Self, extents: &Self::Extents) -> Rect {
        DEFAULT_EPSILON_F64.edge_bbox(edge_start, edge_end, extents)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Tests using Vertex trait methods

    #[test]
    fn test_is_coincident_same_point() {
        let a = [0.5_f32, 0.5];
        let b = [0.5, 0.5];
        assert!(a.is_coincident(&b));
    }

    #[test]
    fn test_is_coincident_within_epsilon() {
        let a = [0.5_f32, 0.5];
        let b = [0.500001, 0.500001];
        assert!(a.is_coincident(&b));
    }

    #[test]
    fn test_is_coincident_outside_epsilon() {
        let a = [0.5_f32, 0.5];
        let b = [0.51, 0.5];
        assert!(!a.is_coincident(&b));
    }

    #[test]
    fn test_is_coincident_far_apart() {
        let a = [0.0_f32, 0.0];
        let b = [1.0, 1.0];
        assert!(!a.is_coincident(&b));
    }

    #[test]
    fn test_is_on_edge_midpoint() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let mid = [0.5, 0.0];
        assert!(mid.is_on_edge(&start, &end));
    }

    #[test]
    fn test_is_on_edge_at_start() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.0, 0.0];
        assert!(p.is_on_edge(&start, &end));
    }

    #[test]
    fn test_is_on_edge_at_end() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [1.0, 0.0];
        assert!(p.is_on_edge(&start, &end));
    }

    #[test]
    fn test_is_on_edge_not_on_line() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.5, 0.1];
        assert!(!p.is_on_edge(&start, &end));
    }

    #[test]
    fn test_is_on_edge_beyond_end() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [1.5, 0.0];
        assert!(!p.is_on_edge(&start, &end));
    }

    #[test]
    fn test_is_on_edge_before_start() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [-0.5, 0.0];
        assert!(!p.is_on_edge(&start, &end));
    }

    #[test]
    fn test_is_on_edge_diagonal() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 1.0];
        let mid = [0.5, 0.5];
        assert!(mid.is_on_edge(&start, &end));
    }

    #[test]
    fn test_is_on_edge_within_epsilon_margin() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 0.0];
        let p = [0.5, 0.000001]; // Very close to the line
        assert!(p.is_on_edge(&start, &end));
    }

    #[test]
    fn test_from_intersection_crossing_lines() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 1.0];
        let b_start = [0.0, 1.0];
        let b_end = [1.0, 0.0];

        let result = <[f32; 2]>::from_intersection(&a_start, &a_end, &b_start, &b_end);
        assert!(result.unwrap().is_coincident(&[0.5, 0.5]));
    }

    #[test]
    fn test_from_intersection_t_junction() {
        let a_start = [0.0_f32, 0.5];
        let a_end = [1.0, 0.5];
        let b_start = [0.5, 0.0];
        let b_end = [0.5, 1.0];

        let result = <[f32; 2]>::from_intersection(&a_start, &a_end, &b_start, &b_end);
        assert!(result.unwrap().is_coincident(&[0.5, 0.5]));
    }

    #[test]
    fn test_from_intersection_parallel() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 0.0];
        let b_start = [0.0, 1.0];
        let b_end = [1.0, 1.0];

        let result = <[f32; 2]>::from_intersection(&a_start, &a_end, &b_start, &b_end);
        assert!(result.is_none());
    }

    #[test]
    fn test_from_intersection_no_overlap() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.6, 0.4];
        let b_end = [1.0, 0.0];

        let result = <[f32; 2]>::from_intersection(&a_start, &a_end, &b_start, &b_end);
        assert!(result.is_none());
    }

    #[test]
    fn test_from_intersection_collinear_edges() {
        // Two collinear edges that share an endpoint
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.5, 0.5];
        let b_end = [1.0, 1.0];

        let result = <[f32; 2]>::from_intersection(&a_start, &a_end, &b_start, &b_end);
        // Collinear edges return None (det is too small)
        assert!(result.is_none());
    }

    #[test]
    fn test_from_intersection_crossing_at_endpoint() {
        // Two edges that share an endpoint but cross
        let a_start = [0.0_f32, 0.0];
        let a_end = [0.5, 0.5];
        let b_start = [0.5, 0.5];
        let b_end = [1.0, 0.0];

        let result = <[f32; 2]>::from_intersection(&a_start, &a_end, &b_start, &b_end);
        assert!(result.unwrap().is_coincident(&[0.5, 0.5]));
    }

    #[test]
    fn test_from_intersection_almost_parallel() {
        let a_start = [0.0_f32, 0.0];
        let a_end = [1.0, 0.0];
        let b_start = [0.0, 0.0];
        let b_end = [1.0, 0.000001]; // Almost parallel

        let result = <[f32; 2]>::from_intersection(&a_start, &a_end, &b_start, &b_end);
        // Should be None because det is too small
        assert!(result.is_none());
    }

    #[test]
    fn test_merge_simple() {
        let a = [0.0_f32, 0.0];
        let b = [1.0, 1.0];
        let merged = a.merged_with(&b);
        assert!(merged.is_coincident(&[0.5, 0.5]));
    }

    #[test]
    fn test_merge_same_point() {
        let a = [0.5_f32, 0.5];
        let b = [0.5, 0.5];
        let merged = a.merged_with(&b);
        assert!(merged.is_coincident(&[0.5, 0.5]));
    }

    #[test]
    fn test_merge_negative_coords() {
        let a = [-1.0_f32, -1.0];
        let b = [1.0, 1.0];
        let merged = a.merged_with(&b);
        assert!(merged.is_coincident(&[0.0, 0.0]));
    }

    #[test]
    fn test_edge_bbox_horizontal() {
        let start = [0.2_f32, 0.5];
        let end = [0.8, 0.5];
        let extents = <[f32; 2]>::extents([(&[0., 0.], &[1., 1.])].into_iter());
        let bbox = <[f32; 2]>::edge_bbox(&start, &end, &extents);

        assert!(bbox.min_x < bbox.max_x);
        assert!(bbox.min_y <= bbox.max_y);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_edge_bbox_vertical() {
        let start = [0.5_f32, 0.2];
        let end = [0.5, 0.8];
        let extents = <[f32; 2]>::extents([(&[0., 0.], &[1., 1.])].into_iter());
        let bbox = <[f32; 2]>::edge_bbox(&start, &end, &extents);

        assert!(bbox.min_x <= bbox.max_x);
        assert!(bbox.min_y < bbox.max_y);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_edge_bbox_diagonal() {
        let start = [0.0_f32, 0.0];
        let end = [1.0, 1.0];
        let extents = <[f32; 2]>::extents([(&[0., 0.], &[1., 1.])].into_iter());
        let bbox = <[f32; 2]>::edge_bbox(&start, &end, &extents);

        assert!(bbox.min_x < bbox.max_x);
        assert!(bbox.min_y < bbox.max_y);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_edge_bbox_point_edge() {
        let start = [0.5_f32, 0.5];
        let end = [0.5, 0.5];
        let extents = <[f32; 2]>::extents([(&[0., 0.], &[1., 1.])].into_iter());
        let bbox = <[f32; 2]>::edge_bbox(&start, &end, &extents);

        assert!(bbox.min_x <= bbox.max_x);
        assert!(bbox.min_y <= bbox.max_y);
        assert!(bbox.overlaps(&bbox));
    }

    #[test]
    fn test_bboxes_overlap() {
        // Edge 1 is horizontal
        let start1 = [0.5_f32, 0.];
        let end1 = [0.5, 1.];

        // Edge 2 is vertical
        let start2 = [0.5_f32, 0.5];
        let end2 = [1., 0.5];

        let extents = <[f32; 2]>::extents([(&[0., 0.], &[1., 1.])].into_iter());
        let bbox1 = <[f32; 2]>::edge_bbox(&start1, &end1, &extents);
        let bbox2 = <[f32; 2]>::edge_bbox(&start2, &end2, &extents);

        assert!(bbox1.overlaps(&bbox2));
    }

    #[test]
    fn test_bboxes_dont_overlap() {
        let start1 = [0.1_f32, 0.1];
        let end1 = [0.2, 0.2];

        let start2 = [0.2_f32, 0.4];
        let end2 = [0.1, 0.5];

        let extents = <[f32; 2]>::extents([(&[0., 0.], &[1., 1.])].into_iter());
        let bbox1 = <[f32; 2]>::edge_bbox(&start1, &end1, &extents);
        let bbox2 = <[f32; 2]>::edge_bbox(&start2, &end2, &extents);

        assert!(!bbox1.overlaps(&bbox2));
    }
}

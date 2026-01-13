use crate::rtree::Rect;
use crate::triangle_kernel::TriangleKernel;
use std::cmp::Ordering;
use std::fmt::Debug;
pub mod polyline;

// Re-export sweep-line types for convenience
pub use crate::sweep_line::{SweepLineChain, SweepLineEvent, SweepLineEventType, SweepLineSegment};

pub trait Kernel: Sized {
    type Vertex: Copy + Ord + Debug;
    type Edge: Edge;
    type Extents;
    type Intersection;
    type SweepLineEdgePortion: Copy + PartialEq + Debug;
    type SweepLineEventPoint: Copy + PartialEq + Debug;
    type TriangleKernel: TriangleKernel;

    /// Check if two vertices are coincident (at the same location)
    fn vertices_coincident(&self, a: Self::Vertex, b: Self::Vertex) -> bool;

    /// Check if two edges are coincident (fully overlap)
    fn edges_coincident(&self, a: Self::Edge, b: Self::Edge) -> bool;

    /// Check if this vertex lies on the edge from edge_start to edge_end
    fn vertex_on_edge(&self, vertex: Self::Vertex, edge: Self::Edge) -> bool;

    /// See if two edges intersect
    fn intersection(&self, a: Self::Edge, b: Self::Edge) -> Option<Self::Intersection>;

    /// Merge two vertices, returning the merged result
    fn merged_vertex(&mut self, a: Self::Vertex, b: Self::Vertex) -> Self::Vertex;

    /// Merge two edges, returning the merged result.
    /// These edges must already share all of their vertices
    /// The two edges that are returned must either be equal or cancel.
    fn merged_edges(&mut self, a: Self::Edge, b: Self::Edge) -> (Self::Edge, Self::Edge);

    /// Creates and returns the vertex for an intersection
    fn intersection_vertex(&mut self, intersection: Self::Intersection) -> Self::Vertex;

    /// Generate an "extents" object from a list of edges,
    /// which is used when calculating the edge_bbox
    fn extents(&self, edges: impl Iterator<Item = Self::Edge>) -> Self::Extents;

    /// Compute the axis-aligned bounding box of an edge
    fn edge_bbox(&self, edge: Self::Edge, extents: &Self::Extents) -> Rect;

    /// Compare the angular order of two vectors from self to a and self to b.
    /// Returns the sign of the cross product: (a - self) x (b - self).
    /// Greater = b is counterclockwise from a (positive cross product)
    /// Equal = collinear (zero cross product)
    /// Less = b is clockwise from a (negative cross product)
    fn sin_cmp(&self, common: Self::Vertex, a: Self::Vertex, b: Self::Vertex) -> Ordering;

    // Returns the vertices which are endpoints for this edge
    fn vertices_for_edge(&self, edge: Self::Edge) -> Few<Self::Vertex>;

    // Replaces instances of old_v with new_v in edge
    // This returns None if the resultant edge is a reflex edge with zero area.
    fn replace_vertex_in_edge(
        &self,
        edge: Self::Edge,
        old_v: Self::Vertex,
        new_v: Self::Vertex,
    ) -> Option<Self::Edge>;

    // Splits an edge at the given vertex, returning two new edges
    // `vertex` must not be equal to or coincident with an endpoint already on `edge`
    fn split_edge(&self, edge: Self::Edge, vertex: Self::Vertex) -> (Self::Edge, Self::Edge);

    fn sweep_line_events_for_edge(
        &self,
        edge: Self::Edge,
    ) -> impl Iterator<Item = SweepLineEvent<Self>>;

    fn sweep_line_event_cmp_bottom_up(
        &self,
        a: &SweepLineEvent<Self>,
        b: &SweepLineEvent<Self>,
    ) -> Ordering;

    fn sweep_line_event_point(&self, event: &SweepLineEvent<Self>) -> Self::SweepLineEventPoint;

    fn sweep_line_segment_cmp(
        &self,
        segment: &SweepLineSegment<Self>,
        event_point: Self::SweepLineEventPoint,
    ) -> Ordering;

    // Triangulation methods
    fn sweep_line_event_point_to_triangle_vertex(
        &self,
        triangle_kernel: &mut Self::TriangleKernel,
        event_point: Self::SweepLineEventPoint,
    ) -> <Self::TriangleKernel as TriangleKernel>::Vertex;

    fn sweep_line_edge_segment_to_triangle_vertices(
        &self,
        triangle_kernel: &mut Self::TriangleKernel,
        segment: &SweepLineSegment<Self>,
    ) -> impl Iterator<Item = <Self::TriangleKernel as TriangleKernel>::Vertex>;

    fn sweep_line_event_cmp_clockwise(
        &self,
        a: &SweepLineEvent<Self>,
        b: &SweepLineEvent<Self>,
    ) -> Ordering;
}

pub trait Edge: Copy + Ord + Debug {
    const MIN: Self;
    const MAX: Self;

    fn reversed(self) -> Self;
}

pub enum Few<T> {
    Zero,
    One(T),
    Two(T, T),
}

impl<T: Copy> Iterator for Few<T> {
    type Item = T;

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

impl Edge for (u32, u32) {
    const MIN: Self = (u32::MIN, u32::MIN);
    const MAX: Self = (u32::MAX, u32::MAX);

    fn reversed(self) -> Self {
        (self.1, self.0)
    }
}

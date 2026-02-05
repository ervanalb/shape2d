use crate::rtree::Rect;
use crate::triangle_kernel::TriangleKernel;
use std::cmp::Ordering;
use std::fmt::Debug;
pub mod line;

use crate::sweep_line::{SweepLineEvent, SweepLineSegment};

pub trait Kernel: Sized {
    type Vertex: Copy + Ord + Debug;
    type Edge: Edge;
    type Extents;
    type MergePoint;
    type MergeCurve;
    type SplitPoint;
    type IntersectionPoint;
    type SweepLineEdgePortion: Copy + PartialEq + Debug;
    type SweepLineEventPoint: Copy + PartialEq + Debug;
    type TriangleKernel: TriangleKernel;
    type CapStyle: Copy + Debug;
    type OffsetAmount: Copy + Debug;

    /// Check if two vertices are coincident (at the same location)
    fn vertices_coincident(&self, a: Self::Vertex, b: Self::Vertex) -> Option<Self::MergePoint>;

    /// Check if two edges are coincident (fully overlap)
    fn edges_coincident(&self, a: Self::Edge, b: Self::Edge) -> EdgeInteraction<Self::MergeCurve>;

    /// Check if this vertex lies on the edge from edge_start to edge_end,
    /// and if so, return the point on the edge nearest to the vertex
    fn split(&self, vertex: Self::Vertex, edge: Self::Edge) -> Option<Self::SplitPoint>;

    /// See if two edges intersect
    fn intersection(&self, a: Self::Edge, b: Self::Edge) -> Option<Self::IntersectionPoint>;

    /// Merge two vertices, returning the merged result
    fn merge_vertices(&mut self, pt: Self::MergePoint) -> Self::Vertex;

    /// Merge two edges, returning the merged result.
    /// These edges must already share all of their vertices
    fn merge_edges(&mut self, curve: Self::MergeCurve) -> Self::Edge;

    /// Creates a new vertex from a split edge
    fn new_split_vertex(&mut self, pt: Self::SplitPoint) -> Self::Vertex;

    /// Creates a new vertex from an intersection
    fn new_intersection_vertex(&mut self, pt: Self::IntersectionPoint) -> Self::Vertex;

    /// Generate an "extents" object from a list of edges,
    /// which is used when calculating the edge_bbox
    fn extents(&self, edges: impl Iterator<Item = Self::Edge>) -> Self::Extents;

    /// Compute the axis-aligned bounding box of an edge
    fn edge_bbox(&self, edge: Self::Edge, extents: &Self::Extents) -> Rect;

    /// Compare the angular order of two vectors from common to a and common to b.
    /// Less = a comes before b in counterclockwise ordering (a x b yields a positive cross product)
    /// Equal = collinear (zero cross product)
    /// Greater = a comes after b in counterclockwise ordering (a x b yields a negative cross product)
    fn sin_cmp(&self, common: Self::Vertex, a: Self::Vertex, b: Self::Vertex) -> Ordering;

    // Returns the vertex that this edge starts at, and the vertex that this edge ends at.
    // Returns None if this edge has no vertices (e.g. a circle.)
    // These two vertices may be the same.
    fn vertices_for_edge(&self, edge: Self::Edge) -> Option<(Self::Vertex, Self::Vertex)>;

    // Replaces instances of old_v with new_v in edge
    // This returns None if the resultant edge is a reflex edge with zero area.
    fn replace_vertex_in_edge(
        &mut self,
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

    fn sweep_line_event_cmp(&self, a: &SweepLineEvent<Self>, b: &SweepLineEvent<Self>) -> Ordering;

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

    // Offset methods
    fn vertex_event_cmp(&self, a: &VertexEvent<Self>, b: &VertexEvent<Self>) -> Ordering;

    fn offset_edge_loops(
        &mut self,
        edges: &[(u32, Self::Edge)],
        offset: Self::OffsetAmount,
        cap_style: Self::CapStyle,
    ) -> Vec<Self::Edge>;
}

pub trait Edge: Copy + Ord + Debug {
    const MIN: Self;
    const MAX: Self;

    fn reversed(self) -> Self;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum EdgeSide {
    Tail,
    Head,
}

impl EdgeSide {
    pub fn other(self) -> Self {
        match self {
            EdgeSide::Tail => EdgeSide::Head,
            EdgeSide::Head => EdgeSide::Tail,
        }
    }

    pub fn select<T>(self, (a, b): (T, T)) -> T {
        match self {
            EdgeSide::Tail => a,
            EdgeSide::Head => b,
        }
    }
}

#[derive(Clone)]
pub struct VertexEvent<K: Kernel> {
    pub(crate) event_type: EdgeSide,
    pub(crate) edge: K::Edge,
}

impl<K: Kernel> std::fmt::Debug for VertexEvent<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("VertexEvent")
            .field("event_type", &self.event_type)
            .field("edge", &self.edge)
            .finish()
    }
}

impl Edge for (u32, u32) {
    const MIN: Self = (u32::MIN, u32::MIN);
    const MAX: Self = (u32::MAX, u32::MAX);

    fn reversed(self) -> Self {
        (self.1, self.0)
    }
}

pub enum EdgeInteraction<T> {
    None,
    Merge(T),  // Curves are in same direction
    Cancel(T), // Curves are in opposite directions
}

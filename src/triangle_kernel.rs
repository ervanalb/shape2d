use std::cmp::Ordering;
use std::fmt::Debug;

/// Triangle kernel trait for managing vertices and triangles during triangulation
pub trait TriangleKernel {
    /// The vertex handle type (e.g., u32 for indexed vertices)
    type Vertex: Copy + Ord + Debug;

    /// The triangle type
    type Triangle: Copy + Ord + Debug;

    /// Create a triangle from three vertex handles
    fn push_triangle(
        &mut self,
        v0: Self::Vertex,
        v1: Self::Vertex,
        v2: Self::Vertex,
    ) -> Self::Triangle;

    /// Compare vertices in sweep-line order (left-to-right, bottom-to-top)
    fn sweep_line_cmp(&self, a: Self::Vertex, b: Self::Vertex) -> Ordering;

    /// Compare the angular order of vectors from common to a and common to b
    /// Returns the sign of the cross product: (a - common) x (b - common)
    /// Greater = b is counterclockwise from a (positive cross product)
    /// Equal = collinear (zero cross product)
    /// Less = b is clockwise from a (negative cross product)
    fn sin_cmp(&self, common: Self::Vertex, a: Self::Vertex, b: Self::Vertex) -> Ordering;
}

/// Concrete implementation of TriangleKernel for f32 coordinates
pub struct F32TriangleKernel {
    vertices: Vec<[f32; 2]>,
}

impl F32TriangleKernel {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
        }
    }

    pub fn push_vertex(&mut self, data: [f32; 2]) -> u32 {
        let index = self.vertices.len() as u32;
        self.vertices.push(data);
        index
    }

    pub fn v(&self, i: u32) -> [f32; 2] {
        self.vertices[i as usize]
    }
}

impl TriangleKernel for F32TriangleKernel {
    type Vertex = u32;
    type Triangle = (u32, u32, u32);

    fn push_triangle(
        &mut self,
        v0: Self::Vertex,
        v1: Self::Vertex,
        v2: Self::Vertex,
    ) -> Self::Triangle {
        (v0, v1, v2)
    }

    fn sweep_line_cmp(&self, a: Self::Vertex, b: Self::Vertex) -> Ordering {
        let a_v = self.v(a);
        let b_v = self.v(b);
        crate::kernel::polyline::sweep_line_cmp_f32(a_v, b_v)
    }

    fn sin_cmp(&self, common: Self::Vertex, a: Self::Vertex, b: Self::Vertex) -> Ordering {
        let common_v = self.v(common);
        let a_v = self.v(a);
        let b_v = self.v(b);
        crate::kernel::polyline::sin_cmp_f32(common_v, a_v, b_v)
    }
}

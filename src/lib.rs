mod clean;
mod hilbert;
mod rect;
mod rtree;
//mod clip;
pub mod geometry;

pub use clean::{DirtyFlag, clean, partial_clean};
pub use rect::Rect;
//pub use clip::clip;
//pub use vertex::{DEFAULT_EPSILON_F32, DEFAULT_EPSILON_F64, Epsilon, Vertex, VertexF32, VertexF64};

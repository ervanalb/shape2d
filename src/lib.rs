mod hilbert;
mod rect;
mod rtree;
mod clean;
mod clip;
pub mod vertex;

pub use rect::Rect;
pub use vertex::{DEFAULT_EPSILON_F32, DEFAULT_EPSILON_F64, Epsilon, Vertex, VertexF32, VertexF64};
pub use clean::{clean, partial_clean, DirtyFlag};
pub use clip::clip;

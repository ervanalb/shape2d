use crate::rtree::RTree;
use crate::{Rect, Vertex};
use std::collections::{BTreeMap, BTreeSet};

type EdgeIndex = u32;
type VertexIndex = u32;

/// Flag indicating whether an edge is known to be clean or might be dirty
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirtyFlag {
    /// Edge is known to not intersect or be coincident with other clean edges
    Clean,
    /// Edge may intersect or be coincident with other edges
    Dirty,
}

/// Information about an edge
#[derive(Debug, Clone)]
struct EdgeInfo {
    start: VertexIndex,
    end: VertexIndex,
    hilbert_value: u32,
    bbox: Rect,
}

/// Actions that can be taken during the cleaning process
#[derive(Debug)]
enum Action<V: Vertex> {
    /// Merge two coincident vertices
    MergeVertices { v1: VertexIndex, v2: VertexIndex },
    /// Split an edge at a vertex
    SplitEdge {
        edge: EdgeIndex,
        split_vertex: VertexIndex,
    },
    /// Split both edges at their intersection point
    SplitBothEdges {
        edge_a: EdgeIndex,
        edge_b: EdgeIndex,
        intersection: V,
    },
}

/// Main spatial index structure
struct SpatialIndex<'a, V: Vertex> {
    vertices: &'a mut Vec<V>,
    edge_info: Vec<EdgeInfo>,
    extents: V::Extents,
    hilbert_to_edges: BTreeSet<(u32, EdgeIndex)>,
    vertex_to_edges: BTreeSet<(VertexIndex, EdgeIndex)>,
    dirty_set: BTreeSet<EdgeIndex>,
    r_tree: RTree,
}

impl<'a, V: Vertex> SpatialIndex<'a, V> {
    /// Build the spatial index from vertices and edges
    fn new(
        vertices: &'a mut Vec<V>,
        edges: impl Iterator<Item = (VertexIndex, VertexIndex, DirtyFlag)>,
    ) -> Option<Self> {
        let mut edge_info = Vec::new();
        let mut hilbert_to_edges = BTreeSet::new();
        let mut vertex_to_edges = BTreeSet::new();
        let mut dirty_set = BTreeSet::new();
        let mut hilbert_to_rect: BTreeMap<u32, Rect> = BTreeMap::new();

        // Loop over edges to fill out our edge_info struct partially
        for (edge_idx, (start, end, dirty_flag)) in edges.enumerate() {
            // Store edge info
            edge_info.push(EdgeInfo {
                start,
                end,
                hilbert_value: 0,
                bbox: Default::default(),
            });

            // Add to dirty set if dirty
            if dirty_flag == DirtyFlag::Dirty {
                dirty_set.insert(edge_idx as u32);
            }
        }

        // Loop over our edges to find the extents
        let extents = V::extents(edge_info.iter().map(|&EdgeInfo { start, end, .. }| {
            (&vertices[start as usize], &vertices[end as usize])
        }))?;

        // Now that we have the extents,
        // fill in the rest of the missing data
        for (edge_idx, edge_info) in edge_info.iter_mut().enumerate() {
            let edge_idx = edge_idx as u32;

            // Calculate bounding box
            let bbox = V::edge_bbox(
                &vertices[edge_info.start as usize],
                &vertices[edge_info.end as usize],
                &extents,
            );

            // Calculate hilbert value
            let hilbert_value = bbox.hilbert_value();

            // Add to hilbert_to_edges
            hilbert_to_edges.insert((hilbert_value, edge_idx));

            // Add to vertex_to_edges
            vertex_to_edges.insert((edge_info.start, edge_idx));
            vertex_to_edges.insert((edge_info.end, edge_idx));

            // Update hilbert_to_rect mapping
            hilbert_to_rect
                .entry(hilbert_value)
                .and_modify(|existing_rect| *existing_rect = existing_rect.combine(&bbox))
                .or_insert(bbox);
        }

        // Build the R-tree
        let r_tree = RTree::build(hilbert_to_rect);

        Some(Self {
            vertices,
            edge_info,
            extents,
            hilbert_to_edges,
            vertex_to_edges,
            dirty_set,
            r_tree,
        })
    }

    /// Insert a new edge into the spatial index
    fn insert(&mut self, start: VertexIndex, end: VertexIndex, hilbert_value: u32) {
        let edge_idx = self.edge_info.len() as EdgeIndex;

        let bbox = V::edge_bbox(
            &self.vertices[start as usize],
            &self.vertices[end as usize],
            &self.extents,
        );

        self.edge_info.push(EdgeInfo {
            start,
            end,
            hilbert_value,
            bbox,
        });

        self.hilbert_to_edges.insert((hilbert_value, edge_idx));
        self.vertex_to_edges.insert((start, edge_idx));
        self.vertex_to_edges.insert((end, edge_idx));
        self.dirty_set.insert(edge_idx);
        self.r_tree.enlarge(hilbert_value, bbox);
    }

    /// Edit an edge's endpoints
    fn edit(&mut self, edge_idx: EdgeIndex, new_start: VertexIndex, new_end: VertexIndex) {
        let info = &self.edge_info[edge_idx as usize];
        let old_start = info.start;
        let old_end = info.end;
        let hilbert_value = info.hilbert_value;

        // Calculate new bounding box
        let new_bbox = V::edge_bbox(
            &self.vertices[new_start as usize],
            &self.vertices[new_end as usize],
            &self.extents,
        );

        // Update edge info
        self.edge_info[edge_idx as usize].start = new_start;
        self.edge_info[edge_idx as usize].end = new_end;
        self.edge_info[edge_idx as usize].bbox = new_bbox;

        // Update vertex_to_edges
        if old_start != new_start {
            self.vertex_to_edges.remove(&(old_start, edge_idx));
            self.vertex_to_edges.insert((new_start, edge_idx));
        }
        if old_end != new_end {
            self.vertex_to_edges.remove(&(old_end, edge_idx));
            self.vertex_to_edges.insert((new_end, edge_idx));
        }

        // Mark as dirty
        self.dirty_set.insert(edge_idx);

        // Enlarge R-tree
        self.r_tree.enlarge(hilbert_value, new_bbox);
    }

    /// Get all edges that reference a given vertex
    fn edges_for_vertex(&self, vertex_idx: VertexIndex) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.vertex_to_edges
            .range((vertex_idx, 0)..=(vertex_idx, u32::MAX))
            .map(|(_, edge_idx)| *edge_idx)
    }

    /// Get all edges with a given hilbert value
    fn edges_for_hilbert(&self, hilbert_value: u32) -> impl Iterator<Item = EdgeIndex> + '_ {
        self.hilbert_to_edges
            .range((hilbert_value, 0)..=(hilbert_value, u32::MAX))
            .map(|(_, edge_idx)| *edge_idx)
    }

    /// Run the cleaning algorithm
    fn clean(&mut self) {
        while let Some(&dirty_edge_idx) = self.dirty_set.iter().next() {
            // Remove from dirty set first
            self.dirty_set.remove(&dirty_edge_idx);

            // Get edge info
            let dirty_info = &self.edge_info[dirty_edge_idx as usize];

            let mut action: Option<Action<V>> = None;

            // Iterate over all potentially intersecting edges
            'itest: for hilbert_value in self.r_tree.search(&dirty_info.bbox) {
                for candidate_idx in self.edges_for_hilbert(hilbert_value) {
                    let candidate_info = self.edge_info[candidate_idx as usize].clone();

                    // Get vertices
                    let dirty_vertices = (dirty_info.start, dirty_info.end);
                    let candidate_vertices = (candidate_info.start, candidate_info.end);

                    // Test 1: Check for coincident vertices
                    // Check all 6 pairs: (a.0, a.1],(a.0, b.0],(a.0, b.1],(a.1, b.0],(a.1, b.1],(b.0, b.1)
                    let pairs = [
                        (dirty_vertices.0, dirty_vertices.1),
                        (dirty_vertices.0, candidate_vertices.0),
                        (dirty_vertices.0, candidate_vertices.1),
                        (dirty_vertices.1, candidate_vertices.0),
                        (dirty_vertices.1, candidate_vertices.1),
                        (candidate_vertices.0, candidate_vertices.1),
                    ];

                    for (v1, v2) in pairs {
                        if v1 != v2
                            && self.vertices[v1 as usize].is_coincident(&self.vertices[v2 as usize])
                        {
                            action = Some(Action::MergeVertices { v1, v2 });
                            break 'itest;
                        }
                    }

                    // Test 2: Check for vertex on edge
                    // Check if either endpoint of dirty_edge lies on candidate_edge
                    for &vertex_idx in &[dirty_info.start, dirty_info.end] {
                        if vertex_idx != candidate_info.start
                            && vertex_idx != candidate_info.end
                            && self.vertices[vertex_idx as usize].is_on_edge(
                                &self.vertices[candidate_info.start as usize],
                                &self.vertices[candidate_info.end as usize],
                            )
                        {
                            action = Some(Action::SplitEdge {
                                edge: candidate_idx,
                                split_vertex: vertex_idx,
                            });
                            break 'itest;
                        }
                    }

                    // Check if either endpoint of candidate_edge lies on dirty_edge
                    for &vertex_idx in &[candidate_info.start, candidate_info.end] {
                        if vertex_idx != dirty_info.start
                            && vertex_idx != dirty_info.end
                            && self.vertices[vertex_idx as usize].is_on_edge(
                                &self.vertices[dirty_info.start as usize],
                                &self.vertices[dirty_info.end as usize],
                            )
                        {
                            action = Some(Action::SplitEdge {
                                edge: dirty_edge_idx,
                                split_vertex: vertex_idx,
                            });
                            break 'itest;
                        }
                    }

                    // Test 3: Check for edge intersection
                    // Skip if edges share a vertex (they can't intersect in the middle)
                    if dirty_info.start != candidate_info.start
                        && dirty_info.start != candidate_info.end
                        && dirty_info.end != candidate_info.start
                        && dirty_info.end != candidate_info.end
                        && let Some(intersection) = V::from_intersection(
                            &self.vertices[dirty_info.start as usize],
                            &self.vertices[dirty_info.end as usize],
                            &self.vertices[candidate_info.start as usize],
                            &self.vertices[candidate_info.end as usize],
                        )
                    {
                        action = Some(Action::SplitBothEdges {
                            edge_a: dirty_edge_idx,
                            edge_b: candidate_idx,
                            intersection,
                        });
                        break 'itest;
                    }
                }
            }

            // Handle the action if one was found
            if let Some(act) = action {
                match act {
                    Action::MergeVertices { v1, v2 } => {
                        // Add merged vertex as a new vertex
                        let merged_vertex =
                            self.vertices[v1 as usize].merged_with(&self.vertices[v2 as usize]);
                        let new_v = self.vertices.len() as VertexIndex;
                        self.vertices.push(merged_vertex);

                        // Update all edges referencing v1 or v2
                        for v in [v1, v2] {
                            while let Some(edge_idx) = { self.edges_for_vertex(v).next() } {
                                let EdgeInfo {
                                    mut start, mut end, ..
                                } = self.edge_info[edge_idx as usize];
                                if start == v {
                                    start = new_v;
                                }
                                if end == v {
                                    end = new_v;
                                }
                                self.edit(edge_idx, start, end);
                            }
                        }
                    }
                    Action::SplitEdge { edge, split_vertex } => {
                        let EdgeInfo {
                            start,
                            end: old_end,
                            hilbert_value,
                            ..
                        } = self.edge_info[edge as usize];

                        // Shorten the original edge
                        self.edit(edge, start, split_vertex);

                        // Create new edge from split point to old endpoint
                        self.insert(split_vertex, old_end, hilbert_value);
                    }
                    Action::SplitBothEdges {
                        edge_a,
                        edge_b,
                        intersection,
                    } => {
                        // Add intersection point as a new vertex
                        let new_vertex_idx = self.vertices.len() as VertexIndex;
                        self.vertices.push(intersection);

                        // Get edge info before modifying edges
                        let EdgeInfo {
                            start: edge_a_start,
                            end: edge_a_end,
                            hilbert_value: edge_a_hilbert,
                            ..
                        } = self.edge_info[edge_a as usize];
                        let EdgeInfo {
                            start: edge_b_start,
                            end: edge_b_end,
                            hilbert_value: edge_b_hilbert,
                            ..
                        } = self.edge_info[edge_b as usize];

                        // Shorten both edges to intersection point
                        self.edit(edge_a, edge_a_start, new_vertex_idx);
                        self.edit(edge_b, edge_b_start, new_vertex_idx);

                        // Create new edges from intersection to old endpoints
                        self.insert(new_vertex_idx, edge_a_end, edge_a_hilbert);
                        self.insert(new_vertex_idx, edge_b_end, edge_b_hilbert);
                    }
                }
            }
        }
    }

    /// Extract the final clean edge list
    fn extract_edges(&self) -> Vec<(u32, u32)> {
        self.hilbert_to_edges
            .iter()
            .map(|(_, edge_idx)| {
                let info = &self.edge_info[*edge_idx as usize];
                (info.start, info.end)
            })
            .collect()
    }
}

/// Partial clean: clean edges with mixed dirty/clean flags
pub fn partial_clean<V: Vertex>(
    vertices: &mut Vec<V>,
    edges: impl Iterator<Item = (VertexIndex, VertexIndex, DirtyFlag)>,
) -> Vec<(u32, u32)> {
    let Some(mut spatial_index) = SpatialIndex::new(vertices, edges) else {
        // No input edges
        return vec![];
    };
    spatial_index.clean();
    spatial_index.extract_edges()
}

/// Full clean: clean all edges (assumes all dirty)
pub fn clean<V: Vertex>(
    vertices: &mut Vec<V>,
    edges: impl Iterator<Item = (VertexIndex, VertexIndex)>,
) -> Vec<(u32, u32)> {
    partial_clean(
        vertices,
        edges.map(|(start, end)| (start, end, DirtyFlag::Dirty)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_two_edges() {
        let mut vertices = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]];

        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut vertices, edges.into_iter());
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_intersecting_edges() {
        let mut vertices = vec![[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]];

        // Two edges that intersect
        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut vertices, edges.into_iter());
        // Should split into 4 edges
        assert_eq!(result.len(), 4);
        // Should create a new vertex at the intersection
        assert_eq!(vertices.len(), 5);
    }

    #[test]
    fn test_vertex_on_edge() {
        let mut vertices = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 0.0], // On the first edge
        ];

        let edges = [(0, 1)];

        let result = clean(&mut vertices, edges.into_iter());
        // Should not split since there's only one edge and the vertex isn't connected
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_coincident_vertices() {
        let mut vertices = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 0.0], // Coincident with vertex 0
            [1.0, 1.0],
        ];

        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut vertices, edges.into_iter());
        // Should merge vertices 0 and 2
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_square_polygon() {
        let mut vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];

        let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];

        let result = clean(&mut vertices, edges.into_iter());
        // Square edges don't intersect, should remain 4 edges
        assert_eq!(result.len(), 4);
        // No new vertices should be created
        assert_eq!(vertices.len(), 4);
    }

    #[test]
    fn test_t_junction() {
        let mut vertices = vec![[0.0, 0.5], [1.0, 0.5], [0.5, 0.0], [0.5, 1.0]];

        // T-junction: horizontal edge intersected by vertical edge
        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut vertices, edges.into_iter());
        // Should split into 4 edges meeting at center
        assert_eq!(result.len(), 4);
        // Should create intersection vertex at (0.5, 0.5)
        assert_eq!(vertices.len(), 5);
    }

    // Integration tests using [f32; 2] vertex type
    mod f32_array_tests {
        use super::*;

        #[test]
        fn test_f32_simple_square() {
            let mut vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];

            let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];

            let result = clean(&mut vertices, edges.into_iter());
            assert_eq!(result.len(), 4);
            assert_eq!(vertices.len(), 4);
        }

        #[test]
        fn test_f32_intersecting_diagonals() {
            let mut vertices = vec![[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]];

            let edges = [(0, 1), (2, 3)];

            let result = clean(&mut vertices, edges.into_iter());
            assert_eq!(result.len(), 4);
            assert_eq!(vertices.len(), 5);

            // Check that intersection point is at [0.5, 0.5]
            let intersection = vertices[4];
            assert!(intersection.is_coincident(&[0.5, 0.5]));
        }

        #[test]
        fn test_f32_coincident_vertices() {
            let mut vertices = vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [0.000001, 0.000001], // Very close to [0.0, 0.0]
                [1.0, 1.0],
            ];

            let edges = [(0, 1), (2, 3)];

            let result = clean(&mut vertices, edges.into_iter());
            assert_eq!(result.len(), 2);
        }

        #[test]
        fn test_f32_t_junction() {
            let mut vertices = vec![[0.0, 0.5], [1.0, 0.5], [0.5, 0.0], [0.5, 1.0]];

            let edges = [(0, 1), (2, 3)];

            let result = clean(&mut vertices, edges.into_iter());
            assert_eq!(result.len(), 4);
            assert_eq!(vertices.len(), 5);

            // Check intersection point
            let intersection = vertices[4];
            assert!(intersection.is_coincident(&[0.5, 0.5]));
        }

        #[test]
        fn test_f32_multiple_intersections() {
            let mut vertices = vec![
                // Horizontal line
                [0.0, 0.5],
                [1.0, 0.5],
                // Vertical line
                [0.5, 0.0],
                [0.5, 1.0],
                // Diagonal line
                [0.0, 0.0],
                [1.0, 1.0],
            ];

            let edges = [(0, 1), (2, 3), (4, 5)];

            let result = clean(&mut vertices, edges.into_iter());
            // Each of the 3 edges should be split at their intersections
            // Horizontal and vertical intersect at [0.5, 0.5]
            // Diagonal intersects both at [0.5, 0.5]
            // So all three meet at one point, creating 6 edges
            assert_eq!(result.len(), 6);
        }

        #[test]
        fn test_f32_complex_polygon() {
            // Star-like pattern
            let mut vertices = vec![
                [0.5, 0.0], // Top
                [0.7, 0.5], // Right
                [0.5, 1.0], // Bottom
                [0.3, 0.5], // Left
                [0.5, 0.5], // Center
            ];

            let edges = [
                (0, 2), // Vertical through center
                (1, 3), // Horizontal through center
            ];

            let result = clean(&mut vertices, edges.into_iter());
            // Both edges pass through center, should be split there
            assert_eq!(result.len(), 4);
        }
    }
}

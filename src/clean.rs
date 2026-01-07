use crate::rtree::RTree;
use crate::{Geometry, Rect};
use std::collections::{BTreeMap, BTreeSet};

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
struct EdgeInfo<G: Geometry> {
    edge: G::Edge,
    hilbert_value: u32,
    bbox: Rect,
}

/// Actions that can be taken during the cleaning process
#[derive(Debug)]
enum Action<G: Geometry> {
    /// Merge two coincident vertices
    MergeVertices { v1: G::Vertex, v2: G::Vertex },
    /// Split an edge at a vertex
    SplitEdge { edge: u32, split_vertex: G::Vertex },
    /// Split both edges at their intersection point
    SplitBothEdges {
        edge_a: u32,
        edge_b: u32,
        intersection: G::Intersection,
    },
}

/// Main spatial index structure
struct SpatialIndex<'a, G: Geometry> {
    geometry: &'a mut G,
    edge_info: Vec<EdgeInfo<G>>,
    extents: G::Extents,
    hilbert_to_edges: BTreeSet<(u32, u32)>,
    vertex_to_edges: BTreeSet<(G::Vertex, u32)>,
    dirty_set: BTreeSet<u32>,
    r_tree: RTree,
}

impl<'a, G: Geometry> SpatialIndex<'a, G> {
    /// Build the spatial index from vertices and edges
    fn new(geometry: &'a mut G, edges: impl Iterator<Item = (G::Edge, DirtyFlag)>) -> Self {
        let mut edge_info = Vec::new();
        let mut hilbert_to_edges = BTreeSet::new();
        let mut vertex_to_edges = BTreeSet::new();
        let mut dirty_set = BTreeSet::new();
        let mut hilbert_to_rect: BTreeMap<u32, Rect> = BTreeMap::new();

        // Loop over edges to fill out our edge_info struct partially
        for (edge_idx, (edge, dirty_flag)) in edges.enumerate() {
            // Store edge info
            edge_info.push(EdgeInfo {
                edge,
                hilbert_value: 0,
                bbox: Default::default(),
            });

            // Add to dirty set if dirty
            if dirty_flag == DirtyFlag::Dirty {
                dirty_set.insert(edge_idx as u32);
            }
        }

        // Loop over our edges to find the extents
        let extents = geometry.extents(edge_info.iter().map(|ei| ei.edge));

        // Now that we have the extents,
        // fill in the rest of the missing data
        for (edge_idx, edge_info) in edge_info.iter_mut().enumerate() {
            let edge_idx = edge_idx as u32;

            // Calculate bounding box and hilbert value
            let bbox = geometry.edge_bbox(edge_info.edge, &extents);
            let hilbert_value = bbox.hilbert_value();

            edge_info.bbox = bbox;
            edge_info.hilbert_value = hilbert_value;

            // Add to hilbert_to_edges
            hilbert_to_edges.insert((hilbert_value, edge_idx));

            // Add to vertex_to_edges
            for vertex in geometry.vertices_for_edge(edge_info.edge) {
                vertex_to_edges.insert((vertex, edge_idx));
            }

            // Update hilbert_to_rect mapping
            hilbert_to_rect
                .entry(hilbert_value)
                .and_modify(|existing_rect| *existing_rect = existing_rect.combine(&bbox))
                .or_insert(bbox);
        }

        // Build the R-tree
        let r_tree = RTree::build(hilbert_to_rect);

        Self {
            geometry,
            edge_info,
            extents,
            hilbert_to_edges,
            vertex_to_edges,
            dirty_set,
            r_tree,
        }
    }

    /// Insert a new edge into the spatial index
    fn insert(&mut self, edge: G::Edge, hilbert_value: u32) {
        let edge_idx = self.edge_info.len() as u32;

        let bbox = self.geometry.edge_bbox(edge, &self.extents);

        self.edge_info.push(EdgeInfo {
            edge,
            hilbert_value,
            bbox,
        });

        self.hilbert_to_edges.insert((hilbert_value, edge_idx));
        for vertex in self.geometry.vertices_for_edge(edge) {
            self.vertex_to_edges.insert((vertex, edge_idx));
        }
        self.dirty_set.insert(edge_idx);
        self.r_tree.enlarge(hilbert_value, bbox);
    }

    /// Edit an edge's endpoints
    fn replace_edge(&mut self, edge_idx: u32, new_edge: G::Edge) {
        let info = &mut self.edge_info[edge_idx as usize];
        let hilbert_value = info.hilbert_value;

        // Remove old endpoints
        for vertex in self.geometry.vertices_for_edge(info.edge) {
            self.vertex_to_edges.remove(&(vertex, edge_idx));
        }

        info.edge = new_edge;

        // Calculate new bounding box
        info.bbox = self.geometry.edge_bbox(info.edge, &self.extents);

        // Add new endpoints
        for vertex in self.geometry.vertices_for_edge(info.edge) {
            self.vertex_to_edges.insert((vertex, edge_idx));
        }

        // Mark as dirty
        self.dirty_set.insert(edge_idx);

        // Enlarge R-tree
        self.r_tree.enlarge(hilbert_value, info.bbox);
    }

    /// Get all edges that reference a given vertex
    fn edges_for_vertex(&self, vertex: G::Vertex) -> impl Iterator<Item = u32> + '_ {
        self.vertex_to_edges
            .range((vertex, 0)..=(vertex, u32::MAX))
            .map(|(_, edge_idx)| *edge_idx)
    }

    /// Get all edges with a given hilbert value
    fn edges_for_hilbert(&self, hilbert_value: u32) -> impl Iterator<Item = u32> + '_ {
        self.hilbert_to_edges
            .range((hilbert_value, 0)..=(hilbert_value, u32::MAX))
            .map(|(_, edge_idx)| *edge_idx)
    }

    /// Run the cleaning algorithm
    fn clean(&mut self) {
        while let Some(dirty_edge_idx) = self.dirty_set.pop_first() {
            // Get edge info
            let dirty_info = &self.edge_info[dirty_edge_idx as usize];

            let mut action: Option<Action<G>> = None;
            'itest: {
                // Test 1: Check for coincident vertices within this edge
                for v1 in self.geometry.vertices_for_edge(dirty_info.edge) {
                    for v2 in self.geometry.vertices_for_edge(dirty_info.edge) {
                        if v1 != v2 && self.geometry.vertices_coincident(v1, v2) {
                            action = Some(Action::MergeVertices { v1, v2 });
                            break 'itest;
                        }
                    }
                }

                // Iterate over all potentially intersecting edges
                for hilbert_value in self.r_tree.search(&dirty_info.bbox) {
                    for candidate_idx in self.edges_for_hilbert(hilbert_value) {
                        let candidate_info = &self.edge_info[candidate_idx as usize];

                        // Test 2: Check for coincident vertices
                        for v1 in self.geometry.vertices_for_edge(dirty_info.edge) {
                            for v2 in self.geometry.vertices_for_edge(candidate_info.edge) {
                                if v1 != v2 && self.geometry.vertices_coincident(v1, v2) {
                                    {
                                        action = Some(Action::MergeVertices { v1, v2 });
                                        break 'itest;
                                    }
                                }
                            }
                        }

                        // Test 3: Check for vertex on edge
                        'v_on_e: for vertex in self.geometry.vertices_for_edge(dirty_info.edge) {
                            for v2 in self.geometry.vertices_for_edge(candidate_info.edge) {
                                if vertex == v2 {
                                    // Don't do vertex-on-edge test
                                    // if this vertex is already an endpoint of the edge
                                    continue 'v_on_e;
                                }
                            }
                            if self.geometry.vertex_on_edge(vertex, candidate_info.edge) {
                                action = Some(Action::SplitEdge {
                                    edge: candidate_idx,
                                    split_vertex: vertex,
                                });
                                break 'itest;
                            }
                        }
                        'v_on_e: for vertex in self.geometry.vertices_for_edge(candidate_info.edge)
                        {
                            for v2 in self.geometry.vertices_for_edge(dirty_info.edge) {
                                if vertex == v2 {
                                    // Don't do vertex-on-edge test
                                    // if this vertex is already an endpoint of the edge
                                    continue 'v_on_e;
                                }
                            }
                            if self.geometry.vertex_on_edge(vertex, dirty_info.edge) {
                                action = Some(Action::SplitEdge {
                                    edge: dirty_edge_idx,
                                    split_vertex: vertex,
                                });
                                break 'itest;
                            }
                        }

                        // Test 4: Check for edge intersection
                        if let Some(intersection) = self
                            .geometry
                            .intersection(dirty_info.edge, candidate_info.edge)
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
            }

            // Handle the action if one was found
            match action {
                Some(Action::MergeVertices { v1, v2 }) => {
                    // Add merged vertex as a new vertex
                    let merged_vertex = self.geometry.merged_vertex(v1, v2);

                    // Update all edges referencing v1 or v2
                    for v in [v1, v2] {
                        while let Some(edge_idx) = { self.edges_for_vertex(v).next() } {
                            let mut new_edge = self.edge_info[edge_idx as usize].edge;
                            self.geometry
                                .replace_vertex_in_edge(&mut new_edge, v1, merged_vertex);
                            self.geometry
                                .replace_vertex_in_edge(&mut new_edge, v2, merged_vertex);
                            self.replace_edge(edge_idx, new_edge);
                        }
                    }
                }
                Some(Action::SplitEdge {
                    edge: edge_idx,
                    split_vertex,
                }) => {
                    let EdgeInfo {
                        edge,
                        hilbert_value,
                        ..
                    } = self.edge_info[edge_idx as usize];

                    let (new_e1, new_e2) = self.geometry.split_edge(edge, split_vertex);

                    self.replace_edge(edge_idx, new_e1);
                    self.insert(new_e2, hilbert_value);
                }
                Some(Action::SplitBothEdges {
                    edge_a: edge_a_idx,
                    edge_b: edge_b_idx,
                    intersection,
                }) => {
                    let EdgeInfo {
                        edge: edge_a,
                        hilbert_value: edge_a_hilbert,
                        ..
                    } = self.edge_info[edge_a_idx as usize];
                    let EdgeInfo {
                        edge: edge_b,
                        hilbert_value: edge_b_hilbert,
                        ..
                    } = self.edge_info[edge_b_idx as usize];

                    let vertex = self.geometry.intersection_vertex(intersection);
                    let (edge_a1, edge_a2) = self.geometry.split_edge(edge_a, vertex);
                    let (edge_b1, edge_b2) = self.geometry.split_edge(edge_b, vertex);

                    self.replace_edge(edge_a_idx, edge_a1);
                    self.replace_edge(edge_b_idx, edge_b1);
                    self.insert(edge_a2, edge_a_hilbert);
                    self.insert(edge_b2, edge_b_hilbert);
                }
                None => {}
            }
        }
    }

    /// Extract the final clean edge list
    fn extract_edges(&self) -> Vec<G::Edge> {
        self.hilbert_to_edges
            .iter()
            .map(|(_, edge_idx)| {
                let info = &self.edge_info[*edge_idx as usize];
                info.edge
            })
            .collect()
    }
}

/// Partial clean: clean edges with mixed dirty/clean flags
pub fn partial_clean<G: Geometry>(
    geometry: &mut G,
    edges: impl Iterator<Item = (G::Edge, DirtyFlag)>,
) -> Vec<G::Edge> {
    let mut spatial_index = SpatialIndex::new(geometry, edges);
    spatial_index.clean();
    spatial_index.extract_edges()
}

/// Full clean: clean all edges (assumes all dirty)
pub fn clean<G: Geometry>(geometry: &mut G, edges: impl Iterator<Item = G::Edge>) -> Vec<G::Edge> {
    partial_clean(geometry, edges.map(|edge| (edge, DirtyFlag::Dirty)))
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

    #[test]
    fn test_multiple_intersections() {
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
    fn test_complex_polygon() {
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

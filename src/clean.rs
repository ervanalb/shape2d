use crate::rtree::RTree;
use crate::{Geometry, Rect};
use std::collections::btree_map::Entry;
use std::collections::{BTreeMap, BTreeSet};
use std::iter;

/// Flag indicating whether an edge is known to be clean or might be dirty
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DirtyFlag {
    /// Edge is known to not intersect or be coincident with other clean edges
    Clean,
    /// Edge may intersect or be coincident with other edges
    Dirty,
}

/// Information about an edge
#[derive(Debug, Clone, Default)]
struct EdgeInfo {
    multiplicity: u32,
    hilbert_value: u32,
    bbox: Rect,
}

/// Actions that can be taken during the cleaning process
#[derive(Debug)]
enum Action<G: Geometry> {
    /// Cancel two opposite edges
    CancelEdges { e1: G::Edge, e2: G::Edge },
    /// Merge two coincident vertices
    MergeVertices { v1: G::Vertex, v2: G::Vertex },
    /// Merge two coincident edges
    MergeEdges { e1: G::Edge, e2: G::Edge },
    /// Split an edge at a vertex
    SplitEdge {
        edge: G::Edge,
        split_vertex: G::Vertex,
    },
    /// Split both edges at their intersection point
    SplitBothEdges {
        e1: G::Edge,
        e2: G::Edge,
        intersection: G::Intersection,
    },
}

/// Main spatial index structure
struct SpatialIndex<'a, G: Geometry> {
    geometry: &'a mut G,
    edge_info: BTreeMap<G::Edge, EdgeInfo>,
    extents: G::Extents,
    hilbert_to_edges: BTreeSet<(u32, G::Edge)>,
    vertex_to_edges: BTreeSet<(G::Vertex, G::Edge)>,
    dirty_set: BTreeSet<G::Edge>,
    r_tree: RTree,
}

impl<'a, G: Geometry> SpatialIndex<'a, G> {
    /// Build the spatial index from vertices and edges
    fn new(geometry: &'a mut G, edges: impl Iterator<Item = (G::Edge, DirtyFlag)>) -> Self {
        let mut edge_info: BTreeMap<G::Edge, EdgeInfo> = BTreeMap::new();
        let mut hilbert_to_edges = BTreeSet::new();
        let mut vertex_to_edges = BTreeSet::new();
        let mut dirty_set = BTreeSet::new();
        let mut hilbert_to_rect: BTreeMap<u32, Rect> = BTreeMap::new();

        // Loop over edges to fill out our edge_info struct partially
        for (edge, dirty_flag) in edges {
            edge_info.entry(edge).or_default().multiplicity += 1;

            // Add to dirty set if dirty
            if dirty_flag == DirtyFlag::Dirty {
                dirty_set.insert(edge);
            }
        }

        // Loop over our edges to find the extents
        let extents = geometry.extents(edge_info.keys().copied());

        // Now that we have the extents,
        // fill in the rest of the missing data
        for (&edge, edge_info) in edge_info.iter_mut() {
            // Calculate bounding box and hilbert value
            let bbox = geometry.edge_bbox(edge, &extents);
            let hilbert_value = bbox.hilbert_value();

            edge_info.bbox = bbox;
            edge_info.hilbert_value = hilbert_value;

            // Add to hilbert_to_edges
            hilbert_to_edges.insert((hilbert_value, edge));

            // Add to vertex_to_edges
            for vertex in geometry.vertices_for_edge(edge) {
                vertex_to_edges.insert((vertex, edge));
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
    fn insert(&mut self, edge: G::Edge, hilbert_value: u32, multiplicity: u32) {
        match self.edge_info.entry(edge) {
            Entry::Occupied(mut e) => {
                e.get_mut().multiplicity += multiplicity;
            }
            Entry::Vacant(e) => {
                let bbox = self.geometry.edge_bbox(edge, &self.extents);
                e.insert(EdgeInfo {
                    multiplicity,
                    hilbert_value,
                    bbox,
                });

                self.hilbert_to_edges.insert((hilbert_value, edge));
                for vertex in self.geometry.vertices_for_edge(edge) {
                    self.vertex_to_edges.insert((vertex, edge));
                }
                self.dirty_set.insert(edge);
                self.r_tree.enlarge(hilbert_value, bbox);
            }
        }
    }

    /// Remove an edge from the spatial index
    fn remove(&mut self, edge: G::Edge) -> EdgeInfo {
        let Entry::Occupied(e) = self.edge_info.entry(edge) else {
            // Edge not found
            panic!("Edge not found in edge_info");
        };
        // Remove this edge completely
        let info = e.remove();
        for vertex in self.geometry.vertices_for_edge(edge) {
            self.vertex_to_edges.remove(&(vertex, edge));
        }
        self.hilbert_to_edges.remove(&(info.hilbert_value, edge));
        self.dirty_set.remove(&edge);
        info
    }

    /// Get all edges that reference a given vertex
    fn edges_for_vertex(&self, vertex: G::Vertex) -> impl Iterator<Item = G::Edge> + '_ {
        self.vertex_to_edges
            .range((vertex, G::MIN_EDGE)..=(vertex, G::MAX_EDGE))
            .map(|&(_, e)| e)
    }

    /// Get all edges with a given hilbert value
    fn edges_for_hilbert(&self, hilbert_value: u32) -> impl Iterator<Item = G::Edge> + '_ {
        self.hilbert_to_edges
            .range((hilbert_value, G::MIN_EDGE)..=(hilbert_value, G::MAX_EDGE))
            .map(|&(_, e)| e)
    }

    /// Run the cleaning algorithm
    fn clean(&mut self) {
        while let Some(dirty_edge) = self.dirty_set.pop_first() {
            // Get edge info
            let dirty_info = &self.edge_info.get(&dirty_edge).unwrap();

            let mut action: Option<Action<G>> = None;
            'itest: {
                // Test 1: Check for coincident vertices within this edge
                for v1 in self.geometry.vertices_for_edge(dirty_edge) {
                    for v2 in self.geometry.vertices_for_edge(dirty_edge) {
                        if v1 != v2 && self.geometry.vertices_coincident(v1, v2) {
                            action = Some(Action::MergeVertices { v1, v2 });
                            break 'itest;
                        }
                    }
                }

                // Iterate over all potentially intersecting edges
                for hilbert_value in self.r_tree.search(&dirty_info.bbox) {
                    for candidate_edge in self.edges_for_hilbert(hilbert_value) {
                        if dirty_edge == candidate_edge {
                            // Don't bother checking this edge against itself
                            continue;
                        }

                        // Test 2: Check if these edges cancel
                        if self.geometry.edges_cancel(dirty_edge, candidate_edge) {
                            action = Some(Action::CancelEdges {
                                e1: dirty_edge,
                                e2: candidate_edge,
                            });
                            break 'itest;
                        }

                        // Test 3: Check if these two edges are fully coincident
                        if self.geometry.edges_coincident(dirty_edge, candidate_edge) {
                            action = Some(Action::MergeEdges {
                                e1: dirty_edge,
                                e2: candidate_edge,
                            });
                            break 'itest;
                        }

                        // Test 4: Check for coincident vertices
                        for v1 in self.geometry.vertices_for_edge(dirty_edge) {
                            for v2 in self.geometry.vertices_for_edge(candidate_edge) {
                                if v1 != v2 && self.geometry.vertices_coincident(v1, v2) {
                                    {
                                        action = Some(Action::MergeVertices { v1, v2 });
                                        break 'itest;
                                    }
                                }
                            }
                        }

                        // Test 5: Check for vertex on edge
                        'v_on_e: for vertex in self.geometry.vertices_for_edge(dirty_edge) {
                            for v2 in self.geometry.vertices_for_edge(candidate_edge) {
                                if vertex == v2 {
                                    // Don't do vertex-on-edge test
                                    // if this vertex is already an endpoint of the edge
                                    continue 'v_on_e;
                                }
                            }
                            if self.geometry.vertex_on_edge(vertex, candidate_edge) {
                                action = Some(Action::SplitEdge {
                                    edge: candidate_edge,
                                    split_vertex: vertex,
                                });
                                break 'itest;
                            }
                        }
                        'v_on_e: for vertex in self.geometry.vertices_for_edge(candidate_edge) {
                            for v2 in self.geometry.vertices_for_edge(dirty_edge) {
                                if vertex == v2 {
                                    // Don't do vertex-on-edge test
                                    // if this vertex is already an endpoint of the edge
                                    continue 'v_on_e;
                                }
                            }
                            if self.geometry.vertex_on_edge(vertex, dirty_edge) {
                                action = Some(Action::SplitEdge {
                                    edge: dirty_edge,
                                    split_vertex: vertex,
                                });
                                break 'itest;
                            }
                        }

                        // Test 6: Check for edge intersection
                        if let Some(intersection) =
                            self.geometry.intersection(dirty_edge, candidate_edge)
                        {
                            action = Some(Action::SplitBothEdges {
                                e1: dirty_edge,
                                e2: candidate_edge,
                                intersection,
                            });
                            break 'itest;
                        }
                    }
                }
            }

            // Handle the action if one was found
            match action {
                Some(Action::CancelEdges { e1, e2 }) => {
                    // Remove both edges
                    let old_1 = self.remove(e1);
                    let old_2 = self.remove(e2);

                    // See if we need to add one back due to a multiplicity imbalance
                    if old_1.multiplicity > old_2.multiplicity {
                        self.insert(
                            e1,
                            old_1.hilbert_value,
                            old_1.multiplicity - old_2.multiplicity,
                        );
                    } else if old_2.multiplicity > old_1.multiplicity {
                        self.insert(
                            e2,
                            old_2.hilbert_value,
                            old_2.multiplicity - old_1.multiplicity,
                        );
                    }
                }
                Some(Action::MergeVertices { v1, v2 }) => {
                    // Add merged vertex as a new vertex
                    let merged_vertex = self.geometry.merged_vertex(v1, v2);

                    // Update all edges referencing v1 or v2
                    for v in [v1, v2] {
                        while let Some(edge) = { self.edges_for_vertex(v).next() } {
                            let new_edge = edge;
                            let old = self.remove(edge);
                            let new_edge =
                                self.geometry
                                    .replace_vertex_in_edge(new_edge, v1, merged_vertex);
                            if let Some(new_edge) = new_edge {
                                let new_edge = self.geometry.replace_vertex_in_edge(
                                    new_edge,
                                    v2,
                                    merged_vertex,
                                );
                                if let Some(new_edge) = new_edge {
                                    self.insert(new_edge, old.hilbert_value, old.multiplicity);
                                }
                            }
                        }
                    }
                }
                Some(Action::MergeEdges { e1, e2 }) => {
                    // Construct merged edges
                    let (new_e1, new_e2) = self.geometry.merged_edges(e1, e2);

                    // Remove old edges
                    let old_e1 = self.remove(e1);
                    let old_e2 = self.remove(e2);

                    // Insert new edges
                    // (Doesn't really matter which hilbert value we use,
                    // or even if we use different ones,
                    // but let's put them into the same hilbert value)
                    self.insert(new_e1, old_e1.hilbert_value, old_e1.multiplicity);
                    self.insert(new_e2, old_e1.hilbert_value, old_e2.multiplicity);
                }
                Some(Action::SplitEdge { edge, split_vertex }) => {
                    let (new_e1, new_e2) = self.geometry.split_edge(edge, split_vertex);
                    let old = self.remove(edge);
                    self.insert(new_e1, old.hilbert_value, old.multiplicity);
                    self.insert(new_e2, old.hilbert_value, old.multiplicity);
                }
                Some(Action::SplitBothEdges {
                    e1: edge_a,
                    e2: edge_b,
                    intersection,
                }) => {
                    let vertex = self.geometry.intersection_vertex(intersection);
                    let (edge_a1, edge_a2) = self.geometry.split_edge(edge_a, vertex);
                    let (edge_b1, edge_b2) = self.geometry.split_edge(edge_b, vertex);

                    let old_a = self.remove(edge_a);
                    self.insert(edge_a1, old_a.hilbert_value, old_a.multiplicity);
                    self.insert(edge_a2, old_a.hilbert_value, old_a.multiplicity);
                    let old_b = self.remove(edge_b);
                    self.insert(edge_b1, old_b.hilbert_value, old_b.multiplicity);
                    self.insert(edge_b2, old_b.hilbert_value, old_b.multiplicity);
                }
                None => {}
            }
        }
    }

    /// Extract the final clean edge list
    fn extract_edges(&self) -> Vec<G::Edge> {
        self.edge_info
            .iter()
            .flat_map(|(&edge, info)| iter::repeat(edge).take(info.multiplicity as usize))
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
    use crate::geometry::MyGeometry;

    #[test]
    fn test_simple_two_edges() {
        let mut geometry = MyGeometry::new(vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]);

        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut geometry, edges.into_iter());
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_intersecting_edges() {
        let mut geometry = MyGeometry::new(vec![[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]]);

        // Two edges that intersect
        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut geometry, edges.into_iter());
        // Should split into 4 edges
        assert_eq!(result.len(), 4);
        // Should create a new vertex at the intersection
        assert_eq!(geometry.vertex_count(), 5);
    }

    #[test]
    fn test_vertex_on_edge() {
        let mut geometry = MyGeometry::new(vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 0.0], // On the first edge
        ]);

        let edges = [(0, 1)];

        let result = clean(&mut geometry, edges.into_iter());
        // Should not split since there's only one edge and the vertex isn't connected
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_coincident_vertices() {
        let mut geometry = MyGeometry::new(vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 0.0], // Coincident with vertex 0
            [1.0, 1.0],
        ]);

        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut geometry, edges.into_iter());
        // Should merge vertices 0 and 2
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_square_polygon() {
        let mut geometry = MyGeometry::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);

        let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];

        let result = clean(&mut geometry, edges.into_iter());
        // Square edges don't intersect, should remain 4 edges
        assert_eq!(result.len(), 4);
        // No new vertices should be created
        assert_eq!(geometry.vertex_count(), 4);
    }

    #[test]
    fn test_t_junction() {
        let mut geometry = MyGeometry::new(vec![[0.0, 0.5], [1.0, 0.5], [0.5, 0.0], [0.5, 1.0]]);

        // T-junction: horizontal edge intersected by vertical edge
        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut geometry, edges.into_iter());
        // Should split into 4 edges meeting at center
        assert_eq!(result.len(), 4);
        // Should create intersection vertex at (0.5, 0.5)
        assert_eq!(geometry.vertex_count(), 5);
    }

    #[test]
    fn test_multiple_intersections() {
        let mut geometry = MyGeometry::new(vec![
            // Horizontal line
            [0.0, 0.5],
            [1.0, 0.5],
            // Vertical line
            [0.5, 0.0],
            [0.5, 1.0],
            // Diagonal line
            [0.0, 0.0],
            [1.0, 1.0],
        ]);

        let edges = [(0, 1), (2, 3), (4, 5)];

        let result = clean(&mut geometry, edges.into_iter());
        // Each of the 3 edges should be split at their intersections
        // Horizontal and vertical intersect at [0.5, 0.5]
        // Diagonal intersects both at [0.5, 0.5]
        // So all three meet at one point, creating 6 edges
        assert_eq!(result.len(), 6);
    }

    #[test]
    fn test_complex_polygon() {
        // Star-like pattern
        let mut geometry = MyGeometry::new(vec![
            [0.5, 0.0], // Top
            [0.7, 0.5], // Right
            [0.5, 1.0], // Bottom
            [0.3, 0.5], // Left
            [0.5, 0.5], // Center
        ]);

        let edges = [
            (0, 2), // Vertical through center
            (1, 3), // Horizontal through center
        ];

        let result = clean(&mut geometry, edges.into_iter());
        // Both edges pass through center, should be split there
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_canceling_edges() {
        let mut geometry = MyGeometry::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);

        // Two edges with same endpoints but opposite directions should cancel
        let edges = [
            (0, 1), // Forward edge
            (1, 0), // Reverse edge (cancels with first)
            (1, 2), // Different edge (should remain)
        ];

        let result = clean(&mut geometry, edges.into_iter());
        // The two canceling edges should be removed, leaving only one
        assert_eq!(result.len(), 1);
        assert!(result.contains(&(1, 2)));
    }

    #[test]
    fn test_merged_vertices_remove_short_edge() {
        let mut geometry = MyGeometry::new(vec![
            [0.0, 0.0],
            [0.0, 0.0], // Coincident with vertex 0
            [1.0, 0.0],
        ]);

        // Edge from vertex 0 to vertex 1 (which are coincident)
        // should disappear when vertices are merged
        let edges = [
            (0, 1), // Very short edge between coincident vertices
            (1, 2), // Regular edge
        ];

        let result = clean(&mut geometry, edges.into_iter());
        // The short edge should be removed during vertex merging
        // Only the regular edge should remain
        assert_eq!(result.len(), 1);
        // The edge should connect the merged vertex to vertex 2
        assert!(result[0].1 == 2 || result[0].0 == 2);
    }
}

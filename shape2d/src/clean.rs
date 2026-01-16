use crate::kernel::{Edge, Kernel};
use crate::rtree::{RTree, Rect};
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
enum Action<K: Kernel> {
    /// Cancel two opposite edges
    CancelEdges { e1: K::Edge, e2: K::Edge },
    /// Merge two coincident vertices
    MergeVertices { v1: K::Vertex, v2: K::Vertex },
    /// Merge two coincident edges
    MergeEdges { e1: K::Edge, e2: K::Edge },
    /// Split an edge at a vertex
    SplitEdge {
        edge: K::Edge,
        split_vertex: K::Vertex,
        pt: K::Point,
    },
    /// Split both edges at their intersection point
    SplitBothEdges {
        e1: K::Edge,
        e2: K::Edge,
        pt: K::Point,
    },
}

impl<K: Kernel> std::fmt::Debug for Action<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Action::CancelEdges { e1, e2 } => f
                .debug_struct("CancelEdges")
                .field("e1", e1)
                .field("e2", e2)
                .finish(),
            Action::MergeVertices { v1, v2 } => f
                .debug_struct("MergeVertices")
                .field("v1", v1)
                .field("v2", v2)
                .finish(),
            Action::MergeEdges { e1, e2 } => f
                .debug_struct("MergeEdges")
                .field("e1", e1)
                .field("e2", e2)
                .finish(),
            Action::SplitEdge {
                edge,
                split_vertex,
                pt: _,
            } => f
                .debug_struct("SplitEdge")
                .field("edge", edge)
                .field("split_vertex", split_vertex)
                .finish(),
            Action::SplitBothEdges { e1, e2, pt: _ } => f
                .debug_struct("SplitBothEdges")
                .field("e1", e1)
                .field("e2", e2)
                .finish(),
        }
    }
}

/// Main spatial index structure
struct SpatialIndex<'a, K: Kernel> {
    kernel: &'a mut K,
    edge_info: BTreeMap<K::Edge, EdgeInfo>,
    extents: K::Extents,
    hilbert_to_edges: BTreeSet<(u32, K::Edge)>,
    vertex_to_edges: BTreeSet<(K::Vertex, K::Edge)>,
    dirty_set: BTreeSet<K::Edge>,
    r_tree: RTree,
}

impl<'a, K: Kernel> SpatialIndex<'a, K> {
    /// Build the spatial index from vertices and edges
    fn new(kernel: &'a mut K, edges: impl Iterator<Item = (K::Edge, DirtyFlag)>) -> Self {
        let mut edge_info: BTreeMap<K::Edge, EdgeInfo> = BTreeMap::new();
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
        let extents = kernel.extents(edge_info.keys().copied());

        // Now that we have the extents,
        // fill in the rest of the missing data
        for (&edge, edge_info) in edge_info.iter_mut() {
            // Calculate bounding box and hilbert value
            let bbox = kernel.edge_bbox(edge, &extents);
            let hilbert_value = bbox.hilbert_value();

            edge_info.bbox = bbox;
            edge_info.hilbert_value = hilbert_value;

            // Add to hilbert_to_edges
            hilbert_to_edges.insert((hilbert_value, edge));

            // Add to vertex_to_edges
            if let Some((start_v, end_v)) = kernel.vertices_for_edge(edge) {
                vertex_to_edges.insert((start_v, edge));
                vertex_to_edges.insert((end_v, edge));
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
            kernel,
            edge_info,
            extents,
            hilbert_to_edges,
            vertex_to_edges,
            dirty_set,
            r_tree,
        }
    }

    /// Insert a new edge into the spatial index
    fn insert(&mut self, edge: K::Edge, hilbert_value: u32, multiplicity: u32) {
        match self.edge_info.entry(edge) {
            Entry::Occupied(mut e) => {
                e.get_mut().multiplicity += multiplicity;
            }
            Entry::Vacant(e) => {
                let bbox = self.kernel.edge_bbox(edge, &self.extents);
                e.insert(EdgeInfo {
                    multiplicity,
                    hilbert_value,
                    bbox,
                });

                self.hilbert_to_edges.insert((hilbert_value, edge));
                if let Some((start_v, end_v)) = self.kernel.vertices_for_edge(edge) {
                    self.vertex_to_edges.insert((start_v, edge));
                    self.vertex_to_edges.insert((end_v, edge));
                }
                self.r_tree.enlarge(hilbert_value, bbox);
                self.dirty_set.insert(edge);
            }
        }
    }

    /// Remove an edge from the spatial index
    fn remove(&mut self, edge: K::Edge) -> EdgeInfo {
        let Entry::Occupied(e) = self.edge_info.entry(edge) else {
            // Edge not found
            panic!("Edge not found in edge_info");
        };
        // Remove this edge completely
        let info = e.remove();
        if let Some((start_v, end_v)) = self.kernel.vertices_for_edge(edge) {
            self.vertex_to_edges.remove(&(start_v, edge));
            self.vertex_to_edges.remove(&(end_v, edge));
        }
        self.hilbert_to_edges.remove(&(info.hilbert_value, edge));
        self.dirty_set.remove(&edge);
        info
    }

    /// Get all edges that reference a given vertex
    fn edges_for_vertex(&self, vertex: K::Vertex) -> impl Iterator<Item = K::Edge> + '_ {
        self.vertex_to_edges
            .range((vertex, K::Edge::MIN)..=(vertex, K::Edge::MAX))
            .map(|&(_, e)| e)
    }

    /// Get all edges with a given hilbert value
    fn edges_for_hilbert(&self, hilbert_value: u32) -> impl Iterator<Item = K::Edge> + '_ {
        self.hilbert_to_edges
            .range((hilbert_value, K::Edge::MIN)..=(hilbert_value, K::Edge::MAX))
            .map(|&(_, e)| e)
    }

    /// Run the cleaning algorithm
    fn clean(&mut self) {
        while let Some(&dirty_edge) = self.dirty_set.first() {
            // Get edge info
            let dirty_info = &self.edge_info.get(&dirty_edge).unwrap();

            let mut action: Option<Action<K>> = None;
            'itest: {
                // Test 1: Check for coincident vertices within this edge
                if let Some((v1, v2)) = self.kernel.vertices_for_edge(dirty_edge) {
                    if v1 != v2 && self.kernel.vertices_coincident(v1, v2) {
                        action = Some(Action::MergeVertices { v1, v2 });
                        break 'itest;
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
                        if dirty_edge.reversed() == candidate_edge {
                            action = Some(Action::CancelEdges {
                                e1: dirty_edge,
                                e2: candidate_edge,
                            });
                            break 'itest;
                        }

                        // Test 3: Check if these two edges are fully coincident
                        if self.kernel.edges_coincident(dirty_edge, candidate_edge) {
                            action = Some(Action::MergeEdges {
                                e1: dirty_edge,
                                e2: candidate_edge,
                            });
                            break 'itest;
                        }

                        // Test 4: Check for coincident vertices
                        if let Some((e1_start, e1_end)) = self.kernel.vertices_for_edge(dirty_edge)
                            && let Some((e2_start, e2_end)) =
                                self.kernel.vertices_for_edge(candidate_edge)
                        {
                            for (v1, v2) in [
                                (e1_start, e2_start),
                                (e1_start, e2_end),
                                (e1_end, e2_start),
                                (e1_end, e2_end),
                            ] {
                                if v1 != v2 && self.kernel.vertices_coincident(v1, v2) {
                                    {
                                        action = Some(Action::MergeVertices { v1, v2 });
                                        break 'itest;
                                    }
                                }
                            }
                        }

                        // Test 5: Check for vertex on edge

                        if let Some((dirty_edge_start, dirty_edge_end)) =
                            self.kernel.vertices_for_edge(dirty_edge)
                            && let Some((candidate_edge_start, candidate_edge_end)) =
                                self.kernel.vertices_for_edge(candidate_edge)
                        {
                            'v_on_e: for vertex in [dirty_edge_start, dirty_edge_end] {
                                for v2 in [candidate_edge_start, candidate_edge_end] {
                                    if vertex == v2 {
                                        // Don't do vertex-on-edge test
                                        // if this vertex is already an endpoint of the edge
                                        continue 'v_on_e;
                                    }
                                }
                                if let Some(pt) = self.kernel.vertex_on_edge(vertex, candidate_edge)
                                {
                                    action = Some(Action::SplitEdge {
                                        edge: candidate_edge,
                                        split_vertex: vertex,
                                        pt,
                                    });
                                    break 'itest;
                                }
                            }
                            'v_on_e: for vertex in [candidate_edge_start, candidate_edge_end] {
                                for v2 in [dirty_edge_start, dirty_edge_end] {
                                    if vertex == v2 {
                                        // Don't do vertex-on-edge test
                                        // if this vertex is already an endpoint of the edge
                                        continue 'v_on_e;
                                    }
                                }
                                if let Some(pt) = self.kernel.vertex_on_edge(vertex, dirty_edge) {
                                    action = Some(Action::SplitEdge {
                                        edge: dirty_edge,
                                        split_vertex: vertex,
                                        pt,
                                    });
                                    break 'itest;
                                }
                            }
                        }

                        // Test 6: Check for edge intersection
                        if let Some(intersection) =
                            self.kernel.intersection(dirty_edge, candidate_edge)
                        {
                            action = Some(Action::SplitBothEdges {
                                e1: dirty_edge,
                                e2: candidate_edge,
                                pt: intersection,
                            });
                            break 'itest;
                        }
                    }
                }

                // If all tests pass, remove this edge from the dirty set
                self.dirty_set.pop_first();
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
                    let merged_vertex = self.kernel.merged_vertex(v1, v2);

                    // Update all edges referencing v1 or v2
                    for v in [v1, v2] {
                        while let Some(edge) = { self.edges_for_vertex(v).next() } {
                            let new_edge = edge;
                            let old = self.remove(edge);
                            let new_edge =
                                self.kernel
                                    .replace_vertex_in_edge(new_edge, v1, merged_vertex);
                            if let Some(new_edge) = new_edge {
                                let new_edge =
                                    self.kernel
                                        .replace_vertex_in_edge(new_edge, v2, merged_vertex);
                                if let Some(new_edge) = new_edge {
                                    self.insert(new_edge, old.hilbert_value, old.multiplicity);
                                }
                            }
                        }
                    }
                }
                Some(Action::MergeEdges { e1, e2 }) => {
                    // Construct merged edges
                    let (new_e1, new_e2) = self.kernel.merged_edges(e1, e2);

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
                Some(Action::SplitEdge {
                    edge,
                    split_vertex: old_split_vertex,
                    pt,
                }) => {
                    // Add split point as a new vertex
                    let new_split_vertex = self.kernel.push_vertex(pt);

                    // Update all edges referencing old_split_vertex
                    while let Some(edge) = { self.edges_for_vertex(old_split_vertex).next() } {
                        let new_edge = edge;
                        let old = self.remove(edge);
                        if let Some(new_edge) = self.kernel.replace_vertex_in_edge(
                            new_edge,
                            old_split_vertex,
                            new_split_vertex,
                        ) {
                            self.insert(new_edge, old.hilbert_value, old.multiplicity);
                        }
                    }

                    // Split the edge
                    let (new_e1, new_e2) = self.kernel.split_edge(edge, new_split_vertex);
                    let old = self.remove(edge);
                    self.insert(new_e1, old.hilbert_value, old.multiplicity);
                    self.insert(new_e2, old.hilbert_value, old.multiplicity);
                }
                Some(Action::SplitBothEdges {
                    e1: edge_a,
                    e2: edge_b,
                    pt: intersection,
                }) => {
                    let vertex = self.kernel.push_vertex(intersection);
                    let (edge_a1, edge_a2) = self.kernel.split_edge(edge_a, vertex);
                    let (edge_b1, edge_b2) = self.kernel.split_edge(edge_b, vertex);

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
    fn extract_edges(&self) -> Vec<K::Edge> {
        self.edge_info
            .iter()
            .flat_map(|(&edge, info)| iter::repeat(edge).take(info.multiplicity as usize))
            .collect()
    }
}

/// Partial clean: clean edges with mixed dirty/clean flags
pub fn partial_clean<K: Kernel>(
    kernel: &mut K,
    edges: impl Iterator<Item = (K::Edge, DirtyFlag)>,
) -> Vec<K::Edge> {
    let mut spatial_index = SpatialIndex::new(kernel, edges);
    spatial_index.clean();
    spatial_index.extract_edges()
}

/// Full clean: clean all edges (assumes all dirty)
pub fn clean<K: Kernel>(kernel: &mut K, edges: impl Iterator<Item = K::Edge>) -> Vec<K::Edge> {
    partial_clean(kernel, edges.map(|edge| (edge, DirtyFlag::Dirty)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::polyline::F32 as Kernel;

    #[test]
    fn test_simple_two_edges() {
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]);

        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut kernel, edges.into_iter());
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_intersecting_edges() {
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]]);

        // Two edges that intersect
        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut kernel, edges.into_iter());
        // Should split into 4 edges
        assert_eq!(result.len(), 4);
        // Should create a new vertex at the intersection
        assert_eq!(kernel.vertices.len(), 5);
    }

    #[test]
    fn test_vertex_on_edge() {
        let mut kernel = Kernel::new(vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 0.0], // On the first edge
        ]);

        let edges = [(0, 1)];

        let result = clean(&mut kernel, edges.into_iter());
        // Should not split since there's only one edge and the vertex isn't connected
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_coincident_vertices() {
        let mut kernel = Kernel::new(vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 0.0], // Coincident with vertex 0
            [1.0, 1.0],
        ]);

        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut kernel, edges.into_iter());
        // Should merge vertices 0 and 2
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_square_polygon() {
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);

        let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];

        let result = clean(&mut kernel, edges.into_iter());
        // Square edges don't intersect, should remain 4 edges
        assert_eq!(result.len(), 4);
        // No new vertices should be created
        assert_eq!(kernel.vertices.len(), 4);
    }

    #[test]
    fn test_t_junction() {
        let mut kernel = Kernel::new(vec![[0.0, 0.5], [1.0, 0.5], [0.5, 0.0], [0.5, 1.0]]);

        // T-junction: horizontal edge intersected by vertical edge
        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut kernel, edges.into_iter());
        // Should split into 4 edges meeting at center
        assert_eq!(result.len(), 4);
        // Should create intersection vertex at (0.5, 0.5)
        assert_eq!(kernel.vertices.len(), 5);
    }

    #[test]
    fn test_multiple_intersections() {
        let mut kernel = Kernel::new(vec![
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

        let result = clean(&mut kernel, edges.into_iter());
        // Each of the 3 edges should be split at their intersections
        // Horizontal and vertical intersect at [0.5, 0.5]
        // Diagonal intersects both at [0.5, 0.5]
        // So all three meet at one point, creating 6 edges
        dbg!(&result);
        assert_eq!(result.len(), 6);
    }

    #[test]
    fn test_complex_polygon() {
        // Star-like pattern
        let mut kernel = Kernel::new(vec![
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

        let result = clean(&mut kernel, edges.into_iter());
        // Both edges pass through center, should be split there
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_canceling_edges() {
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);

        // Two edges with same endpoints but opposite directions should cancel
        let edges = [
            (0, 1), // Forward edge
            (1, 0), // Reverse edge (cancels with first)
            (1, 2), // Different edge (should remain)
        ];

        let result = clean(&mut kernel, edges.into_iter());
        // The two canceling edges should be removed, leaving only one
        assert_eq!(result.len(), 1);
        assert!(result.contains(&(1, 2)));
    }

    #[test]
    fn test_merged_vertices_remove_short_edge() {
        let mut kernel = Kernel::new(vec![
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

        let result = clean(&mut kernel, edges.into_iter());
        // The short edge should be removed during vertex merging
        // Only the regular edge should remain
        assert_eq!(result.len(), 1);
        // The edge should connect the merged vertex to vertex 2
        assert!(result[0].1 == 2 || result[0].0 == 2);
    }

    #[test]
    fn test_complete_clean() {
        let mut kernel = Kernel::new_with_epsilon(
            vec![
                [0.0, 0.0],
                [2.0, 0.0],
                [2.0, 2.0],
                [1.0, 1.0],
                [3.5337768, 0.99069494],
                [3.118609, 2.8236988],
            ],
            0.1,
        );

        let edges = vec![(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)];

        let edges = clean(&mut kernel, edges.iter().copied());
        // Cleaning twice shouldn't change the result
        let edges2 = clean(&mut kernel, edges.iter().copied());
        assert_eq!(edges, edges2);
    }

    #[test]
    fn test_halting() {
        let mut kernel = Kernel::new_with_epsilon(
            vec![
                [0.04894346, 3.309017],
                [1.0, 3.6618032],
                [0.8928008, 2.950759],
                [2.0, 3.4],
            ],
            0.7,
        );

        let edges = vec![(0, 1), (1, 3), (2, 3)];

        clean(&mut kernel, edges.iter().copied());
        // Clean shouldn't run forever
    }
}

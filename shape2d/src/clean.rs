use crate::kernel::{Edge, EdgeInteraction, Kernel};
use crate::rtree::{RTree, Rect};
use std::cmp::Ordering;
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
    MergeVertices {
        v1: K::Vertex,
        v2: K::Vertex,
        pt: K::MergePoint,
    },
    /// Merge two coincident edges
    MergeEdges {
        e1: K::Edge,
        e2: K::Edge,
        cancel: bool,
        curve: K::MergeCurve,
    },
    /// Split an edge at a vertex
    SplitEdge {
        edge: K::Edge,
        split_vertex: K::Vertex,
        pt: K::SplitPoint,
    },
    /// Split both edges at their intersection point
    SplitBothEdges {
        e1: K::Edge,
        e2: K::Edge,
        pt: K::IntersectionPoint,
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
            Action::MergeVertices { v1, v2, pt: _ } => f
                .debug_struct("MergeVertices")
                .field("v1", v1)
                .field("v2", v2)
                .finish(),
            Action::MergeEdges {
                e1,
                e2,
                cancel,
                curve: _,
            } => f
                .debug_struct("MergeEdges")
                .field("e1", e1)
                .field("e2", e2)
                .field("cancel", cancel)
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
struct SpatialIndex<K: Kernel> {
    edge_info: BTreeMap<K::Edge, EdgeInfo>,
    extents: K::Extents,
    hilbert_to_edges: BTreeSet<(u32, K::Edge)>,
    vertex_to_edges: BTreeSet<(K::Vertex, K::Edge)>,
    dirty_set: BTreeSet<K::Edge>,
    r_tree: RTree,
}

impl<K: Kernel> SpatialIndex<K> {
    /// Build the spatial index from vertices and edges
    fn new(kernel: &mut K, edges: impl Iterator<Item = (K::Edge, DirtyFlag)>) -> Self {
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
            edge_info,
            extents,
            hilbert_to_edges,
            vertex_to_edges,
            dirty_set,
            r_tree,
        }
    }

    fn rebuild(&mut self, kernel: &mut K) {
        let edge_info = std::mem::replace(&mut self.edge_info, Default::default());
        let dirty_set = std::mem::replace(&mut self.dirty_set, Default::default());

        let edges = edge_info.iter().flat_map(|(&edge, info)| {
            iter::repeat((
                edge,
                match dirty_set.contains(&edge) {
                    true => DirtyFlag::Dirty,
                    false => DirtyFlag::Clean,
                },
            ))
            .take(info.multiplicity as usize)
        });

        *self = Self::new(kernel, edges);
    }

    /// Insert a new edge into the spatial index
    fn insert(&mut self, kernel: &mut K, edge: K::Edge, hilbert_value: u32, multiplicity: u32) {
        match self.edge_info.entry(edge) {
            Entry::Occupied(mut e) => {
                e.get_mut().multiplicity += multiplicity;
            }
            Entry::Vacant(e) => {
                let bbox = kernel.edge_bbox(edge, &self.extents);
                e.insert(EdgeInfo {
                    multiplicity,
                    hilbert_value,
                    bbox,
                });

                self.hilbert_to_edges.insert((hilbert_value, edge));
                if let Some((start_v, end_v)) = kernel.vertices_for_edge(edge) {
                    self.vertex_to_edges.insert((start_v, edge));
                    self.vertex_to_edges.insert((end_v, edge));
                }
                self.r_tree.enlarge(hilbert_value, bbox);
                self.dirty_set.insert(edge);
            }
        }
    }

    /// Remove an edge from the spatial index
    fn remove(&mut self, kernel: &mut K, edge: K::Edge) -> EdgeInfo {
        let Entry::Occupied(e) = self.edge_info.entry(edge) else {
            // Edge not found
            panic!("Edge not found in edge_info");
        };
        // Remove this edge completely
        let info = e.remove();
        if let Some((start_v, end_v)) = kernel.vertices_for_edge(edge) {
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

    fn clean(&mut self, kernel: &mut K) {
        let entities = self.edge_info.len();
        let mut i = 0;
        println!("# of entities: {}", entities);
        'cleaning: loop {
            for _ in 0..=entities {
                self.clean_one_thing(kernel);
                i += 1;
                if self.dirty_set.is_empty() {
                    break 'cleaning;
                }
            }
            // After cleaning 2x number of entities, rebuild the Rtree
            println!("Rebuild!");
            self.rebuild(kernel);
        }
        println!("Cleaned, took {} steps", i);
    }

    /// Run the cleaning algorithm
    fn clean_one_thing(&mut self, kernel: &mut K) {
        while let Some(&dirty_edge) = self.dirty_set.first() {
            // Get edge info
            let dirty_info = &self.edge_info.get(&dirty_edge).unwrap();

            let mut action: Option<Action<K>> = None;
            'itest: {
                // Test 1: Check for coincident vertices within this edge
                if let Some((v1, v2)) = kernel.vertices_for_edge(dirty_edge) {
                    if v1 != v2
                        && let Some(pt) = kernel.vertices_coincident(v1, v2)
                    {
                        action = Some(Action::MergeVertices { pt, v1, v2 });
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
                        match kernel.edges_coincident(dirty_edge, candidate_edge) {
                            EdgeInteraction::Merge(curve) => {
                                action = Some(Action::MergeEdges {
                                    e1: dirty_edge,
                                    e2: candidate_edge,
                                    cancel: false,
                                    curve,
                                });
                                break 'itest;
                            }
                            EdgeInteraction::Cancel(curve) => {
                                action = Some(Action::MergeEdges {
                                    e1: dirty_edge,
                                    e2: candidate_edge,
                                    cancel: true,
                                    curve,
                                });
                                break 'itest;
                            }
                            EdgeInteraction::None => {}
                        }

                        // Test 4: Check for coincident vertices
                        if let Some((e1_start, e1_end)) = kernel.vertices_for_edge(dirty_edge)
                            && let Some((e2_start, e2_end)) =
                                kernel.vertices_for_edge(candidate_edge)
                        {
                            for (v1, v2) in [
                                (e1_start, e2_start),
                                (e1_start, e2_end),
                                (e1_end, e2_start),
                                (e1_end, e2_end),
                            ] {
                                if v1 != v2
                                    && let Some(pt) = kernel.vertices_coincident(v1, v2)
                                {
                                    {
                                        action = Some(Action::MergeVertices { v1, v2, pt });
                                        break 'itest;
                                    }
                                }
                            }
                        }

                        // Test 5: Check for vertex on edge

                        if let Some((dirty_edge_start, dirty_edge_end)) =
                            kernel.vertices_for_edge(dirty_edge)
                            && let Some((candidate_edge_start, candidate_edge_end)) =
                                kernel.vertices_for_edge(candidate_edge)
                        {
                            'v_on_e: for vertex in [dirty_edge_start, dirty_edge_end] {
                                for v2 in [candidate_edge_start, candidate_edge_end] {
                                    if vertex == v2 {
                                        // Don't do vertex-on-edge test
                                        // if this vertex is already an endpoint of the edge
                                        continue 'v_on_e;
                                    }
                                }
                                if let Some(pt) = kernel.split(vertex, candidate_edge) {
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
                                if let Some(pt) = kernel.split(vertex, dirty_edge) {
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
                        if let Some(intersection) = kernel.intersection(dirty_edge, candidate_edge)
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

            let Some(action) = action else {
                continue;
            };
            match action {
                Action::CancelEdges { e1, e2 } => {
                    // Remove both edges
                    let old_1 = self.remove(kernel, e1);
                    let old_2 = self.remove(kernel, e2);

                    // See if we need to add one back due to a multiplicity imbalance
                    if old_1.multiplicity > old_2.multiplicity {
                        self.insert(
                            kernel,
                            e1,
                            old_1.hilbert_value,
                            old_1.multiplicity - old_2.multiplicity,
                        );
                    } else if old_2.multiplicity > old_1.multiplicity {
                        self.insert(
                            kernel,
                            e2,
                            old_2.hilbert_value,
                            old_2.multiplicity - old_1.multiplicity,
                        );
                    }
                }
                Action::MergeVertices { v1, v2, pt } => {
                    // Add merged vertex as a new vertex
                    let merged_vertex = kernel.merge_vertices(pt);

                    // Update all edges referencing v1 or v2
                    for v in [v1, v2] {
                        while let Some(edge) = { self.edges_for_vertex(v).next() } {
                            let new_edge = edge;
                            let old = self.remove(kernel, edge);
                            let new_edge =
                                kernel.replace_vertex_in_edge(new_edge, v1, merged_vertex);
                            if let Some(new_edge) = new_edge {
                                let new_edge =
                                    kernel.replace_vertex_in_edge(new_edge, v2, merged_vertex);
                                if let Some(new_edge) = new_edge {
                                    self.insert(
                                        kernel,
                                        new_edge,
                                        old.hilbert_value,
                                        old.multiplicity,
                                    );
                                }
                            }
                        }
                    }
                }
                Action::MergeEdges {
                    e1,
                    e2,
                    cancel,
                    curve,
                } => {
                    // Construct merged edge
                    let new_e = kernel.merge_edges(curve);

                    // Remove old edges
                    let old_e1 = self.remove(kernel, e1);
                    let old_e2 = self.remove(kernel, e2);

                    // Insert new edges based on old multiplicities,
                    // and whether they add or cancel.
                    // Note: It doesn't really matter which old edge we choose
                    // to take the hilbert value from
                    match cancel {
                        false => {
                            self.insert(
                                kernel,
                                new_e,
                                old_e1.hilbert_value,
                                old_e1.multiplicity + old_e2.multiplicity,
                            );
                        }
                        true => {
                            match old_e1.multiplicity.cmp(&old_e2.multiplicity) {
                                Ordering::Less => {
                                    // Reverse the new edge
                                    self.insert(
                                        kernel,
                                        new_e.reversed(),
                                        old_e1.hilbert_value,
                                        old_e2.multiplicity - old_e1.multiplicity,
                                    );
                                }
                                Ordering::Equal => {
                                    // Full cancellation, no edges to insert
                                }
                                Ordering::Greater => {
                                    self.insert(
                                        kernel,
                                        new_e.reversed(),
                                        old_e1.hilbert_value,
                                        old_e1.multiplicity - old_e2.multiplicity,
                                    );
                                }
                            }
                        }
                    }
                }
                Action::SplitEdge {
                    edge,
                    split_vertex: old_split_vertex,
                    pt,
                } => {
                    // Add split point as a new vertex
                    let new_split_vertex = kernel.new_split_vertex(pt);

                    // Update all edges referencing old_split_vertex
                    while let Some(edge) = { self.edges_for_vertex(old_split_vertex).next() } {
                        let new_edge = edge;
                        let old = self.remove(kernel, edge);
                        if let Some(new_edge) = kernel.replace_vertex_in_edge(
                            new_edge,
                            old_split_vertex,
                            new_split_vertex,
                        ) {
                            self.insert(kernel, new_edge, old.hilbert_value, old.multiplicity);
                        }
                    }

                    // Split the edge
                    let (new_e1, new_e2) = kernel.split_edge(edge, new_split_vertex);
                    let old = self.remove(kernel, edge);
                    self.insert(kernel, new_e1, old.hilbert_value, old.multiplicity);
                    self.insert(kernel, new_e2, old.hilbert_value, old.multiplicity);
                }
                Action::SplitBothEdges {
                    e1: edge_a,
                    e2: edge_b,
                    pt,
                } => {
                    let vertex = kernel.new_intersection_vertex(pt);
                    let (edge_a1, edge_a2) = kernel.split_edge(edge_a, vertex);
                    let (edge_b1, edge_b2) = kernel.split_edge(edge_b, vertex);

                    let old_a = self.remove(kernel, edge_a);
                    self.insert(kernel, edge_a1, old_a.hilbert_value, old_a.multiplicity);
                    self.insert(kernel, edge_a2, old_a.hilbert_value, old_a.multiplicity);
                    let old_b = self.remove(kernel, edge_b);
                    self.insert(kernel, edge_b1, old_b.hilbert_value, old_b.multiplicity);
                    self.insert(kernel, edge_b2, old_b.hilbert_value, old_b.multiplicity);
                }
            }
            break; // This function only handles one action
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
    spatial_index.clean(kernel);
    spatial_index.extract_edges()
}

// Full clean: clean all edges (assumes all dirty)
pub fn clean<K: Kernel>(kernel: &mut K, edges: impl Iterator<Item = K::Edge>) -> Vec<K::Edge> {
    partial_clean(kernel, edges.map(|edge| (edge, DirtyFlag::Dirty)))
}

#[cfg(test)]
mod tests {
    use crate::kernel::line::{BasicKernelF32, BasicKernelF32WithCustomEpsilon};

    use super::*;

    #[test]
    fn test_simple_two_edges() {
        let mut kernel = BasicKernelF32 {
            points: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]],
        };

        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut kernel, edges.into_iter());
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_intersecting_edges() {
        let mut kernel = BasicKernelF32 {
            points: vec![[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]],
        };

        // Two edges that intersect
        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut kernel, edges.into_iter());
        // Should split into 4 edges
        assert_eq!(result.len(), 4);
        // Should create a new vertex at the intersection
        assert_eq!(kernel.points.len(), 5);
    }

    #[test]
    fn test_vertex_on_edge() {
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [0.5, 0.0], // On the first edge
            ],
        };

        let edges = [(0, 1)];

        let result = clean(&mut kernel, edges.into_iter());
        // Should not split since there's only one edge and the vertex isn't connected
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_coincident_vertices() {
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 0.0], // Coincident with vertex 0
                [1.0, 1.0],
            ],
        };

        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut kernel, edges.into_iter());
        // Should merge vertices 0 and 2
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_square_polygon() {
        let mut kernel = BasicKernelF32 {
            points: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        };

        let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];

        let result = clean(&mut kernel, edges.into_iter());
        // Square edges don't intersect, should remain 4 edges
        assert_eq!(result.len(), 4);
        // No new vertices should be created
        assert_eq!(kernel.points.len(), 4);
    }

    #[test]
    fn test_t_junction() {
        let mut kernel = BasicKernelF32 {
            points: vec![[0.0, 0.5], [1.0, 0.5], [0.5, 0.0], [0.5, 1.0]],
        };

        // T-junction: horizontal edge intersected by vertical edge
        let edges = [(0, 1), (2, 3)];

        let result = clean(&mut kernel, edges.into_iter());
        // Should split into 4 edges meeting at center
        assert_eq!(result.len(), 4);
        // Should create intersection vertex at (0.5, 0.5)
        assert_eq!(kernel.points.len(), 5);
    }

    #[test]
    fn test_multiple_intersections() {
        let mut kernel = BasicKernelF32 {
            points: vec![
                // Horizontal line
                [0.0, 0.5],
                [1.0, 0.5],
                // Vertical line
                [0.5, 0.0],
                [0.5, 1.0],
                // Diagonal line
                [0.0, 0.0],
                [1.0, 1.0],
            ],
        };

        let edges = [(0, 1), (2, 3), (4, 5)];

        let result = clean(&mut kernel, edges.into_iter());
        // Each of the 3 edges should be split at their intersections
        // Horizontal and vertical intersect at [0.5, 0.5]
        // Diagonal intersects both at [0.5, 0.5]
        // So all three meet at one point, creating 6 edges
        assert_eq!(result.len(), 6);
    }

    #[test]
    fn test_complex_polygon() {
        // Star-like pattern
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.5, 0.0], // Top
                [0.7, 0.5], // Right
                [0.5, 1.0], // Bottom
                [0.3, 0.5], // Left
                [0.5, 0.5], // Center
            ],
        };

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
        let mut kernel = BasicKernelF32 {
            points: vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
        };

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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 0.0],
                [0.0, 0.0], // Coincident with vertex 0
                [1.0, 0.0],
            ],
        };

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
        let mut kernel = BasicKernelF32WithCustomEpsilon {
            points: vec![
                [0.0, 0.0],
                [2.0, 0.0],
                [2.0, 2.0],
                [1.0, 1.0],
                [3.5337768, 0.99069494],
                [3.118609, 2.8236988],
            ],
            epsilon: 0.1,
        };

        let edges = vec![(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)];

        let edges = clean(&mut kernel, edges.iter().copied());
        // Cleaning twice shouldn't change the result
        let edges2 = clean(&mut kernel, edges.iter().copied());
        assert_eq!(edges, edges2);
    }

    #[test]
    fn test_halting() {
        let mut kernel = BasicKernelF32WithCustomEpsilon {
            points: vec![
                [0.04894346, 3.309017],
                [1.0, 3.6618032],
                [0.8928008, 2.950759],
                [2.0, 3.4],
            ],
            epsilon: 0.7,
        };

        let edges = vec![(0, 1), (1, 3), (2, 3)];

        clean(&mut kernel, edges.iter().copied());
        // Clean shouldn't run forever
    }

    #[test]
    fn test_halting2() {
        let mut kernel = BasicKernelF32WithCustomEpsilon {
            points: vec![
                [1.0, 4.0],
                [1.2351141, 3.3236067],
                [1.9510565, 3.309017],
                [1.3804226, 2.8763933],
                [1.5877852, 2.190983],
                [1.0, 2.6],
                [0.41221482, 2.190983],
                [0.61957735, 2.8763933],
                [0.04894346, 3.309017],
                [0.76488584, 3.3236067],
                [-0.2, 1.7],
                [2.3, 1.7],
                [2.3, 4.2],
                [-0.2, 4.2],
                [1.5681425, 1.1981344],
                [3.0028763, 2.0102787],
                [2.543337, 3.3672395],
                [1.7937732, 1.6999999],
                [2.2999997, 2.8259892],
                [-0.19999997, 2.3000002],
                [1.7937732, 2.3],
                [0.40000004, 4.2],
                [0.40000004, 1.7],
                [0.036718935, 3.9088924],
                [0.75266135, 3.9234822],
                [-0.16207832, 2.3647285],
                [0.04528421, 3.0501387],
                [0.2570896, 2.3982692],
                [-0.31354427, 2.8308928],
                [0.19814801, 3.5206046],
                [0.43326217, 4.196998],
                [1.3427081, 2.107505],
                [0.75492287, 1.6984882],
                [1.5667379, 4.196998],
                [1.801852, 3.5206046],
                [1.2473387, 3.9234822],
                [1.963281, 3.9088924],
                [1.9547157, 3.0501387],
                [2.1620784, 2.3647285],
                [1.2725751, 1.7202837],
                [2.7073088, 2.532428],
                [1.2450773, 1.6984881],
                [0.657292, 2.107505],
                [2.3410113, 1.4539706],
                [2.1153808, 0.952105],
                [2.3135443, 2.8308928],
                [1.7429104, 2.3982692],
                [1.6999997, 2.8259892],
                [1.6999999, 4.2],
                [2.3, 3.6],
                [-0.2, 3.6],
                [3.0905752, 3.12121],
                [2.8472378, 2.5799599],
                [2.4345798, 1.8178233],
                [1.9750407, 3.174784],
                [-0.36487952, 2.8745625],
                [-0.41099325, 2.9237142],
                [-0.45130366, 2.9777272],
                [-0.48530215, 3.0359206],
                [-0.5125597, 3.0975595],
                [-0.5327325, 3.1618667],
                [-0.5455658, 3.2280307],
                [-0.55089784, 3.2952163],
                [-0.54866135, 3.362576],
                [-0.5388844, 3.4292603],
                [-0.5216905, 3.4944272],
                [-0.49729657, 3.5572546],
                [-0.46601033, 3.6169498],
                [-0.42822665, 3.6727595],
                [-0.38442215, 3.72398],
                [-0.3351496, 3.7699645],
                [-0.2810307, 3.8101327],
                [-0.2227484, 3.8439782],
                [-0.16103795, 3.8710737],
                [-0.09667799, 3.8910775],
                [-0.030480623, 3.903737],
                [0.69752705, 1.6631603],
                [0.63653123, 1.6344922],
                [0.5727051, 1.6128457],
                [0.50685394, 1.5984939],
                [0.4398086, 1.591618],
                [0.37241518, 1.5923045],
                [0.30552393, 1.600545],
                [0.23997883, 1.6162355],
                [0.17660691, 1.6391779],
                [0.11620784, 1.6690829],
                [0.05954376, 1.7055728],
                [0.0073294938, 1.7481877],
                [-0.039776117, 1.7963895],
                [-0.081178606, 1.84957],
                [-0.11635572, 1.9070585],
                [-0.14486349, 1.9681294],
                [-0.1663422, 2.0320122],
                [-0.18052095, 2.0979009],
                [-0.18722075, 2.164964],
                [-0.18635702, 2.2323554],
                [-0.17794085, 2.2992246],
                [0.45893115, 4.2593155],
                [0.49142712, 4.318361],
                [0.5303401, 4.3733892],
                [0.57517904, 4.4237065],
                [0.62537825, 4.4686775],
                [0.6803043, 4.507735],
                [0.73926413, 4.5403857],
                [0.8015138, 4.5662184],
                [0.866268, 4.5849066],
                [0.93270946, 4.596215],
                [1.0, 4.6],
                [1.0672905, 4.596215],
                [1.1337321, 4.5849066],
                [1.1984862, 4.5662184],
                [1.2607359, 4.5403857],
                [1.3196957, 4.507735],
                [1.3746219, 4.4686775],
                [1.424821, 4.4237065],
                [1.4696598, 4.3733892],
                [1.5085728, 4.318361],
                [1.5410688, 4.2593155],
                [2.177941, 2.2992249],
                [2.1863573, 2.2323554],
                [2.1872208, 2.164964],
                [2.180521, 2.0979009],
                [2.1663423, 2.0320122],
                [2.1448636, 1.9681294],
                [2.116356, 1.9070585],
                [2.081179, 1.84957],
                [2.0397763, 1.7963893],
                [1.9926707, 1.7481875],
                [1.9404564, 1.7055728],
                [1.8837922, 1.6690828],
                [1.8233931, 1.6391778],
                [1.7600213, 1.6162355],
                [1.6944762, 1.600545],
                [1.6275849, 1.5923045],
                [1.5601914, 1.5916178],
                [1.4931462, 1.5984938],
                [1.427295, 1.6128457],
                [1.3634688, 1.6344922],
                [1.3024731, 1.6631601],
                [1.860081, 2.2963247],
                [1.9255763, 2.2853441],
                [1.9894571, 2.2671928],
                [2.0509408, 2.242093],
                [2.1092737, 2.2103522],
                [2.1637416, 2.1723592],
                [2.2136774, 2.1285796],
                [2.2584689, 2.0795496],
                [2.2975676, 2.0258698],
                [2.3304946, 1.968198],
                [2.3568463, 1.9072405],
                [2.3763, 1.8437443],
                [2.3886175, 1.778487],
                [2.3936477, 1.7122682],
                [2.391329, 1.6458992],
                [2.38169, 1.5801929],
                [2.3648486, 1.5159544],
                [2.0304806, 3.903737],
                [2.0966778, 3.8910775],
                [2.161038, 3.8710737],
                [2.2227483, 3.8439782],
                [2.2810307, 3.8101327],
                [2.3351495, 3.7699645],
                [2.384422, 3.72398],
                [2.4282265, 3.6727595],
                [2.46601, 3.6169498],
                [2.4972963, 3.5572546],
                [2.5216904, 3.4944272],
                [2.5388842, 3.4292603],
                [2.5486612, 3.3625762],
                [2.5508976, 3.2952163],
                [2.5455656, 3.2280307],
                [2.5327322, 3.1618667],
                [2.5125597, 3.0975595],
                [2.485302, 3.0359206],
                [2.4513035, 2.9777274],
                [2.410993, 2.9237142],
                [2.3648794, 2.8745627],
                [2.815915, 2.5196702],
                [2.7779772, 2.4633083],
                [2.7339108, 2.4115968],
                [2.6842806, 2.3651984],
                [2.629723, 2.3247085],
                [2.570938, 2.290646],
                [2.508679, 2.2634478],
                [2.443744, 2.2434623],
                [2.376966, 2.2309463],
                [2.3092012, 2.22606],
                [2.2413185, 2.2288656],
                [2.1741881, 2.239328],
                [2.108671, 2.2573125],
                [2.0456069, 2.2825885],
                [1.9858048, 2.3148322],
                [1.9300312, 2.35363],
                [1.8790014, 2.3984842],
                [1.8333697, 2.44882],
                [1.7937213, 2.503992],
                [1.7605643, 2.5632927],
                [1.734324, 2.6259618],
                [1.7153368, 2.6911955],
                [1.7038463, 2.7581575],
            ],
            epsilon: 1e-5,
        };

        let edges = vec![
            (19, 20),
            (21, 22),
            (23, 24),
            (25, 26),
            (27, 28),
            (29, 30),
            (31, 32),
            (33, 34),
            (35, 36),
            (37, 38),
            (39, 40),
            (41, 42),
            (43, 44),
            (45, 46),
            (47, 48),
            (49, 50),
            (51, 52),
            (53, 54),
            (22, 10),
            (10, 19),
            (50, 13),
            (13, 21),
            (28, 55),
            (55, 56),
            (56, 57),
            (57, 58),
            (58, 59),
            (59, 60),
            (60, 61),
            (61, 62),
            (62, 63),
            (63, 64),
            (64, 65),
            (65, 66),
            (66, 67),
            (67, 68),
            (68, 69),
            (69, 70),
            (70, 71),
            (71, 72),
            (72, 73),
            (73, 74),
            (74, 75),
            (75, 23),
            (32, 76),
            (76, 77),
            (77, 78),
            (78, 79),
            (79, 80),
            (80, 81),
            (81, 82),
            (82, 83),
            (83, 84),
            (84, 85),
            (85, 86),
            (86, 87),
            (87, 88),
            (88, 89),
            (89, 90),
            (90, 91),
            (91, 92),
            (92, 93),
            (93, 94),
            (94, 95),
            (95, 96),
            (96, 25),
            (26, 7),
            (7, 27),
            (24, 9),
            (9, 29),
            (42, 5),
            (5, 31),
            (30, 97),
            (97, 98),
            (98, 99),
            (99, 100),
            (100, 101),
            (101, 102),
            (102, 103),
            (103, 104),
            (104, 105),
            (105, 106),
            (106, 107),
            (107, 108),
            (108, 109),
            (109, 110),
            (110, 111),
            (111, 112),
            (112, 113),
            (113, 114),
            (114, 115),
            (115, 116),
            (116, 117),
            (117, 33),
            (34, 1),
            (1, 35),
            (46, 3),
            (3, 37),
            (44, 14),
            (14, 39),
            (38, 118),
            (118, 119),
            (119, 120),
            (120, 121),
            (121, 122),
            (122, 123),
            (123, 124),
            (124, 125),
            (125, 126),
            (126, 127),
            (127, 128),
            (128, 129),
            (129, 130),
            (130, 131),
            (131, 132),
            (132, 133),
            (133, 134),
            (134, 135),
            (135, 136),
            (136, 137),
            (137, 138),
            (138, 41),
            (20, 139),
            (139, 140),
            (140, 141),
            (141, 142),
            (142, 143),
            (143, 144),
            (144, 145),
            (145, 146),
            (146, 147),
            (147, 148),
            (148, 149),
            (149, 150),
            (150, 151),
            (151, 152),
            (152, 153),
            (153, 154),
            (154, 155),
            (155, 43),
            (36, 156),
            (156, 157),
            (157, 158),
            (158, 159),
            (159, 160),
            (160, 161),
            (161, 162),
            (162, 163),
            (163, 164),
            (164, 165),
            (165, 166),
            (166, 167),
            (167, 168),
            (168, 169),
            (169, 170),
            (170, 171),
            (171, 172),
            (172, 173),
            (173, 174),
            (174, 175),
            (175, 176),
            (176, 45),
            (52, 177),
            (177, 178),
            (178, 179),
            (179, 180),
            (180, 181),
            (181, 182),
            (182, 183),
            (183, 184),
            (184, 185),
            (185, 186),
            (186, 187),
            (187, 188),
            (188, 189),
            (189, 190),
            (190, 191),
            (191, 192),
            (192, 193),
            (193, 194),
            (194, 195),
            (195, 196),
            (196, 197),
            (197, 198),
            (198, 199),
            (199, 47),
            (48, 12),
            (12, 49),
            (54, 16),
            (16, 51),
            (40, 15),
            (15, 53),
        ];
        clean(&mut kernel, edges.iter().copied());
        // Clean shouldn't run forever
    }
}

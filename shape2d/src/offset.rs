use std::collections::BTreeMap;

use crate::{
    ClippingError, clean, clip,
    kernel::{Edge, EdgeSide, Kernel, VertexEvent},
};

/// Error type for offset operations.
/// Passing in cleaned, clipped, manifold data to offset() should never cause an error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OffsetError {
    /// Invalid topology in the input geometry
    /// For example:
    /// * non-manifold geometry (open loops)
    /// * areas with winding number that is not 0 or 1
    Topology,
}

impl From<ClippingError> for OffsetError {
    fn from(value: ClippingError) -> OffsetError {
        match value {
            ClippingError::Topology => OffsetError::Topology,
        }
    }
}

pub fn offset_raw<K: Kernel>(
    kernel: &mut K,
    edges: impl Iterator<Item = K::Edge>,
    offset_amount: K::OffsetAmount,
    cap_style: K::CapStyle,
) -> Result<Vec<K::Edge>, OffsetError> {
    // Go through all the edges.
    // If the edge has endpoints, add an event for each endpoint to a list,
    // which will be used to establish topological connectivity.
    let mut events = Vec::new();
    let mut lone_edges = Vec::new();
    for edge in edges {
        if kernel.vertices_for_edge(edge).is_some() {
            events.push(VertexEvent::<K> {
                event_type: EdgeSide::Tail,
                edge,
            });
            events.push(VertexEvent::<K> {
                event_type: EdgeSide::Head,
                edge,
            });
        } else {
            lone_edges.push(edge);
        }
    }

    // Sort the events, which will group vertices together,
    // and order events around each vertex by incidence angle.
    events.sort_by(|a, b| kernel.vertex_event_cmp(a, b));

    // We will now pair up events of the form (Head, Tail)
    // to establish the connectivity of the data.
    // Each pair will be pushed to the edge_junctions list.
    let mut edge_junctions = BTreeMap::new();
    for vertex_events in events.chunk_by(|a, b| {
        a.event_type
            .select(kernel.vertices_for_edge(a.edge).unwrap())
            == b.event_type
                .select(kernel.vertices_for_edge(b.edge).unwrap())
    }) {
        match vertex_events[0].event_type {
            EdgeSide::Tail => {
                // Edges are sorted CCW around the vertex.
                // We want to find edges in [outgoing, incoming] pairs
                // since the area between these is filled.
                let (pairs, remainder) = vertex_events.as_chunks::<2>();
                if remainder.len() > 0 {
                    return Err(OffsetError::Topology);
                }
                for [outgoing, incoming] in pairs {
                    if !matches!(
                        (outgoing.event_type, incoming.event_type),
                        (EdgeSide::Tail, EdgeSide::Head)
                    ) {
                        return Err(OffsetError::Topology);
                    }
                    edge_junctions.insert(incoming.edge, outgoing.edge);
                }
            }
            EdgeSide::Head => {
                // If our first edge is incoming,
                // we skip it so that we can get the parity right,
                // and handle it at the end.
                let (pairs, remainder) = vertex_events[1..].as_chunks::<2>();
                if remainder.len() != 1 {
                    return Err(OffsetError::Topology);
                }
                for [outgoing, incoming] in pairs {
                    if !matches!(
                        (outgoing.event_type, incoming.event_type),
                        (EdgeSide::Tail, EdgeSide::Head)
                    ) {
                        return Err(OffsetError::Topology);
                    }
                    edge_junctions.insert(incoming.edge, outgoing.edge);
                }
                // Handle the last + first
                let [outgoing, incoming] = [&remainder[0], &vertex_events[0]];
                if !matches!(
                    (outgoing.event_type, incoming.event_type),
                    (EdgeSide::Tail, EdgeSide::Head)
                ) {
                    return Err(OffsetError::Topology);
                }
                edge_junctions.insert(incoming.edge, outgoing.edge);
            }
        }
    }

    let mut edge_loops = vec![];
    let mut i = 0_u32;

    for edge in lone_edges.into_iter() {
        edge_loops.push((i, edge));
        i += 1;
    }

    while let Some((start, mut cur)) = edge_junctions.pop_first() {
        edge_loops.push((i, start));
        while cur != start {
            edge_loops.push((i, cur));
            cur = edge_junctions.remove(&cur).ok_or(OffsetError::Topology)?;
        }
        i += 1;
    }

    // Process edge loops, offsetting each appropriately.
    Ok(kernel.offset_edge_loops(&edge_loops, offset_amount, cap_style))

    /*
    // Process edges with endpoints, offsetting each appropriately,
    // and tracking where they came from so they can be re-linked.
    let mut original_to_offset = BTreeMap::new();
    for edge in events.iter().filter_map(|event| match event.event_type {
        EdgeSide::Tail => Some(event.edge),
        EdgeSide::Head => None,
    }) {
        let offset_edges = kernel.offset_edge(
            edge,
            offset_amount,
            *prev_edge.get(&edge).unwrap(),
            *next_edge.get(&edge).unwrap(),
        );

        original_to_offset.insert(edge, offset_edge);
        result_edges.push(offset_edge);
    }

    // Process endpoints, adding caps as appropriate.
    //let mut fake_outgoing_edges = BTreeMap::new();
    //let mut fake_incoming_edges = BTreeMap::new();
    for (original_incoming_edge, original_outgoing_edge) in edge_junctions.into_iter() {
        let &offset_incoming_edge = original_to_offset.get(&original_incoming_edge).unwrap();
        let &offset_outgoing_edge = original_to_offset.get(&original_outgoing_edge).unwrap();

        let (_, original_vertex) = kernel.vertices_for_edge(original_incoming_edge).unwrap();
        if let Some(intersection) = kernel.intersection(offset_incoming_edge, offset_outgoing_edge)
        {
            let vertex = kernel.push_vertex(intersection);
            let (edge_a1, _) = kernel.split_edge(offset_incoming_edge, vertex);
            let (_, edge_b2) = kernel.split_edge(offset_outgoing_edge, vertex);

            result_edges.push(edge_a1);
            result_edges.push(edge_b2);
            original_to_offset.insert(original_incoming_edge, edge_a1);
            original_to_offset.insert(original_outgoing_edge, edge_b2);

            result_edges.retain(|&e| e != offset_incoming_edge && e != offset_outgoing_edge);
        } else
        {
            kernel.cap_corner(
                offset_incoming_edge,
                offset_outgoing_edge,
                original_vertex,
                offset_amount,
                cap_style,
                |e, _| {
                    result_edges.push(e);
                    //if fake {
                    //    fake_outgoing_edges.insert(offset_incoming_edge, e);
                    //    fake_incoming_edges.insert(offset_outgoing_edge, e);
                    //}
                },
            );
        }
    }
    let mut edges_to_remove = vec![];
    for (&edge, &incoming) in fake_incoming_edges.iter() {
        let Some(&outgoing) = fake_outgoing_edges.get(&edge) else {
            continue;
        };
        if kernel.intersection(incoming, outgoing).is_some() {
            edges_to_remove.push((incoming, edge, outgoing));
        }
    }
    for (incoming, edge, outgoing) in edges_to_remove {
        let (start, _) = kernel.vertices_for_edge(incoming).unwrap();
        let (_, end) = kernel.vertices_for_edge(outgoing).unwrap();

        result_edges.retain(|&e| e != incoming && e != outgoing && e != edge);

        let replacement_edge = kernel.line_between(start, end);
        result_edges.push(replacement_edge);

        let mut insertion = None;
        for (&next_edge, &next_incoming) in fake_incoming_edges.iter() {
            if next_incoming == outgoing {
                insertion = Some(next_edge);
                break;
            }
        }
        if let Some(next_edge) = insertion {
            fake_incoming_edges.insert(next_edge, replacement_edge);
        }
        let mut insertion = None;
        for (&prev_edge, &prev_outgoing) in fake_outgoing_edges.iter() {
            if prev_outgoing == incoming {
                insertion = Some(prev_edge);
            }
        }
        if let Some(prev_edge) = insertion {
            fake_outgoing_edges.insert(prev_edge, replacement_edge);
        }
    }
    */
}

pub fn offset<K: Kernel>(
    kernel: &mut K,
    edges: impl Iterator<Item = K::Edge>,
    offset_amount: K::OffsetAmount,
    cap_style: K::CapStyle,
) -> Result<Vec<K::Edge>, OffsetError> {
    let edges = offset_raw(kernel, edges, offset_amount, cap_style)?;
    let edges = clean(kernel, edges.into_iter());
    let edges = clip(kernel, edges.into_iter(), |w| w > 0)?;
    Ok(edges)
}

pub fn offset_segments_raw<K: Kernel>(
    kernel: &mut K,
    edges: impl Iterator<Item = K::Edge>,
    offset_amount: K::OffsetAmount,
    cap_style: K::CapStyle,
) -> Result<Vec<K::Edge>, OffsetError> {
    offset_raw(
        kernel,
        edges.flat_map(|e| [e, e.reversed()]),
        offset_amount,
        cap_style,
    )
}

pub fn offset_segments<K: Kernel>(
    kernel: &mut K,
    edges: impl Iterator<Item = K::Edge>,
    offset_amount: K::OffsetAmount,
    cap_style: K::CapStyle,
) -> Result<Vec<K::Edge>, OffsetError> {
    let edges = offset_segments_raw(kernel, edges, offset_amount, cap_style)?;
    let edges = clean(kernel, edges.into_iter());
    let edges = clip(kernel, edges.into_iter(), |w| w > 0)?;
    Ok(edges)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::polyline::{CapStyleF32, F32 as Kernel};

    #[test]
    fn test_offset_square_erodes_to_nothing() {
        // Create a 2x2 square centered at origin
        let mut kernel =
            Kernel::new_with_vertices(vec![[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]);

        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        // Erode by -2.0 (larger than the square's half-width)
        // This should completely erode the square
        let cap_style = CapStyleF32::Bevel;
        let result = offset(&mut kernel, edges.into_iter(), -2.0, cap_style).unwrap();

        // The result should be empty (square fully eroded)
        assert!(result.is_empty(), "Square should be completely eroded");
    }

    #[test]
    fn test_offset_square_erodes_partially() {
        // Create a 4x4 square centered at origin
        let mut kernel =
            Kernel::new_with_vertices(vec![[-2.0, -2.0], [2.0, -2.0], [2.0, 2.0], [-2.0, 2.0]]);

        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        // Erode by -0.5
        let cap_style = CapStyleF32::Arc { tolerance: 0.01 };
        let result = offset(&mut kernel, edges.into_iter(), -0.5, cap_style).unwrap();

        // The result should be a square
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_offset_square_expands() {
        // Create a 2x2 square centered at origin
        let mut kernel =
            Kernel::new_with_vertices(vec![[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]);

        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        // Expand by +0.5
        let cap_style = CapStyleF32::Bevel;
        let result = offset(&mut kernel, edges.into_iter(), 0.5, cap_style).unwrap();

        // Should have 4 edges + 4 caps
        assert_eq!(result.len(), 8,);
    }

    #[test]
    fn test_offset_triangle_expands() {
        // Create a triangle
        let mut kernel = Kernel::new_with_vertices(vec![[0.0, 2.0], [-2.0, -2.0], [2.0, -2.0]]);

        let edges = vec![(0, 1), (1, 2), (2, 0)];

        // Expand by +0.5
        let cap_style = CapStyleF32::Bevel;
        let result = offset(&mut kernel, edges.into_iter(), 0.5, cap_style).unwrap();

        // The result should have 3 edges + 3 caps
        assert!(!result.is_empty(), "Expanded triangle should have edges");
        assert_eq!(result.len(), 6);
    }

    #[test]
    fn test_offset_zero_amount() {
        // Create a simple square
        let mut kernel =
            Kernel::new_with_vertices(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);

        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        // Offset by 0 should return similar topology
        let cap_style = CapStyleF32::Bevel;
        let result = offset(&mut kernel, edges.into_iter(), 0.0, cap_style).unwrap();

        // Should have edges (though vertices may have changed)
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_offset_restores_approximately() {
        // Create a large square
        let mut kernel =
            Kernel::new_with_vertices(vec![[-5.0, -5.0], [5.0, -5.0], [5.0, 5.0], [-5.0, 5.0]]);

        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        let cap_style = CapStyleF32::Arc { tolerance: 0.01 };

        // Expand by +1.0
        let expanded = offset(&mut kernel, edges.into_iter(), 1.0, cap_style).unwrap();

        // Erode back by -1.1
        let restored = offset(&mut kernel, expanded.into_iter(), -1.1, cap_style).unwrap();

        // The result should have 4 edges (though not identical to original due to rounding)
        assert_eq!(restored.len(), 4);
    }

    #[test]
    fn test_conjoined_triangles() {
        // Create a large square
        let mut kernel =
            Kernel::new_with_vertices(vec![[0., 0.], [-1., -1.], [-1., -2.], [1., 2.], [-1., 1.]]);

        let edges = vec![(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)];

        // offset_raw should not return a topology error
        offset_raw(&mut kernel, edges.into_iter(), 0.1, CapStyleF32::Bevel).unwrap();
    }

    #[test]
    fn test_offset_segments() {
        let mut kernel = Kernel::new_with_vertices(vec![[0., 0.], [1., 1.]]);
        let edges = vec![(0, 1)];

        let edges =
            offset_segments_raw(&mut kernel, edges.into_iter(), 0.1, CapStyleF32::Bevel).unwrap();

        // Offset of a line with beveled edges should be a rectangle
        assert_eq!(edges.len(), 4);
    }
}

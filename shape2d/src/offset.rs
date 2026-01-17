use std::collections::BTreeMap;

use crate::{
    ClippingError, clean, clip,
    kernel::{EdgeSide, Kernel, VertexEvent},
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
    cap_style: &K::CapStyle,
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
    let mut edge_junctions = Vec::with_capacity(events.len() / 2);
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
                    edge_junctions.push((incoming.edge, outgoing.edge));
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
                    edge_junctions.push((incoming.edge, outgoing.edge));
                }
                // Handle the last + first
                let [outgoing, incoming] = [&remainder[0], &vertex_events[0]];
                if !matches!(
                    (outgoing.event_type, incoming.event_type),
                    (EdgeSide::Tail, EdgeSide::Head)
                ) {
                    return Err(OffsetError::Topology);
                }
                edge_junctions.push((incoming.edge, outgoing.edge));
            }
        }
    }

    let mut result_edges = Vec::new();

    // Process lone edges, offsetting each appropriately.
    for edge in lone_edges.into_iter() {
        result_edges.push(kernel.offset_edge(edge));
    }

    // Process edges with endpoints, offsetting each appropriately,
    // and tracking where they came from so they can be re-linked.
    let mut original_to_offset = BTreeMap::new();
    for edge in events.iter().filter_map(|event| match event.event_type {
        EdgeSide::Tail => Some(event.edge),
        EdgeSide::Head => None,
    }) {
        let offset_edge = kernel.offset_edge(edge);
        original_to_offset.insert(edge, offset_edge);
        result_edges.push(offset_edge);
    }

    // Process endpoints, adding caps as appropriate.
    for (original_incoming_edge, original_outgoing_edge) in edge_junctions.into_iter() {
        let &offset_incoming_edge = original_to_offset.get(&original_incoming_edge).unwrap();
        let &offset_outgoing_edge = original_to_offset.get(&original_outgoing_edge).unwrap();

        let (_, original_vertex) = kernel.vertices_for_edge(original_incoming_edge).unwrap();

        result_edges.push(kernel.cap_edges(
            offset_incoming_edge,
            offset_outgoing_edge,
            original_vertex,
            cap_style,
        ));
    }

    Ok(result_edges)
}

pub fn offset<K: Kernel>(
    kernel: &mut K,
    edges: impl Iterator<Item = K::Edge>,
    cap_style: &K::CapStyle,
) -> Result<Vec<K::Edge>, OffsetError> {
    let edges = offset_raw(kernel, edges, cap_style)?;
    let edges = clean(kernel, edges.into_iter());
    let edges = clip(kernel, edges.into_iter(), |w| w > 0)?;
    Ok(edges)
}

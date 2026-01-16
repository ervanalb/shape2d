use crate::kernel::Kernel;

/*
pub fn offset_raw<K: Kernel>(kernel: &mut K, edges: impl Iterator<Item = K::Edge>) -> Vec<K::Edge> {
    let mut events = Vec::new();
    for edge in edges {
        for event in geometry.sweep_line_events_for_edge(edge) {
            events.push(event);
        }
    }
    events.sort_by(|a, b| geometry.sweep_line_event_cmp_bottom_up(a, b));

    // Iterate over events for each vertex
    for vertex_events in events
        .chunk_by(|a, b| geometry.sweep_line_event_point(a) == geometry.sweep_line_event_point(b))
    {
        let pt = geometry.sweep_line_event_point(&vertex_events[0]);
        let vertex = geometry.sweep_line_event_point_to_triangle_vertex(triangle_kernel, pt);

        let end_events_count =
            vertex_events.partition_point(|a| matches!(a.event_type, SweepLineEventType::End));

        // Handle pairs of end events
        let (end_event_pairs, end_events_remainder) =
            vertex_events[0..end_events_count].as_chunks::<2>();

        let (start_event_pairs, start_events_remainder) =
            vertex_events[end_events_count..].as_chunks::<2>();

        if start_events_remainder.len() != end_events_remainder.len() {
            return Err(TriangulationError::Topology);
        }

        let leftover_event = if start_events_remainder.len() == 1 {
            let e1 = &end_events_remainder[0];
            let e2 = &start_events_remainder[0];

            let is_top = match (e1.segment.chain, e2.segment.chain) {
                (SweepLineChain::Top, SweepLineChain::Top) => true,

                (SweepLineChain::Bottom, SweepLineChain::Bottom) => false,
                _ => return Err(TriangulationError::Topology),
            };
            Some((e1, e2, is_top))
        } else {
            None
        };

        // Handle pairs of end events
        for [e1, e2] in end_event_pairs {
            match (e1.segment.chain, e2.segment.chain) {
                (SweepLineChain::Bottom, SweepLineChain::Top) => {
                    // End(Bottom) + End(Top) = end vertex
                    handle_end_vertex(
                        geometry,
                        triangle_kernel,
                        pt,
                        vertex,
                        e1,
                        e2,
                        &mut monotone_events,
                        &mut status,
                    )?;
                }
                (SweepLineChain::Top, SweepLineChain::Bottom) => {
                    // End(Top) + End(Bottom) = merge vertex
                    handle_merge_vertex(
                        geometry,
                        triangle_kernel,
                        pt,
                        vertex,
                        e1,
                        e2,
                        &mut monotone_events,
                        &mut status,
                    )?;
                }
                _ => {
                    return Err(TriangulationError::Topology);
                }
            }
        }

        // If the leftover pair is a top vertex,
        // handle it now, between end pairs and start pairs.
        if let Some((e1, e2, true)) = leftover_event {
            // End(Top) + Start(Top) = top vertex
            handle_top_vertex(
                geometry,
                triangle_kernel,
                pt,
                vertex,
                e1,
                e2,
                &mut monotone_events,
                &mut status,
            )?;
        }

        // Handle pairs of start events
        for [e1, e2] in start_event_pairs {
            match (e1.segment.chain, e2.segment.chain) {
                (SweepLineChain::Bottom, SweepLineChain::Top) => {
                    // Start(Top) + Start(Bottom) = start vertex
                    handle_start_vertex(
                        geometry,
                        triangle_kernel,
                        pt,
                        vertex,
                        e1,
                        e2,
                        &mut monotone_events,
                        &mut status,
                        &mut component_allocator,
                    );
                }
                (SweepLineChain::Top, SweepLineChain::Bottom) => {
                    // Start(Bottom) + Start(Top) = split vertex
                    handle_split_vertex(
                        geometry,
                        triangle_kernel,
                        pt,
                        vertex,
                        e1,
                        e2,
                        &mut monotone_events,
                        &mut status,
                        &mut component_allocator,
                    )?;
                }
                _ => {
                    return Err(TriangulationError::Topology);
                }
            }
        }

        // If the leftover pair is a bottom vertex,
        // handle it last
        if let Some((e1, e2, false)) = leftover_event {
            // End(Bottom) + Start(Bottom) = bottom vertex
            handle_bottom_vertex(
                geometry,
                triangle_kernel,
                pt,
                vertex,
                e1,
                e2,
                &mut monotone_events,
                &mut status,
            )?;
        }
    }

    Ok(monotone_events)
}
*/

use std::cmp::Ordering;

use crate::kernel::{Kernel, SweepLineChain, SweepLineEvent, SweepLineEventType, TriangleVertex};

/// Main triangulation function
///
/// Takes a geometry kernel and edges, and returns a triangle mesh
///
/// # Arguments
/// * `geometry` - The geometric kernel
/// * `edges` - Iterator over edges that bound the region to triangulate
///
/// # Returns
/// A tuple of (vertices, triangles) where:
/// - vertices is a Vec of G::TriangleVertex
/// - triangles is a Vec of [u32; 3] indices into the vertices
pub fn triangulate<G: Kernel>(
    geometry: &mut G,
    edges: impl Iterator<Item = G::Edge>,
) -> (Vec<G::TriangleVertex>, Vec<[u32; 3]>) {
    // Stage 1: Partition into monotone pieces
    let (vertices, monotone_events) = partition_into_monotone(geometry, edges);

    // Stage 2: Triangulate each monotone component
    let triangles = triangulate_monotone(vertices.as_slice(), monotone_events);

    (vertices, triangles)
}

/// Helper type to track vertex characterization
#[derive(Debug, Clone, Copy, PartialEq)]
enum HelperType {
    Top,
    Bottom,
    Merge {
        upper_component: MonotoneComponentIndex,
    },
}

type MonotoneComponentIndex = u32;

/// An entry in the status structure for stage 1 (partitioning into monotone components)
#[derive(Debug, Clone)]
struct StatusEntry<G: Kernel> {
    edge: G::Edge,
    segment: G::SweepLineEdgeSegment,
    chain: SweepLineChain,
    monotone_component: MonotoneComponentIndex,
    helper_vertex: u32,
    helper_type: HelperType,
}

/// Event for stage 2 (monotone triangulation)
#[derive(Debug, Clone, Copy)]
struct MonotoneEvent {
    vertex_index: u32,
    monotone_component: MonotoneComponentIndex,
    chain: SweepLineChain,
}

struct MonotoneComponentAllocator(MonotoneComponentIndex);

impl MonotoneComponentAllocator {
    fn new() -> Self {
        Self(0)
    }
    fn allocate(&mut self) -> MonotoneComponentIndex {
        let i = self.0;
        self.0 += 1;
        i
    }
}

struct IndexedList<T>(Vec<T>);

impl<T> IndexedList<T> {
    fn new() -> Self {
        Self(vec![])
    }
    // TODO: consider avoiding pushing a vertex more than once
    // (store a mapping and have a push_if_new() function?)
    // (or maybe pushing a vertex more than once is important
    // for preserving edge-specific metadata?)
    fn push(&mut self, value: T) -> u32 {
        let i = self.0.len() as u32;
        self.0.push(value);
        i
    }
}

impl<T> From<IndexedList<T>> for Vec<T> {
    fn from(value: IndexedList<T>) -> Self {
        value.0
    }
}

/// Stage 1: Partition into monotone pieces
fn partition_into_monotone<G: Kernel>(
    geometry: &mut G,
    edges: impl Iterator<Item = G::Edge>,
) -> (Vec<G::TriangleVertex>, Vec<MonotoneEvent>) {
    let mut vertices = IndexedList::new();
    let mut monotone_events = Vec::new();
    let mut status: Vec<StatusEntry<G>> = Vec::new();
    let mut monotone_component_allocator = MonotoneComponentAllocator::new();

    // Build event queue with events that share a vertex sorted clockwise
    let mut events = Vec::new();
    for edge in edges {
        for event in geometry.sweep_line_events_for_edge(edge) {
            events.push(event);
        }
    }
    events.sort_by(|a, b| geometry.sweep_line_event_cmp_clockwise(a, b));

    let (event_pairs, remainder) = events.as_chunks::<2>();
    if remainder.len() > 0 {
        panic!("Topology error--vertex with odd number of edges");
    }

    // Process events in pairs around each vertex
    for [e1, e2] in event_pairs {
        let v = geometry.sweep_line_event_point(&e1);
        let v2 = geometry.sweep_line_event_point(&e2);
        if v != v2 {
            panic!("Topology error--vertex with odd number of edges");
        }
        process_event_pair(
            geometry,
            v,
            e1,
            e2,
            &mut vertices,
            &mut monotone_events,
            &mut status,
            &mut monotone_component_allocator,
        );
    }

    (vertices.into(), monotone_events)
}

/// Process a pair of events at a vertex
fn process_event_pair<G: Kernel>(
    geometry: &mut G,
    vertex: G::SweepLineEventPoint,
    event_a: &SweepLineEvent<G>,
    event_b: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
    monotone_component_allocator: &mut MonotoneComponentAllocator,
) {
    // Categorize the vertex based on the event pair

    match (
        event_a.event_type,
        event_a.chain,
        event_b.event_type,
        event_b.chain,
    ) {
        (
            SweepLineEventType::Start,
            SweepLineChain::Top,
            SweepLineEventType::Start,
            SweepLineChain::Bottom,
        ) => {
            // Start(Top) + Start(Bottom) = start vertex
            handle_start_vertex(
                geometry,
                vertex,
                event_a,
                event_b,
                vertices,
                monotone_events,
                status,
                monotone_component_allocator,
            );
        }
        (
            SweepLineEventType::Start,
            SweepLineChain::Bottom,
            SweepLineEventType::Start,
            SweepLineChain::Top,
        ) => {
            // Start(Bottom) + Start(Top) = split vertex
            handle_split_vertex(
                geometry,
                vertex,
                event_a,
                event_b,
                vertices,
                monotone_events,
                status,
                monotone_component_allocator,
            );
        }
        (
            SweepLineEventType::End,
            SweepLineChain::Bottom,
            SweepLineEventType::End,
            SweepLineChain::Top,
        ) => {
            // End(Bottom) + End(Top) = end vertex
            handle_end_vertex(
                geometry,
                vertex,
                event_a,
                event_b,
                vertices,
                monotone_events,
                status,
            );
        }
        (
            SweepLineEventType::End,
            SweepLineChain::Top,
            SweepLineEventType::End,
            SweepLineChain::Bottom,
        ) => {
            // End(Top) + End(Bottom) = merge vertex
            handle_merge_vertex(
                geometry,
                vertex,
                event_a,
                event_b,
                vertices,
                monotone_events,
                status,
            );
        }
        (
            SweepLineEventType::End,
            SweepLineChain::Bottom,
            SweepLineEventType::Start,
            SweepLineChain::Bottom,
        ) => {
            // End(Bottom) + Start(Bottom) = bottom vertex
            handle_bottom_vertex(
                geometry,
                vertex,
                event_a,
                event_b,
                vertices,
                monotone_events,
                status,
            );
        }
        (
            SweepLineEventType::End,
            SweepLineChain::Top,
            SweepLineEventType::Start,
            SweepLineChain::Top,
        ) => {
            // End(Top) + Start(Top) = top vertex
            handle_top_vertex(
                geometry,
                vertex,
                event_a,
                event_b,
                vertices,
                monotone_events,
                status,
            );
        }
        _ => {
            panic!("Topology error--invalid arrangement of edges around vertex")
        }
    }
}

/// Handle start vertex
fn handle_start_vertex<G: Kernel>(
    geometry: &mut G,
    vertex: G::SweepLineEventPoint,
    _upper_event: &SweepLineEvent<G>,
    lower_event: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
    monotone_component_allocator: &mut MonotoneComponentAllocator,
) {
    let monotone_component = monotone_component_allocator.allocate();

    // Output vertex V with component I, chain Bottom
    let triangle_vertex = geometry.sweep_line_event_point_to_triangle_vertex(vertex);
    let vertex_index = vertices.push(triangle_vertex);
    monotone_events.push(MonotoneEvent {
        vertex_index,
        monotone_component,
        chain: SweepLineChain::Bottom,
    });

    // Insert lower segment into status
    let pos = find_insertion_point(geometry, status, vertex);

    status.insert(
        pos,
        StatusEntry {
            edge: lower_event.edge,
            segment: lower_event.segment,
            chain: lower_event.chain,
            monotone_component,
            helper_vertex: vertex_index,
            helper_type: HelperType::Top,
        },
    );
}

/// Handle end vertex
fn handle_end_vertex<G: Kernel>(
    geometry: &mut G,
    vertex: G::SweepLineEventPoint,
    lower_event: &SweepLineEvent<G>,
    upper_event: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
) {
    // Search status & remove lower segment, noting helper H & component index I
    let pos = find_status_entry(
        geometry,
        status,
        vertex,
        lower_event.edge,
        lower_event.segment,
    );
    let entry = status.remove(pos);
    let helper_type = entry.helper_type;
    let component_i = entry.monotone_component;

    // Output lower segment with index I, chain Bottom
    output_edge_segment(
        geometry,
        lower_event,
        vertices,
        monotone_events,
        component_i,
        SweepLineChain::Bottom,
    );

    // Output V with index I, chain Top (component end vertex)
    let vertex_index = vertices.push(geometry.sweep_line_event_point_to_triangle_vertex(vertex));
    monotone_events.push(MonotoneEvent {
        vertex_index,
        monotone_component: component_i,
        chain: SweepLineChain::Top,
    });

    // Check if helper was a merge vertex
    if let HelperType::Merge {
        upper_component: component_j,
    } = helper_type
    {
        // Output upper segment with index J, chain Top
        output_edge_segment(
            geometry,
            upper_event,
            vertices,
            monotone_events,
            component_j,
            SweepLineChain::Top,
        );

        // Output a second V with index J, chain Top (component end vertex)
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: component_j,
            chain: SweepLineChain::Top,
        });
    } else {
        // Output upper segment with index I, chain Top
        output_edge_segment(
            geometry,
            upper_event,
            vertices,
            monotone_events,
            component_i,
            SweepLineChain::Top,
        );
    }
}

/// Handle split vertex
fn handle_split_vertex<G: Kernel>(
    geometry: &mut G,
    vertex: G::SweepLineEventPoint,
    upper_event: &SweepLineEvent<G>,
    _lower_event: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
    monotone_component_allocator: &mut MonotoneComponentAllocator,
) {
    // Search status to find segment directly below this vertex, noting its index I and helper H
    let i = find_status_below(geometry, status, vertex);
    let segment_below = &mut status[i];
    let helper_vertex = segment_below.helper_vertex;
    let helper_type = segment_below.helper_type;
    let component_i = segment_below.monotone_component;

    let vertex_index = vertices.push(geometry.sweep_line_event_point_to_triangle_vertex(vertex));

    match helper_type {
        HelperType::Merge {
            upper_component: component_j,
        } => {
            // Output V with index J, chain Bottom
            monotone_events.push(MonotoneEvent {
                vertex_index,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            });

            // Output V with index I, chain Top
            monotone_events.push(MonotoneEvent {
                vertex_index,
                monotone_component: component_i,
                chain: SweepLineChain::Top,
            });

            // Edit the status to replace the existing helper H with V and set the helper type to Top
            segment_below.helper_vertex = vertex_index;
            segment_below.helper_type = HelperType::Top;

            // Insert upper segment into status with index J, helper V, helper type Bottom
            let pos = find_insertion_point(geometry, status, vertex);
            status.insert(
                pos,
                StatusEntry {
                    edge: upper_event.edge,
                    segment: upper_event.segment,
                    chain: upper_event.chain,
                    monotone_component: component_j,
                    helper_vertex: vertex_index,
                    helper_type: HelperType::Bottom,
                },
            );
        }
        HelperType::Bottom => {
            // Allocate new index J
            let component_j = monotone_component_allocator.allocate();

            // Output H with index J, chain Bottom (component start vertex)
            monotone_events.push(MonotoneEvent {
                vertex_index: helper_vertex,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            });

            // Output V with index J, chain Top
            monotone_events.push(MonotoneEvent {
                vertex_index,
                monotone_component: component_j,
                chain: SweepLineChain::Top,
            });

            // Output V with index I, chain Bottom
            monotone_events.push(MonotoneEvent {
                vertex_index,
                monotone_component: component_i,
                chain: SweepLineChain::Bottom,
            });

            // Edit status to replace existing index I with J
            // Edit status to replace existing helper H with V and set helper type to Top
            segment_below.monotone_component = component_j;
            segment_below.helper_vertex = vertex_index;
            segment_below.helper_type = HelperType::Top;

            // Insert upper segment into status with index I, helper V, helper type Bottom
            let pos = find_insertion_point(geometry, status, vertex);
            status.insert(
                pos,
                StatusEntry {
                    edge: upper_event.edge,
                    segment: upper_event.segment,
                    chain: upper_event.chain,
                    monotone_component: component_i,
                    helper_vertex: vertex_index,
                    helper_type: HelperType::Bottom,
                },
            );
        }
        HelperType::Top => {
            // Allocate new index J
            let component_j = monotone_component_allocator.allocate();

            // Output H with index J, chain Bottom (component start vertex)
            monotone_events.push(MonotoneEvent {
                vertex_index: helper_vertex,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            });

            // Output V with index J, chain Bottom
            monotone_events.push(MonotoneEvent {
                vertex_index,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            });

            // Output V with index I, chain Top
            monotone_events.push(MonotoneEvent {
                vertex_index,
                monotone_component: component_i,
                chain: SweepLineChain::Top,
            });

            // Edit status to replace existing helper H with V and set helper type to Top
            segment_below.helper_vertex = vertex_index;
            segment_below.helper_type = HelperType::Top;

            // Insert upper segment into status with index J, helper V, helper type Bottom
            let pos = find_insertion_point(geometry, status, vertex);
            status.insert(
                pos,
                StatusEntry {
                    edge: upper_event.edge,
                    segment: upper_event.segment,
                    chain: upper_event.chain,
                    monotone_component: component_j,
                    helper_vertex: vertex_index,
                    helper_type: HelperType::Bottom,
                },
            );
        }
    }
}

/// Handle merge vertex
fn handle_merge_vertex<G: Kernel>(
    geometry: &mut G,
    vertex: G::SweepLineEventPoint,
    lower_event: &SweepLineEvent<G>,
    upper_event: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
) {
    // Search status & remove upper segment, noting its index I2 and helper H2
    let pos = find_status_entry(
        geometry,
        status,
        vertex,
        upper_event.edge,
        upper_event.segment,
    );
    let upper_segment = status.remove(pos);
    let helper2_type = upper_segment.helper_type;
    let component_i2 = upper_segment.monotone_component;

    // Search status to find segment directly below this vertex, noting its index I and helper H
    let i = find_status_below(geometry, status, vertex);
    let segment_below = &mut status[i];
    let helper_h_type = segment_below.helper_type;
    let component_i = segment_below.monotone_component;

    let vertex_index = vertices.push(geometry.sweep_line_event_point_to_triangle_vertex(vertex));

    // Handle lower segment
    if let HelperType::Merge {
        upper_component: component_j,
    } = helper_h_type
    {
        // Output lower segment with index J, chain Top
        output_edge_segment(
            geometry,
            lower_event,
            vertices,
            monotone_events,
            component_j,
            SweepLineChain::Top,
        );

        // Output V with index J, chain Top (component end vertex)
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: component_j,
            chain: SweepLineChain::Top,
        });
    } else {
        // Output lower segment with index I, chain Top
        output_edge_segment(
            geometry,
            lower_event,
            vertices,
            monotone_events,
            component_i,
            SweepLineChain::Top,
        );
    }

    // (both cases) Output V with index I, chain Top
    monotone_events.push(MonotoneEvent {
        vertex_index,
        monotone_component: component_i,
        chain: SweepLineChain::Top,
    });

    // Output upper segment with index I2, chain Bottom
    output_edge_segment(
        geometry,
        upper_event,
        vertices,
        monotone_events,
        component_i2,
        SweepLineChain::Bottom,
    );

    if let HelperType::Merge {
        upper_component: component_j2,
    } = helper2_type
    {
        // Output V with index I2, chain Top (component end vertex)
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: component_i2,
            chain: SweepLineChain::Top,
        });

        // Output V with index J2, chain Bottom
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: component_j2,
            chain: SweepLineChain::Bottom,
        });

        // Edit status to set existing helper H to V (merge helper index = J2)
        segment_below.helper_vertex = vertex_index;
        segment_below.helper_type = HelperType::Merge {
            upper_component: component_j2,
        };
    } else {
        // Output V with index I2, chain Bottom
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: component_i2,
            chain: SweepLineChain::Bottom,
        });

        // Edit status to set existing helper H to V (merge helper index = I2)
        segment_below.helper_vertex = vertex_index;
        segment_below.helper_type = HelperType::Merge {
            upper_component: component_i2,
        };
    }
}

/// Handle bottom vertex
fn handle_bottom_vertex<G: Kernel>(
    geometry: &mut G,
    vertex: G::SweepLineEventPoint,
    left_event: &SweepLineEvent<G>,
    right_event: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
) {
    // Search status & remove left segment, noting its index I and helper H
    let event_point = geometry.sweep_line_event_point(left_event);
    let pos = find_status_entry(
        geometry,
        status,
        event_point,
        left_event.edge,
        left_event.segment,
    );
    let entry = status.remove(pos);
    let helper_type = entry.helper_type;
    let component = entry.monotone_component;

    // Output left segment with index I, chain Bottom
    output_edge_segment(
        geometry,
        left_event,
        vertices,
        monotone_events,
        component,
        SweepLineChain::Bottom,
    );

    let vertex_index = vertices.push(geometry.sweep_line_event_point_to_triangle_vertex(vertex));

    // Check if helper is a merge vertex
    if let HelperType::Merge { upper_component } = helper_type {
        // Output V with index I, chain Top (component end vertex)
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: component,
            chain: SweepLineChain::Top,
        });

        // Output V with index J, chain Bottom
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: upper_component,
            chain: SweepLineChain::Bottom,
        });

        // Insert right segment into status with index J and helper V (helper type: Bottom)
        let pos = find_insertion_point(geometry, status, event_point);
        status.insert(
            pos,
            StatusEntry {
                edge: right_event.edge,
                segment: right_event.segment,
                chain: right_event.chain,
                monotone_component: upper_component,
                helper_vertex: vertex_index,
                helper_type: HelperType::Bottom,
            },
        );
    } else {
        // Output V with index I, chain Bottom
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: component,
            chain: SweepLineChain::Bottom,
        });

        // Insert right segment into status with index I and helper V (helper type: Bottom)
        let pos = find_insertion_point(geometry, status, event_point);
        status.insert(
            pos,
            StatusEntry {
                edge: right_event.edge,
                segment: right_event.segment,
                chain: right_event.chain,
                monotone_component: component,
                helper_vertex: vertex_index,
                helper_type: HelperType::Bottom,
            },
        );
    }
}

/// Handle top vertex
fn handle_top_vertex<G: Kernel>(
    geometry: &mut G,
    vertex: G::SweepLineEventPoint,
    left_event: &SweepLineEvent<G>,
    _right_event: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
) {
    // Search status to find segment directly below this vertex, noting its index I and helper H
    let below_pos = find_status_below(geometry, status, vertex);
    let helper_type = status[below_pos].helper_type;
    let component = status[below_pos].monotone_component;

    let vertex_index = vertices.push(geometry.sweep_line_event_point_to_triangle_vertex(vertex));

    // Check if helper is a merge vertex
    if let HelperType::Merge { upper_component } = helper_type {
        // Output left segment with index J, chain Top
        output_edge_segment(
            geometry,
            left_event,
            vertices,
            monotone_events,
            upper_component,
            SweepLineChain::Top,
        );

        // Output V with index J, chain Top (component end vertex)
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: upper_component,
            chain: SweepLineChain::Top,
        });
    } else {
        // Output left segment with index I, chain Top
        output_edge_segment(
            geometry,
            left_event,
            vertices,
            monotone_events,
            component,
            SweepLineChain::Top,
        );
    }

    // (both cases) Output V with index I, chain Top
    monotone_events.push(MonotoneEvent {
        vertex_index,
        monotone_component: component,
        chain: SweepLineChain::Top,
    });

    // Edit status to set existing helper H to V (helper type: Top)
    status[below_pos].helper_vertex = vertex_index;
    status[below_pos].helper_type = HelperType::Top;
}

/// Find insertion point in status for a new segment
fn find_insertion_point<G: Kernel>(
    geometry: &G,
    status: &[StatusEntry<G>],
    event_point: G::SweepLineEventPoint,
) -> usize {
    status.partition_point(|entry| {
        let ord =
            geometry.sweep_line_segment_cmp(entry.edge, entry.segment, entry.chain, event_point);
        matches!(ord, Ordering::Less | Ordering::Equal)
    })
}

/// Find a specific status entry by edge and segment
fn find_status_entry<G: Kernel>(
    geometry: &G,
    status: &[StatusEntry<G>],
    event_point: G::SweepLineEventPoint,
    edge: G::Edge,
    segment: G::SweepLineEdgeSegment,
) -> usize {
    let pos = status.partition_point(|entry| {
        let ord =
            geometry.sweep_line_segment_cmp(entry.edge, entry.segment, entry.chain, event_point);
        matches!(ord, Ordering::Less)
    });

    // Linear search around pos to find matching edge and segment
    for i in pos..status.len() {
        if status[i].edge == edge && status[i].segment == segment {
            return i;
        }
    }
    for i in (0..pos).rev() {
        if status[i].edge == edge && status[i].segment == segment {
            return i;
        }
    }

    panic!("Status entry not found");
}

/// Find the status entry directly below a vertex
fn find_status_below<G: Kernel>(
    geometry: &G,
    status: &[StatusEntry<G>],
    event_point: G::SweepLineEventPoint,
) -> usize {
    let pos = status.partition_point(|entry| {
        let ord =
            geometry.sweep_line_segment_cmp(entry.edge, entry.segment, entry.chain, event_point);
        matches!(ord, Ordering::Less)
    });

    if pos == 0 {
        panic!("No segment below vertex");
    }
    pos - 1
}

/// Output an edge segment by discretizing it into triangle vertices
fn output_edge_segment<G: Kernel>(
    geometry: &mut G,
    event: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    component: u32,
    chain: SweepLineChain,
) {
    for triangle_vertex in geometry.sweep_line_edge_segment_to_triangle_vertices(
        event.edge,
        event.segment,
        event.chain,
    ) {
        let vertex_index = vertices.push(triangle_vertex);
        monotone_events.push(MonotoneEvent {
            vertex_index,
            monotone_component: component,
            chain,
        });
    }
}

/// Stage 2: Triangulate monotone components
fn triangulate_monotone<T: TriangleVertex>(
    vertices: &[T],
    mut events: Vec<MonotoneEvent>,
) -> Vec<[u32; 3]> {
    let mut triangles = Vec::new();

    // Sort events by component, then by sweep-line order
    events.sort_by(|a, b| {
        a.monotone_component
            .cmp(&b.monotone_component)
            .then_with(|| {
                vertices[a.vertex_index as usize].sweep_line_cmp(&vertices[b.vertex_index as usize])
            })
    });

    let mut stack: Vec<MonotoneEvent> = Vec::new();

    // Process each monotone component
    for events in events.chunk_by(|a, b| a.monotone_component == b.monotone_component) {
        stack.clear();

        let mut events_iter = events.into_iter();

        // Push first two vertices onto stack
        stack.push(*events_iter.next().unwrap());
        stack.push(
            *events_iter
                .next()
                .expect("Component has fewer than three vertices"),
        );

        for &v in events {
            // Peek the top of the stack
            let &top = stack.last().unwrap();

            if v.chain != top.chain {
                // Different chains: pop all and form triangles
                loop {
                    // Pop a vertex A from the stack
                    let a = stack.pop().unwrap();

                    // Peek a second vertex B from the stack. If this fails, the stack is empty and we are done.
                    let Some(b) = stack.last() else {
                        break;
                    };

                    // We are guaranteed to be able to form all triangles
                    triangles.push(match v.chain {
                        SweepLineChain::Bottom => [v.vertex_index, a.vertex_index, b.vertex_index],
                        SweepLineChain::Top => [v.vertex_index, b.vertex_index, a.vertex_index],
                    });
                }
                // Stack is now empty
                // Push T back onto the stack, then push V to the stack
                stack.push(top);
                stack.push(v);
            } else {
                // Same chain: form triangles while winding is positive
                let a = loop {
                    // Pop a vertex A from the stack
                    // Peek a second vertex B from the stack. If this fails, the stack is empty and we are done.
                    // If V is on the Bottom chain, form triangle V - B - A
                    // If V is on the Top chain, form triangle V - A - B
                    // If this triangle has positive winding (using `sin_cmp`) then output it.
                    // If this triangle has zero or negative winding, then break out of the loop
                    // After looping, push A back onto the stack, then push V on to the stack

                    // Pop a vertex A from the stack
                    let a = stack.pop().unwrap();

                    // Peek a second vertex B from the stack. If this fails, the stack is empty and we are done.
                    let Some(b) = stack.last() else {
                        break a;
                    };

                    let triangle = match v.chain {
                        SweepLineChain::Bottom => [v.vertex_index, b.vertex_index, a.vertex_index],
                        SweepLineChain::Top => [v.vertex_index, a.vertex_index, b.vertex_index],
                    };

                    // We might not be able to form all triangles
                    match vertices[triangle[0] as usize].sin_cmp(
                        &vertices[triangle[1] as usize],
                        &vertices[triangle[2] as usize],
                    ) {
                        Ordering::Greater => {
                            // Triangle OK -- output it
                            triangles.push(triangle);
                        }
                        Ordering::Less | Ordering::Equal => {
                            // Triangle invalid -- stop iteration
                            break a;
                        }
                    }
                };
                stack.push(a);
                stack.push(v);
            }
        }

        if stack.len() != 2 {
            panic!("Finished triangulating with wrong number of vertices leftover");
        }
    }

    triangles
}

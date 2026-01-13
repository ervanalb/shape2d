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

impl<T: std::fmt::Debug> IndexedList<T> {
    fn new() -> Self {
        Self(vec![])
    }
    fn push(&mut self, value: T) -> u32 {
        let i = self.0.len() as u32;
        println!("Push vertex {value:?} as {i}");
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
            let pt = geometry.sweep_line_event_point(&event);
            events.push((pt, event));
        }
    }
    // TODO: could make use of already-computed event point in sweep_line_event_cmp
    events.sort_by(|(_, a), (_, b)| geometry.sweep_line_event_cmp_clockwise(a, b));

    // Iterate over events for each vertex
    for vertex_events in events.chunk_by(|(p1, _), (p2, _)| p1 == p2) {
        let (pt, _) = vertex_events[0];
        let vertex_index = vertices.push(geometry.sweep_line_event_point_to_triangle_vertex(pt));

        let (event_pairs, remainder) = vertex_events.as_chunks::<2>();
        if remainder.len() > 0 {
            panic!("Topology error--vertex with odd number of edges");
        }

        for [(_, e1), (_, e2)] in event_pairs {
            process_event_pair(
                geometry,
                pt,
                vertex_index,
                e1,
                e2,
                &mut vertices,
                &mut monotone_events,
                &mut status,
                &mut monotone_component_allocator,
            );
        }
    }

    (vertices.into(), monotone_events)
}

/// Process a pair of events at a vertex
fn process_event_pair<G: Kernel>(
    geometry: &mut G,
    vertex: G::SweepLineEventPoint,
    vertex_index: u32,
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
            println!("vertex {:?} is a start vertex", vertex);
            handle_start_vertex(
                geometry,
                vertex,
                vertex_index,
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
            println!("vertex {:?} is a split vertex", vertex);
            handle_split_vertex(
                geometry,
                vertex,
                vertex_index,
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
            println!("vertex {:?} is a end vertex", vertex);
            handle_end_vertex(
                geometry,
                vertex,
                vertex_index,
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
            println!("vertex {:?} is a merge vertex", vertex);
            handle_merge_vertex(
                geometry,
                vertex,
                vertex_index,
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
            println!("vertex {:?} is a bottom vertex", vertex);
            handle_bottom_vertex(
                geometry,
                vertex,
                vertex_index,
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
            println!("vertex {:?} is a top vertex", vertex);
            handle_top_vertex(
                geometry,
                vertex,
                vertex_index,
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
    vertex_index: u32,
    _upper_event: &SweepLineEvent<G>,
    lower_event: &SweepLineEvent<G>,
    _vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
    monotone_component_allocator: &mut MonotoneComponentAllocator,
) {
    let monotone_component = monotone_component_allocator.allocate();

    // Output vertex V with component I, chain Bottom
    monotone_events.push(dbg!(MonotoneEvent {
        vertex_index,
        monotone_component,
        chain: SweepLineChain::Bottom,
    }));

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
    vertex_index: u32,
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
    let segment_below = status.remove(pos);
    let helper_vertex = segment_below.helper_vertex;
    let helper_type = segment_below.helper_type;
    let component_i = segment_below.monotone_component;

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
    monotone_events.push(dbg!(MonotoneEvent {
        vertex_index,
        monotone_component: component_i,
        chain: SweepLineChain::Top,
    }));

    // Check if helper was a merge vertex
    if let HelperType::Merge {
        upper_component: component_j,
    } = helper_type
    {
        // Output H with index J, chain Bottom
        // (deferred)
        if helper_vertex != vertex_index {
            monotone_events.push(dbg!(MonotoneEvent {
                vertex_index: helper_vertex,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            }));
        }

        // Output V with index J, chain Top (component end vertex)
        monotone_events.push(dbg!(MonotoneEvent {
            vertex_index,
            monotone_component: component_j,
            chain: SweepLineChain::Top,
        }));

        // Output upper segment with index J, chain Top
        output_edge_segment(
            geometry,
            upper_event,
            vertices,
            monotone_events,
            component_j,
            SweepLineChain::Top,
        );
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
    vertex_index: u32,
    upper_event: &SweepLineEvent<G>,
    _lower_event: &SweepLineEvent<G>,
    _vertices: &mut IndexedList<G::TriangleVertex>,
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

    match helper_type {
        HelperType::Merge {
            upper_component: component_j,
        } => {
            if helper_vertex != vertex_index {
                // Output H with index J, chain Bottom
                // (deferred)
                monotone_events.push(dbg!(MonotoneEvent {
                    vertex_index: helper_vertex,
                    monotone_component: component_j,
                    chain: SweepLineChain::Bottom,
                }));
            }
            // Output V with index J, chain Bottom
            monotone_events.push(dbg!(MonotoneEvent {
                vertex_index,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            }));

            // XXX
            // Output V with index I, chain Top
            if helper_vertex != vertex_index {
                monotone_events.push(dbg!(MonotoneEvent {
                    vertex_index,
                    monotone_component: component_i,
                    chain: SweepLineChain::Top,
                }));
            }

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
            monotone_events.push(dbg!(MonotoneEvent {
                vertex_index: helper_vertex,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            }));

            if helper_vertex != vertex_index {
                // Output V with index J, chain Top
                monotone_events.push(dbg!(MonotoneEvent {
                    vertex_index,
                    monotone_component: component_j,
                    chain: SweepLineChain::Top,
                }));

                // Output V with index I, chain Bottom
                monotone_events.push(dbg!(MonotoneEvent {
                    vertex_index,
                    monotone_component: component_i,
                    chain: SweepLineChain::Bottom,
                }));
            }

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
            monotone_events.push(dbg!(MonotoneEvent {
                vertex_index: helper_vertex,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            }));

            if helper_vertex != vertex_index {
                // Output V with index J, chain Bottom
                monotone_events.push(dbg!(MonotoneEvent {
                    vertex_index,
                    monotone_component: component_j,
                    chain: SweepLineChain::Bottom,
                }));

                // Output V with index I, chain Top
                monotone_events.push(dbg!(MonotoneEvent {
                    vertex_index,
                    monotone_component: component_i,
                    chain: SweepLineChain::Top,
                }));
            }

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
    vertex_index: u32,
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
    let helper2_vertex = upper_segment.helper_vertex;
    let helper2_type = upper_segment.helper_type;
    let component_i2 = upper_segment.monotone_component;

    // Search status to find segment directly below this vertex, noting its index I and helper H
    let i = find_status_below(geometry, status, vertex);
    let segment_below = &mut status[i];
    let helper_h_vertex = segment_below.helper_vertex;
    let helper_h_type = segment_below.helper_type;
    let component_i = segment_below.monotone_component;

    // Handle lower segment
    println!("Handling segment below");
    if let HelperType::Merge {
        upper_component: component_j,
    } = helper_h_type
    {
        println!("Helper {helper_h_vertex:?} is a merge vertex");
        // Output H with index J, chain Bottom
        // (deferred)
        if helper_h_vertex != vertex_index {
            monotone_events.push(dbg!(MonotoneEvent {
                vertex_index: helper_h_vertex,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            }));
        }

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
        monotone_events.push(dbg!(MonotoneEvent {
            vertex_index,
            monotone_component: component_j,
            chain: SweepLineChain::Top,
        }));
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
    if helper_h_vertex != vertex_index {
        monotone_events.push(dbg!(MonotoneEvent {
            vertex_index,
            monotone_component: component_i,
            chain: SweepLineChain::Top,
        }));
    }

    // Output upper segment with index I2, chain Bottom
    output_edge_segment(
        geometry,
        upper_event,
        vertices,
        monotone_events,
        component_i2,
        SweepLineChain::Bottom,
    );

    println!("Handling upper segment");
    if let HelperType::Merge {
        upper_component: component_j2,
    } = helper2_type
    {
        println!("Helper {helper2_vertex:?} is a merge vertex");
        // These would have sorted in the other order
        debug_assert!(vertex_index != helper2_vertex);

        // Output H with index J2, chain Bottom
        // (deferred)
        monotone_events.push(dbg!(MonotoneEvent {
            vertex_index: helper2_vertex,
            monotone_component: component_j2,
            chain: SweepLineChain::Bottom,
        }));

        // Output V with index I2, chain Top (component end vertex)
        monotone_events.push(dbg!(MonotoneEvent {
            vertex_index,
            monotone_component: component_i2,
            chain: SweepLineChain::Top,
        }));

        // Defer outputting V with index J2, chain Bottom,
        // since it might be an end node, and we therefore might want to put it on chain Top

        // Edit status to set existing helper H to V (merge helper index = J2)
        segment_below.helper_vertex = vertex_index;
        segment_below.helper_type = HelperType::Merge {
            upper_component: component_j2,
        };
    } else {
        // Defer outputting V with index I2, chain Bottom,
        // since it might be an end node, and we therefore might want to put it on chain Top

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
    vertex_index: u32,
    left_event: &SweepLineEvent<G>,
    right_event: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
) {
    // Search status & remove left segment, noting its index I and helper H
    let pos = find_status_entry(
        geometry,
        status,
        vertex,
        left_event.edge,
        left_event.segment,
    );
    let entry = status.remove(pos);
    let helper_vertex = entry.helper_vertex;
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

    // Check if helper is a merge vertex
    if let HelperType::Merge {
        upper_component: component_j,
    } = helper_type
    {
        if helper_vertex != vertex_index {
            // Output H with index J, chain Bottom
            // (deferred)
            monotone_events.push(dbg!(MonotoneEvent {
                vertex_index: helper_vertex,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            }));

            // Output V with index I, chain Top (component end vertex)
            monotone_events.push(dbg!(MonotoneEvent {
                vertex_index,
                monotone_component: component,
                chain: SweepLineChain::Top,
            }));
        }

        // Output V with index J, chain Bottom
        monotone_events.push(dbg!(MonotoneEvent {
            vertex_index,
            monotone_component: component_j,
            chain: SweepLineChain::Bottom,
        }));

        // Insert right segment into status with index J and helper V (helper type: Bottom)
        let pos = find_insertion_point(geometry, status, vertex);
        status.insert(
            pos,
            StatusEntry {
                edge: right_event.edge,
                segment: right_event.segment,
                chain: right_event.chain,
                monotone_component: component_j,
                helper_vertex: vertex_index,
                helper_type: HelperType::Bottom,
            },
        );
    } else {
        // Output V with index I, chain Bottom
        monotone_events.push(dbg!(MonotoneEvent {
            vertex_index,
            monotone_component: component,
            chain: SweepLineChain::Bottom,
        }));

        // Insert right segment into status with index I and helper V (helper type: Bottom)
        let pos = find_insertion_point(geometry, status, vertex);
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
    vertex_index: u32,
    left_event: &SweepLineEvent<G>,
    _right_event: &SweepLineEvent<G>,
    vertices: &mut IndexedList<G::TriangleVertex>,
    monotone_events: &mut Vec<MonotoneEvent>,
    status: &mut Vec<StatusEntry<G>>,
) {
    // Search status to find segment directly below this vertex, noting its index I and helper H
    let below_pos = find_status_below(geometry, status, vertex);
    let helper_vertex = status[below_pos].helper_vertex;
    let helper_type = status[below_pos].helper_type;
    let component = status[below_pos].monotone_component;

    // Check if helper is a merge vertex
    if let HelperType::Merge {
        upper_component: component_j,
    } = helper_type
    {
        // Output H with index J, chain Bottom
        // (deferred)
        if helper_vertex != vertex_index {
            monotone_events.push(dbg!(MonotoneEvent {
                vertex_index: helper_vertex,
                monotone_component: component_j,
                chain: SweepLineChain::Bottom,
            }));
        }

        // Output V with index J, chain Top (component end vertex)
        monotone_events.push(dbg!(MonotoneEvent {
            vertex_index,
            monotone_component: component_j,
            chain: SweepLineChain::Top,
        }));

        // Output left segment with index J, chain Top
        output_edge_segment(
            geometry,
            left_event,
            vertices,
            monotone_events,
            component_j,
            SweepLineChain::Top,
        );
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
    monotone_events.push(dbg!(MonotoneEvent {
        vertex_index,
        monotone_component: component,
        chain: SweepLineChain::Top,
    }));

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
// TODO: factor this out, along with the similar code in clip.rs, into a sweep_line module
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
// TODO: factor this out, along with the similar code in clip.rs, into a sweep_line module
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
        monotone_events.push(dbg!(MonotoneEvent {
            vertex_index,
            monotone_component: component,
            chain,
        }));
    }
}

/// Stage 2: Triangulate monotone components
fn triangulate_monotone<T: TriangleVertex + std::fmt::Debug>(
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

    println!("Vertices are: {:?}", vertices);
    println!("Monotone components are: {:?}", events);

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

        for &v in events_iter {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::polyline::F32 as Kernel;

    /// Helper to verify triangle winding and that all triangles reference valid vertices
    fn verify_triangulation(vertices: &[[f32; 2]], triangles: &[[u32; 3]]) {
        for &[i0, i1, i2] in triangles {
            assert!((i0 as usize) < vertices.len());
            assert!((i1 as usize) < vertices.len());
            assert!((i2 as usize) < vertices.len());

            // Check that triangle has positive area (counter-clockwise winding)
            let v0 = vertices[i0 as usize];
            let v1 = vertices[i1 as usize];
            let v2 = vertices[i2 as usize];

            let cross = (v1[0] - v0[0]) * (v2[1] - v0[1]) - (v1[1] - v0[1]) * (v2[0] - v0[0]);
            assert!(
                cross > 0.0,
                "Triangle [{}, {}, {}] has negative or zero winding: {:?} {:?} {:?}",
                i0,
                i1,
                i2,
                v0,
                v1,
                v2
            );
        }
    }

    #[test]
    fn test_simple_triangle() {
        // Simple triangle - has start, top, and end vertices
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [2.0, 0.0], [1.0, 1.0]]);
        let edges = vec![(0, 1), (1, 2), (2, 0)];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 1, "Triangle should produce 1 triangle");
    }

    #[test]
    fn test_simple_quad() {
        // Simple convex quad - has start, top, bottom, and end vertices
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 2, "Quad should produce 2 triangles");
    }

    #[test]
    fn test_regular_polygons() {
        use std::f32::consts::TAU;

        for sides in [3, 4, 5, 6, 7, 8, 9, 10] {
            eprintln!("Testing {}-gon", sides);

            let verts = (0..sides)
                .map(|i| {
                    let angle = (i as f32) * TAU / (sides as f32);
                    [angle.cos(), angle.sin()]
                })
                .collect();
            let edges = (0..sides).map(|i| (i, (i + 1) % sides));

            let mut kernel = Kernel::new(verts);
            let (vertices, triangles) = triangulate(&mut kernel, edges);

            verify_triangulation(&vertices, &triangles);
            assert_eq!(triangles.len(), (sides - 2) as usize);
        }
    }

    #[test]
    fn test_top_and_bottom_vertices() {
        // 4--3
        // |   \
        // |    2--1
        // |       |
        // 5-------0
        //
        // Classifications:
        // 5 = start
        // 4 = top
        // 3 = top
        // 2 = top
        // 0 = bottom
        // 1 = end
        let mut kernel = Kernel::new(vec![
            [3.0, 0.0], // 0
            [3.0, 2.0], // 1
            [2.0, 2.0], // 2
            [1.0, 3.0], // 3
            [0.0, 3.0], // 4
            [0.0, 0.0], // 5
        ]);
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 4);
    }

    #[test]
    fn test_split_vertex() {
        //   /--2
        //  /  /
        // 3  1
        //  \  \
        //   \--0
        //
        // Classifications:
        // 3 = start
        // 1 = split
        // 0 = end
        // 2 = end
        let mut kernel = Kernel::new(vec![
            [2.0, 0.0], // 0
            [1., 1.],   // 1
            [2.0, 2.0], // 2
            [0.0, 1.0], // 3
        ]);
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 2,);
    }

    #[test]
    fn test_start_split_vertex() {
        //   3-2
        //   //
        //   4
        //   \\
        //   0-1
        //
        // Classifications:
        // 4 = start
        // 4 = split
        // 0 = bottom
        // 3 = top
        // 1 = end
        // 2 = end
        let mut kernel = Kernel::new(vec![
            [1.0, 0.0], // 0
            [2.0, 0.0], // 1
            [2.0, 2.0], // 2
            [1.0, 2.0], // 3
            [0.0, 1.0], // 4
        ]);
        let edges = vec![(0, 1), (1, 4), (4, 2), (2, 3), (3, 4), (4, 0)];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 2,);
    }

    #[test]
    fn test_merge_vertex() {
        // Polygon with a merge vertex (concave pointing right)
        // Creates a notch that merges monotone regions
        //
        // 2--\
        //  \  \
        //   3  1
        //  /  /
        // 0--/
        //
        // Classifications:
        // 0 = start
        // 2 = start
        // 3 = merge
        // 1 = end
        let mut kernel = Kernel::new(vec![
            [0.0, 0.0], // 0
            [2.0, 1.0], // 1
            [0.0, 2.0], // 2
            [1.0, 1.0], // 3
        ]);
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 2);
    }

    #[test]
    fn test_merge_end_vertex() {
        //   3-2
        //    \\
        //     4
        //    //
        //   0-1
        //
        // Classifications:
        // 0 = start
        // 3 = start
        // 1 = bottom
        // 2 = top
        // 4 = merge
        // 4 = end
        let mut kernel = Kernel::new(vec![
            [0.0, 0.0], // 0
            [1.0, 0.0], // 1
            [1.0, 2.0], // 2
            [0.0, 2.0], // 3
            [2.0, 1.0], // 4
        ]);
        let edges = vec![(0, 1), (1, 4), (4, 2), (2, 3), (3, 4), (4, 0)];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 2,);
    }

    #[test]
    fn test_split_then_merge() {
        // 4------3
        //  \    /
        //   5  2
        //  /    \
        // 0------1
        //
        // Classifications:
        // 0 = start
        // 4 = start
        // 5 = split
        // 2 = merge
        // 1 = end
        // 3 = end
        let mut kernel = Kernel::new(vec![
            [0.0, 0.0], // 0
            [3.0, 0.0], // 1
            [2.0, 1.0], // 2
            [3.0, 2.0], // 3
            [0.0, 2.0], // 4
            [1.0, 1.0], // 5
        ]);
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 4,);
    }

    #[test]
    fn test_split_merge() {
        // 3---2
        //  \ /
        //   4
        //  / \
        // 0---1
        //
        // Classifications:
        // 0 = start
        // 3 = start
        // 4 = merge
        // 4 = split
        // 1 = end
        // 2 = end
        let mut kernel = Kernel::new(vec![
            [0.0, 0.0], // 0
            [2.0, 0.0], // 1
            [2.0, 2.0], // 2
            [0.0, 2.0], // 3
            [1.0, 1.0], // 4
        ]);
        let edges = vec![(0, 1), (1, 4), (4, 2), (2, 3), (3, 4), (4, 0)];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 2,);
    }

    #[test]
    fn test_split_and_top_and_bottom() {
        // This tests splitting with a top helper & bottom helper
        //
        //     /-6
        //    / /
        //   7 5-----4
        //  /       /
        // 0-----1 3
        //        \ \
        //         \-2
        //
        // Classifications:
        // 0 = start
        // 7 = top
        // 5 = split (helper=7)
        // 1 = bottom
        // 6 = end
        // 3 = split (helper=1)
        // 2 = end
        // 4 = end
        let mut kernel = Kernel::new(vec![
            [0.0, 1.0], // 0
            [3.0, 1.0], // 1
            [5.0, 0.0], // 2
            [4.0, 1.0], // 3
            [5.0, 2.0], // 4
            [2.0, 2.0], // 5
            [3.0, 3.0], // 6
            [1.0, 2.0], // 7
        ]);
        let edges = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 0),
        ];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 6);
    }

    #[test]
    fn test_split_top_and_split_bottom() {
        // This tests splitting with a coincident top helper & coincident bottom helper
        //
        //    7-6
        //    //
        //   5-----4
        //  /     /
        // 0-----1
        //        \\
        //        2-3
        //
        // Classifications:
        // 0 = start
        // 7 = top
        // 5 = split (helper=7)
        // 1 = bottom
        // 6 = end
        // 3 = split (helper=1)
        // 2 = end
        // 4 = end
        let mut kernel = Kernel::new(vec![
            [0.0, 1.0], // 0
            [4.0, 1.0], // 1
            [5.0, 0.0], // 2
            [6.0, 0.0], // 3
            [5.0, 2.0], // 4
            [1.0, 2.0], // 5
            [3.0, 3.0], // 6
            [2.0, 3.0], // 7
        ]);
        let edges = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 1),
            (1, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 5),
            (5, 0),
        ];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 4);
    }

    #[test]
    fn test_combs() {
        // 6--\
        // |   5
        // |  /
        // | 4
        // |  \
        // |   3
        // |  /
        // | 2
        // |  \
        // |   1
        // 0--/
        let right_comb = [
            [0.0, 0.0],
            [2.0, 1.0],
            [1.0, 2.0],
            [2.0, 3.0],
            [1.0, 4.0],
            [2.0, 5.0],
            [0.0, 6.0],
        ];
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0)];

        let up_comb = right_comb.map(|[x, y]| [-y, x]);
        let left_comb = right_comb.map(|[x, y]| [-x, -y]);
        let down_comb = right_comb.map(|[x, y]| [y, -x]);

        for (desc, comb) in [
            ("right", right_comb),
            ("up", up_comb),
            ("left", left_comb),
            ("down", down_comb),
        ] {
            eprintln!("Testing comb: {}", desc);
            let mut kernel = Kernel::new(comb.to_vec());

            let (vertices, triangles) = triangulate(&mut kernel, edges.iter().copied());

            verify_triangulation(&vertices, &triangles);
            assert_eq!(triangles.len(), 5);
        }
    }

    #[test]
    fn test_degenerate_combs() {
        // 6
        // |\
        // | 5
        // |/
        // 4
        // |\
        // | 3
        // |/
        // 2
        // |\
        // | 1
        // |/
        // 0
        let right_comb = [
            [0.0, 0.0],
            [1.0, 1.0],
            [0.0, 2.0],
            [1.0, 3.0],
            [0.0, 4.0],
            [1.0, 5.0],
            [0.0, 6.0],
        ];
        let edges = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 4),
            (4, 2),
            (2, 0),
        ];

        let up_comb = right_comb.map(|[x, y]| [-y, x]);
        let left_comb = right_comb.map(|[x, y]| [-x, -y]);
        let down_comb = right_comb.map(|[x, y]| [y, -x]);

        for (desc, comb) in [
            ("right", right_comb),
            ("up", up_comb),
            ("left", left_comb),
            ("down", down_comb),
        ] {
            eprintln!("Testing degenerate comb: {}", desc);
            let mut kernel = Kernel::new(comb.to_vec());

            let (vertices, triangles) = triangulate(&mut kernel, edges.iter().copied());

            verify_triangulation(&vertices, &triangles);
            assert_eq!(triangles.len(), 3);
        }
    }

    #[test]
    fn test_square_with_holes() {
        // Square with three triangular holes (not touching edges)
        let mut kernel = Kernel::new(vec![
            // Outer square (counter-clockwise)
            [0.0, 0.0],   // 0
            [10.0, 0.0],  // 1
            [10.0, 10.0], // 2
            [0.0, 10.0],  // 3
            // Triangle hole 1
            [2.0, 2.0], // 4
            [3.0, 4.0], // 5
            [4.0, 2.0], // 6
            // Triangle hole 2
            [4.0, 7.0], // 7
            [5.0, 9.0], // 8
            [6.0, 7.0], // 9
            // Triangle hole 3
            [6.0, 2.0], // 10
            [7.0, 4.0], // 11
            [8.0, 2.0], // 12
        ]);

        let edges = vec![
            // Outer square
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            // Triangle hole 1
            (4, 5),
            (5, 6),
            (6, 4),
            // Triangle hole 2
            (7, 8),
            (8, 9),
            (9, 7),
            // Triangle hole 3
            (10, 11),
            (11, 12),
            (12, 10),
        ];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 17);
    }

    #[test]
    fn test_square_with_edge_holes() {
        // 6---5---4
        // |~~/ \~~|
        // |~B---A~|
        // |/|~~~|\|
        // 7 |~~~| 3
        // |\|~~~|/|
        // |~8---9~|
        // |~~\ /~~|
        // 0---1---2
        //
        // Outer boundary: 0-1-2-3-4-5-6-7
        // Four triangular holes connecting them
        let mut kernel = Kernel::new(vec![
            // Outer boundary
            [0.0, 0.0], // 0 - bottom-left corner
            [2.0, 0.0], // 1 - bottom edge midpoint
            [4.0, 0.0], // 2 - bottom-right corner
            [4.0, 2.0], // 3 - right edge midpoint
            [4.0, 4.0], // 4 - top-right corner
            [2.0, 4.0], // 5 - top edge midpoint
            [0.0, 4.0], // 6 - top-left corner
            [0.0, 2.0], // 7 - left edge midpoint
            // Inner square
            [1.0, 1.0], // 8 - bottom of inner square
            [3.0, 1.0], // 9 - right of inner square
            [3.0, 3.0], // A - top of inner square
            [1.0, 3.0], // B - left of inner square
        ]);

        let edges = vec![
            // Outer perimeter (counter-clockwise)
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 0),
            // Triangle hole 1
            (1, 8),
            (8, 9),
            (9, 1),
            // Triangle hole 2
            (3, 9),
            (9, 10),
            (10, 3),
            // Triangle hole 3
            (5, 10),
            (10, 11),
            (11, 5),
            // Triangle hole 4
            (7, 11),
            (11, 8),
            (8, 7),
        ];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());

        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 10);
    }

    #[test]
    fn test_star_cross() {
        //    5-4
        //    \~/
        //  6\ v /3
        //  |~>8<~|
        //  7/ ^ \2
        //    /~\
        //    0-1

        let mut kernel = Kernel::new(vec![
            [1.0, 0.0], // 0
            [3.0, 0.0], // 1
            [4.0, 1.0], // 2
            [4.0, 3.0], // 3
            [3.0, 4.0], // 4
            [1.0, 4.0], // 5
            [0.0, 3.0], // 6
            [0.0, 1.0], // 7
            [2.0, 2.0], // 8
        ]);

        let edges = vec![
            // Bottom triangle
            (0, 1),
            (1, 8),
            (8, 0),
            // Right triangle
            (2, 3),
            (3, 8),
            (8, 2),
            // Top triangle
            (4, 5),
            (5, 8),
            (8, 4),
            // Left triangle
            (6, 7),
            (7, 8),
            (8, 6),
        ];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());
        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 4);
    }

    #[test]
    fn test_star_cross_2() {
        // Same as star_cross but rotated 45 degrees
        let mut kernel = Kernel::new(vec![
            [0.0, 1.0],
            [1.0, 0.0],
            [3.0, 0.0],
            [4.0, 1.0],
            [4.0, 3.0],
            [3.0, 4.0],
            [1.0, 4.0],
            [0.0, 3.0],
            [2.0, 2.0],
        ]);

        let edges = vec![
            // Bottom-left triangle
            (0, 1),
            (1, 8),
            (8, 0),
            // Bottom-right triangle
            (2, 3),
            (3, 8),
            (8, 2),
            // Top-right triangle
            (4, 5),
            (5, 8),
            (8, 4),
            // Top-left triangle
            (6, 7),
            (7, 8),
            (8, 6),
        ];

        let (vertices, triangles) = triangulate(&mut kernel, edges.into_iter());
        verify_triangulation(&vertices, &triangles);
        assert_eq!(triangles.len(), 4);
    }
}

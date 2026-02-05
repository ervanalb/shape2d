use std::cmp::Ordering;

use crate::kernel::Kernel;
use crate::sweep_line::{
    SweepLineChain, SweepLineEvent, SweepLineEventType, SweepLineStatus, SweepLineStatusEntry,
};
use crate::triangle_kernel::TriangleKernel;

/// Error type for triangulation operations.
/// Passing in cleaned data to triangulate() should never cause an error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TriangulationError {
    /// Invalid topology in the input geometry
    /// For example:
    /// * non-manifold geometry (open loops)
    /// * areas with winding number that is not 0 or 1
    Topology,
    /// Error during triangulation algorithm
    /// For example:
    /// * regions with zero area
    Triangulation,
}

/// Main triangulation function
///
/// Takes a geometry kernel, triangle kernel, and edges, and returns a triangle mesh
///
/// # Arguments
/// * `kernel` - The geometric kernel
/// * `triangle_kernel` - The triangle kernel for managing vertices
/// * `edges` - Iterator over edges that bound the region to triangulate
///
/// # Returns
/// A Result containing a vector of triangles
///
/// # Errors
/// Returns a TriangulationError if the input has invalid topology or
/// if the triangulation algorithm encounters an unexpected state
pub fn triangulate<K: Kernel>(
    kernel: &K,
    triangle_kernel: &mut K::TriangleKernel,
    edges: impl Iterator<Item = K::Edge>,
) -> Result<Vec<<K::TriangleKernel as TriangleKernel>::Triangle>, TriangulationError> {
    // Stage 1: Partition into monotone pieces
    let monotone_events = partition_into_monotone(kernel, triangle_kernel, edges)?;

    // Stage 2: Triangulate each monotone component
    let triangles = triangulate_monotone(triangle_kernel, monotone_events)?;

    Ok(triangles)
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

/// Algorithm-specific data for status entries in stage 1 (partitioning into monotone components)
#[derive(Debug, Clone)]
struct StatusData<V> {
    component: MonotoneComponentIndex,
    helper_vertex: V,
    helper_type: HelperType,
}

/// Event for stage 2 (monotone triangulation)
#[derive(Debug, Clone, Copy)]
struct MonotoneEvent<V> {
    vertex: V,
    component: MonotoneComponentIndex,
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

/// Stage 1: Partition into monotone pieces
fn partition_into_monotone<K: Kernel>(
    kernel: &K,
    triangle_kernel: &mut K::TriangleKernel,
    edges: impl Iterator<Item = K::Edge>,
) -> Result<Vec<MonotoneEvent<<K::TriangleKernel as TriangleKernel>::Vertex>>, TriangulationError> {
    let mut monotone_events = Vec::new();
    let mut status =
        SweepLineStatus::<K, StatusData<<K::TriangleKernel as TriangleKernel>::Vertex>>::new();
    let mut component_allocator = MonotoneComponentAllocator::new();

    // Build event queue with events that share a vertex
    let mut events = Vec::new();
    for edge in edges {
        for event in kernel.sweep_line_events_for_edge(edge) {
            events.push(event);
        }
    }
    events.sort_by(|a, b| kernel.sweep_line_event_cmp(a, b));

    // Iterate over events for each vertex
    for vertex_events in
        events.chunk_by(|a, b| kernel.sweep_line_event_point(a) == kernel.sweep_line_event_point(b))
    {
        let pt = kernel.sweep_line_event_point(&vertex_events[0]);
        let vertex = kernel.sweep_line_event_point_to_triangle_vertex(triangle_kernel, pt);

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
                        kernel,
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
                        kernel,
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
                kernel,
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
                        kernel,
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
                        kernel,
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
                kernel,
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

/// Handle start vertex
fn handle_start_vertex<K: Kernel>(
    kernel: &K,
    _triangle_kernel: &mut K::TriangleKernel,
    pt: K::SweepLineEventPoint,
    vertex: <K::TriangleKernel as TriangleKernel>::Vertex,
    lower_event: &SweepLineEvent<K>,
    _upper_event: &SweepLineEvent<K>,
    monotone_events: &mut Vec<MonotoneEvent<<K::TriangleKernel as TriangleKernel>::Vertex>>,
    status: &mut SweepLineStatus<K, StatusData<<K::TriangleKernel as TriangleKernel>::Vertex>>,
    component_allocator: &mut MonotoneComponentAllocator,
) {
    let component = component_allocator.allocate();

    // Output vertex V with component I, chain Bottom
    monotone_events.push(MonotoneEvent {
        vertex,
        component,
        chain: SweepLineChain::Bottom,
    });

    // Insert lower segment into status
    status.insert(
        kernel,
        pt,
        SweepLineStatusEntry::new(
            lower_event.segment,
            StatusData {
                component,
                helper_vertex: vertex,
                helper_type: HelperType::Top,
            },
        ),
    );
}

/// Handle end vertex
fn handle_end_vertex<K: Kernel>(
    kernel: &K,
    triangle_kernel: &mut K::TriangleKernel,
    pt: K::SweepLineEventPoint,
    vertex: <K::TriangleKernel as TriangleKernel>::Vertex,
    lower_event: &SweepLineEvent<K>,
    upper_event: &SweepLineEvent<K>,
    monotone_events: &mut Vec<MonotoneEvent<<K::TriangleKernel as TriangleKernel>::Vertex>>,
    status: &mut SweepLineStatus<K, StatusData<<K::TriangleKernel as TriangleKernel>::Vertex>>,
) -> Result<(), TriangulationError> {
    // Search status & remove lower segment, noting helper H & component index I
    let segment_below = status
        .remove(kernel, pt, &lower_event.segment)
        .ok_or(TriangulationError::Topology)?;

    // Output lower segment with index I, chain Bottom
    output_edge_segment(
        kernel,
        triangle_kernel,
        lower_event,
        monotone_events,
        segment_below.data.component,
        SweepLineChain::Bottom,
    );

    // Output V with index I, chain Top (component end vertex)
    monotone_events.push(MonotoneEvent {
        vertex,
        component: segment_below.data.component,
        chain: SweepLineChain::Top,
    });

    // Check if helper was a merge vertex
    if let HelperType::Merge { upper_component } = segment_below.data.helper_type {
        // Output H with index J, chain Bottom
        // (deferred)
        if segment_below.data.helper_vertex != vertex {
            monotone_events.push(MonotoneEvent {
                vertex: segment_below.data.helper_vertex,
                component: upper_component,
                chain: SweepLineChain::Bottom,
            });
        }

        // Output V with index J, chain Top (component end vertex)
        monotone_events.push(MonotoneEvent {
            vertex,
            component: upper_component,
            chain: SweepLineChain::Top,
        });

        // Output upper segment with index J, chain Top
        output_edge_segment(
            kernel,
            triangle_kernel,
            upper_event,
            monotone_events,
            upper_component,
            SweepLineChain::Top,
        );
    } else {
        // Output upper segment with index I, chain Top
        output_edge_segment(
            kernel,
            triangle_kernel,
            upper_event,
            monotone_events,
            segment_below.data.component,
            SweepLineChain::Top,
        );
    }

    Ok(())
}

/// Handle split vertex
fn handle_split_vertex<K: Kernel>(
    kernel: &K,
    _triangle_kernel: &mut K::TriangleKernel,
    pt: K::SweepLineEventPoint,
    vertex: <K::TriangleKernel as TriangleKernel>::Vertex,
    _lower_event: &SweepLineEvent<K>,
    upper_event: &SweepLineEvent<K>,
    monotone_events: &mut Vec<MonotoneEvent<<K::TriangleKernel as TriangleKernel>::Vertex>>,
    status: &mut SweepLineStatus<K, StatusData<<K::TriangleKernel as TriangleKernel>::Vertex>>,
    component_allocator: &mut MonotoneComponentAllocator,
) -> Result<(), TriangulationError> {
    // Search status to find segment directly below this vertex, noting its index I and helper H
    let segment_below = status
        .get_below_mut(kernel, pt)
        .ok_or(TriangulationError::Topology)?;

    match segment_below.data.helper_type {
        HelperType::Merge { upper_component } => {
            if segment_below.data.helper_vertex != vertex {
                // Output H with index J, chain Bottom
                // (deferred)
                monotone_events.push(MonotoneEvent {
                    vertex: segment_below.data.helper_vertex,
                    component: upper_component,
                    chain: SweepLineChain::Bottom,
                });
            }
            // Output V with index J, chain Bottom
            monotone_events.push(MonotoneEvent {
                vertex,
                component: upper_component,
                chain: SweepLineChain::Bottom,
            });

            // Output V with index I, chain Top
            if segment_below.data.helper_vertex != vertex {
                monotone_events.push(MonotoneEvent {
                    vertex,
                    component: segment_below.data.component,
                    chain: SweepLineChain::Top,
                });
            }

            // Edit the status to replace the existing helper H with V and set the helper type to Top
            segment_below.data.helper_vertex = vertex;
            segment_below.data.helper_type = HelperType::Top;

            // Insert upper segment into status with index J, helper V, helper type Bottom
            status.insert(
                kernel,
                pt,
                SweepLineStatusEntry::new(
                    upper_event.segment,
                    StatusData {
                        component: upper_component,
                        helper_vertex: vertex,
                        helper_type: HelperType::Bottom,
                    },
                ),
            );
        }
        HelperType::Bottom => {
            // Allocate new index J
            let new_component = component_allocator.allocate();

            // Output H with index J, chain Bottom (component start vertex)
            monotone_events.push(MonotoneEvent {
                vertex: segment_below.data.helper_vertex,
                component: new_component,
                chain: SweepLineChain::Bottom,
            });

            if segment_below.data.helper_vertex != vertex {
                // Output V with index J, chain Top
                monotone_events.push(MonotoneEvent {
                    vertex,
                    component: new_component,
                    chain: SweepLineChain::Top,
                });

                // Output V with index I, chain Bottom
                monotone_events.push(MonotoneEvent {
                    vertex,
                    component: segment_below.data.component,
                    chain: SweepLineChain::Bottom,
                });
            }

            // Edit status to replace existing index I with J
            // Edit status to replace existing helper H with V and set helper type to Top
            let segment_below_original_component = segment_below.data.component;
            segment_below.data.component = new_component;
            segment_below.data.helper_vertex = vertex;
            segment_below.data.helper_type = HelperType::Top;

            // Insert upper segment into status with index I, helper V, helper type Bottom
            status.insert(
                kernel,
                pt,
                SweepLineStatusEntry::new(
                    upper_event.segment,
                    StatusData {
                        component: segment_below_original_component,
                        helper_vertex: vertex,
                        helper_type: HelperType::Bottom,
                    },
                ),
            );
        }
        HelperType::Top => {
            // Allocate new index J
            let new_component = component_allocator.allocate();

            // Output H with index J, chain Bottom (component start vertex)
            monotone_events.push(MonotoneEvent {
                vertex: segment_below.data.helper_vertex,
                component: new_component,
                chain: SweepLineChain::Bottom,
            });

            if segment_below.data.helper_vertex != vertex {
                // Output V with index J, chain Bottom
                monotone_events.push(MonotoneEvent {
                    vertex,
                    component: new_component,
                    chain: SweepLineChain::Bottom,
                });

                // Output V with index I, chain Top
                monotone_events.push(MonotoneEvent {
                    vertex,
                    component: segment_below.data.component,
                    chain: SweepLineChain::Top,
                });
            }

            // Edit status to replace existing helper H with V and set helper type to Top
            segment_below.data.helper_vertex = vertex;
            segment_below.data.helper_type = HelperType::Top;

            // Insert upper segment into status with index J, helper V, helper type Bottom
            status.insert(
                kernel,
                pt,
                SweepLineStatusEntry::new(
                    upper_event.segment,
                    StatusData {
                        component: new_component,
                        helper_vertex: vertex,
                        helper_type: HelperType::Bottom,
                    },
                ),
            );
        }
    }

    Ok(())
}

/// Handle merge vertex
fn handle_merge_vertex<K: Kernel>(
    kernel: &K,
    triangle_kernel: &mut K::TriangleKernel,
    pt: K::SweepLineEventPoint,
    vertex: <K::TriangleKernel as TriangleKernel>::Vertex,
    lower_event: &SweepLineEvent<K>,
    upper_event: &SweepLineEvent<K>,
    monotone_events: &mut Vec<MonotoneEvent<<K::TriangleKernel as TriangleKernel>::Vertex>>,
    status: &mut SweepLineStatus<K, StatusData<<K::TriangleKernel as TriangleKernel>::Vertex>>,
) -> Result<(), TriangulationError> {
    // Search status & remove upper segment, noting its index I2 and helper H2
    let (segment_below, upper_segment) = status
        .get_below_mut_and_remove(kernel, pt, &upper_event.segment)
        .ok_or(TriangulationError::Topology)?;
    let segment_below = segment_below.ok_or(TriangulationError::Topology)?;

    // Handle lower segment
    if let HelperType::Merge {
        upper_component: segment_below_helper_upper_component,
    } = segment_below.data.helper_type
    {
        // Output H with index J, chain Bottom
        // (deferred)
        if segment_below.data.helper_vertex != vertex {
            monotone_events.push(MonotoneEvent {
                vertex: segment_below.data.helper_vertex,
                component: segment_below_helper_upper_component,
                chain: SweepLineChain::Bottom,
            });
        }

        // Output lower segment with index J, chain Top
        output_edge_segment(
            kernel,
            triangle_kernel,
            lower_event,
            monotone_events,
            segment_below_helper_upper_component,
            SweepLineChain::Top,
        );

        // Output V with index J, chain Top (component end vertex)
        monotone_events.push(MonotoneEvent {
            vertex,
            component: segment_below_helper_upper_component,
            chain: SweepLineChain::Top,
        });
    } else {
        // Output lower segment with index I, chain Top
        output_edge_segment(
            kernel,
            triangle_kernel,
            lower_event,
            monotone_events,
            segment_below.data.component,
            SweepLineChain::Top,
        );
    }

    if segment_below.data.helper_vertex != vertex {
        // (both cases) Output V with index I, chain Top
        monotone_events.push(MonotoneEvent {
            vertex,
            component: segment_below.data.component,
            chain: SweepLineChain::Top,
        });
    }

    if let HelperType::Merge {
        upper_component: upper_segment_helper_upper_component,
    } = upper_segment.data.helper_type
    {
        debug_assert!(vertex != upper_segment.data.helper_vertex);

        // Output H with index J2, chain Bottom
        // (deferred)
        monotone_events.push(MonotoneEvent {
            vertex: upper_segment.data.helper_vertex,
            component: upper_segment_helper_upper_component,
            chain: SweepLineChain::Bottom,
        });

        // Output V with index I2, chain Top (component end vertex)
        monotone_events.push(MonotoneEvent {
            vertex,
            component: upper_segment.data.component,
            chain: SweepLineChain::Top,
        });

        // Defer outputting V with index J2, chain Bottom,
        // since it might be an end node, and we therefore might want to put it on chain Top

        // Edit status to set existing helper H to V (merge helper index = J2)
        segment_below.data.helper_vertex = vertex;
        segment_below.data.helper_type = HelperType::Merge {
            upper_component: upper_segment_helper_upper_component,
        };
    } else {
        // Defer outputting V with index I2, chain Bottom,
        // since it might be an end node, and we therefore might want to put it on chain Top

        // Edit status to set existing helper H to V (merge helper index = I2)
        segment_below.data.helper_vertex = vertex;
        segment_below.data.helper_type = HelperType::Merge {
            upper_component: upper_segment.data.component,
        };
    }

    // Output upper segment with index I2, chain Bottom
    output_edge_segment(
        kernel,
        triangle_kernel,
        upper_event,
        monotone_events,
        upper_segment.data.component,
        SweepLineChain::Bottom,
    );

    Ok(())
}

/// Handle bottom vertex
fn handle_bottom_vertex<K: Kernel>(
    kernel: &K,
    triangle_kernel: &mut K::TriangleKernel,
    pt: K::SweepLineEventPoint,
    vertex: <K::TriangleKernel as TriangleKernel>::Vertex,
    left_event: &SweepLineEvent<K>,
    right_event: &SweepLineEvent<K>,
    monotone_events: &mut Vec<MonotoneEvent<<K::TriangleKernel as TriangleKernel>::Vertex>>,
    status: &mut SweepLineStatus<K, StatusData<<K::TriangleKernel as TriangleKernel>::Vertex>>,
) -> Result<(), TriangulationError> {
    // Search status & remove left segment, noting its index I and helper H
    let left_segment = status
        .remove(kernel, pt, &left_event.segment)
        .ok_or(TriangulationError::Topology)?;

    // Output left segment with index I, chain Bottom
    output_edge_segment(
        kernel,
        triangle_kernel,
        left_event,
        monotone_events,
        left_segment.data.component,
        SweepLineChain::Bottom,
    );

    // Check if helper is a merge vertex
    if let HelperType::Merge { upper_component } = left_segment.data.helper_type {
        if left_segment.data.helper_vertex != vertex {
            // Output H with index J, chain Bottom
            // (deferred)
            monotone_events.push(MonotoneEvent {
                vertex: left_segment.data.helper_vertex,
                component: upper_component,
                chain: SweepLineChain::Bottom,
            });

            // Output V with index I, chain Top (component end vertex)
            monotone_events.push(MonotoneEvent {
                vertex,
                component: left_segment.data.component,
                chain: SweepLineChain::Top,
            });
        }

        // Output V with index J, chain Bottom
        monotone_events.push(MonotoneEvent {
            vertex,
            component: upper_component,
            chain: SweepLineChain::Bottom,
        });

        // Insert right segment into status with index J and helper V (helper type: Bottom)
        status.insert(
            kernel,
            pt,
            SweepLineStatusEntry::new(
                right_event.segment,
                StatusData {
                    component: upper_component,
                    helper_vertex: vertex,
                    helper_type: HelperType::Bottom,
                },
            ),
        );
    } else {
        // Output V with index I, chain Bottom
        monotone_events.push(MonotoneEvent {
            vertex,
            component: left_segment.data.component,
            chain: SweepLineChain::Bottom,
        });

        // Insert right segment into status with index I and helper V (helper type: Bottom)
        status.insert(
            kernel,
            pt,
            SweepLineStatusEntry::new(
                right_event.segment,
                StatusData {
                    component: left_segment.data.component,
                    helper_vertex: vertex,
                    helper_type: HelperType::Bottom,
                },
            ),
        );
    }

    Ok(())
}

/// Handle top vertex
fn handle_top_vertex<K: Kernel>(
    kernel: &K,
    triangle_kernel: &mut K::TriangleKernel,
    pt: K::SweepLineEventPoint,
    vertex: <K::TriangleKernel as TriangleKernel>::Vertex,
    left_event: &SweepLineEvent<K>,
    _right_event: &SweepLineEvent<K>,
    monotone_events: &mut Vec<MonotoneEvent<<K::TriangleKernel as TriangleKernel>::Vertex>>,
    status: &mut SweepLineStatus<K, StatusData<<K::TriangleKernel as TriangleKernel>::Vertex>>,
) -> Result<(), TriangulationError> {
    // Search status to find segment directly below this vertex, noting its index I and helper H
    let segment_below = status
        .get_below_mut(kernel, pt)
        .ok_or(TriangulationError::Topology)?;

    // Check if helper is a merge vertex
    if let HelperType::Merge { upper_component } = segment_below.data.helper_type {
        // Output H with index J, chain Bottom
        // (deferred)
        if segment_below.data.helper_vertex != vertex {
            monotone_events.push(MonotoneEvent {
                vertex: segment_below.data.helper_vertex,
                component: upper_component,
                chain: SweepLineChain::Bottom,
            });
        }

        // Output V with index J, chain Top (component end vertex)
        monotone_events.push(MonotoneEvent {
            vertex,
            component: upper_component,
            chain: SweepLineChain::Top,
        });

        // Output left segment with index J, chain Top
        output_edge_segment(
            kernel,
            triangle_kernel,
            left_event,
            monotone_events,
            upper_component,
            SweepLineChain::Top,
        );
    } else {
        // Output left segment with index I, chain Top
        output_edge_segment(
            kernel,
            triangle_kernel,
            left_event,
            monotone_events,
            segment_below.data.component,
            SweepLineChain::Top,
        );
    }

    // (both cases) Output V with index I, chain Top
    if segment_below.data.helper_vertex != vertex {
        monotone_events.push(MonotoneEvent {
            vertex,
            component: segment_below.data.component,
            chain: SweepLineChain::Top,
        });
    }

    // Edit status to set existing helper H to V (helper type: Top)
    segment_below.data.helper_vertex = vertex;
    segment_below.data.helper_type = HelperType::Top;

    Ok(())
}

/// Output an edge segment by discretizing it into triangle vertices
fn output_edge_segment<K: Kernel>(
    kernel: &K,
    triangle_kernel: &mut K::TriangleKernel,
    event: &SweepLineEvent<K>,
    monotone_events: &mut Vec<MonotoneEvent<<K::TriangleKernel as TriangleKernel>::Vertex>>,
    component: u32,
    chain: SweepLineChain,
) {
    for vertex in
        kernel.sweep_line_edge_segment_to_triangle_vertices(triangle_kernel, &event.segment)
    {
        monotone_events.push(MonotoneEvent {
            vertex,
            component,
            chain,
        });
    }
}

/// Stage 2: Triangulate monotone components
fn triangulate_monotone<TK: TriangleKernel>(
    triangle_kernel: &mut TK,
    mut events: Vec<MonotoneEvent<TK::Vertex>>,
) -> Result<Vec<TK::Triangle>, TriangulationError> {
    let mut triangles = Vec::new();

    // Sort events by component, then by sweep-line order
    events.sort_by(|a, b| {
        a.component
            .cmp(&b.component)
            .then_with(|| triangle_kernel.sweep_line_cmp(a.vertex, b.vertex))
    });

    let mut stack: Vec<MonotoneEvent<TK::Vertex>> = Vec::new();

    // Process each monotone component
    for events in events.chunk_by(|a, b| a.component == b.component) {
        stack.clear();

        let mut events_iter = events.into_iter();

        // Push first two vertices onto stack
        stack.push(
            *events_iter
                .next()
                .ok_or(TriangulationError::Triangulation)?,
        );
        stack.push(
            *events_iter
                .next()
                .ok_or(TriangulationError::Triangulation)?,
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
                    let (v0, v1, v2) = match v.chain {
                        SweepLineChain::Bottom => (v.vertex, a.vertex, b.vertex),
                        SweepLineChain::Top => (v.vertex, b.vertex, a.vertex),
                    };
                    triangles.push(triangle_kernel.push_triangle(v0, v1, v2));
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

                    let (v0, v1, v2) = match v.chain {
                        SweepLineChain::Bottom => (v.vertex, b.vertex, a.vertex),
                        SweepLineChain::Top => (v.vertex, a.vertex, b.vertex),
                    };

                    // We might not be able to form all triangles
                    match triangle_kernel.sin_cmp(v0, v1, v2) {
                        Ordering::Less => {
                            // Triangle OK -- output it
                            triangles.push(triangle_kernel.push_triangle(v0, v1, v2));
                        }
                        Ordering::Greater | Ordering::Equal => {
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
            return Err(TriangulationError::Triangulation);
        }
    }

    Ok(triangles)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::line::BasicKernelF32;
    use crate::triangle_kernel::TriangleKernelF32;

    /// Helper to verify triangle winding and that all triangles reference valid vertices
    fn verify_triangulation(triangle_kernel: &TriangleKernelF32, triangles: &[[u32; 3]]) {
        for &[i0, i1, i2] in triangles {
            let v0 = triangle_kernel.pt(i0);
            let v1 = triangle_kernel.pt(i1);
            let v2 = triangle_kernel.pt(i2);

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
        // Simple triangle - has start, top, and end points
        let mut kernel = BasicKernelF32 {
            points: vec![[0.0, 0.0], [2.0, 0.0], [1.0, 1.0]],
        };
        let edges = vec![(0, 1), (1, 2), (2, 0)];
        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
        assert_eq!(triangles.len(), 1, "Triangle should produce 1 triangle");
    }

    #[test]
    fn test_simple_quad() {
        // Simple convex quad - has start, top, bottom, and end points
        let mut kernel = BasicKernelF32 {
            points: vec![[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]],
        };
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
        assert_eq!(triangles.len(), 2, "Quad should produce 2 triangles");
    }

    #[test]
    fn test_regular_polygons() {
        use std::f32::consts::TAU;

        for sides in [3, 4, 5, 6, 7, 8, 9, 10] {
            eprintln!("Testing {}-gon", sides);

            let verts: Vec<_> = (0..sides)
                .map(|i| {
                    let angle = (i as f32) * TAU / (sides as f32);
                    [angle.cos(), angle.sin()]
                })
                .collect();
            let edges = (0..sides).map(|i| (i, (i + 1) % sides));

            let mut kernel = BasicKernelF32 { points: verts };
            let mut triangle_kernel = TriangleKernelF32::new();

            let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges).unwrap();

            verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [3.0, 0.0], // 0
                [3.0, 2.0], // 1
                [2.0, 2.0], // 2
                [1.0, 3.0], // 3
                [0.0, 3.0], // 4
                [0.0, 0.0], // 5
            ],
        };
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [2.0, 0.0], // 0
                [1., 1.],   // 1
                [2.0, 2.0], // 2
                [0.0, 1.0], // 3
            ],
        };
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [1.0, 0.0], // 0
                [2.0, 0.0], // 1
                [2.0, 2.0], // 2
                [1.0, 2.0], // 3
                [0.0, 1.0], // 4
            ],
        };
        let edges = vec![(0, 1), (1, 4), (4, 2), (2, 3), (3, 4), (4, 0)];

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 0.0], // 0
                [2.0, 1.0], // 1
                [0.0, 2.0], // 2
                [1.0, 1.0], // 3
            ],
        };
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 0.0], // 0
                [1.0, 0.0], // 1
                [1.0, 2.0], // 2
                [0.0, 2.0], // 3
                [2.0, 1.0], // 4
            ],
        };
        let edges = vec![(0, 1), (1, 4), (4, 2), (2, 3), (3, 4), (4, 0)];

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 0.0], // 0
                [3.0, 0.0], // 1
                [2.0, 1.0], // 2
                [3.0, 2.0], // 3
                [0.0, 2.0], // 4
                [1.0, 1.0], // 5
            ],
        };
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 0.0], // 0
                [2.0, 0.0], // 1
                [2.0, 2.0], // 2
                [0.0, 2.0], // 3
                [1.0, 1.0], // 4
            ],
        };
        let edges = vec![(0, 1), (1, 4), (4, 2), (2, 3), (3, 4), (4, 0)];

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 1.0], // 0
                [3.0, 1.0], // 1
                [5.0, 0.0], // 2
                [4.0, 1.0], // 3
                [5.0, 2.0], // 4
                [2.0, 2.0], // 5
                [3.0, 3.0], // 6
                [1.0, 2.0], // 7
            ],
        };
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

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 1.0], // 0
                [4.0, 1.0], // 1
                [5.0, 0.0], // 2
                [6.0, 0.0], // 3
                [5.0, 2.0], // 4
                [1.0, 2.0], // 5
                [3.0, 3.0], // 6
                [2.0, 3.0], // 7
            ],
        };
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

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
            let mut kernel = BasicKernelF32 {
                points: comb.to_vec(),
            };

            let mut triangle_kernel = TriangleKernelF32::new();

            let triangles =
                triangulate(&mut kernel, &mut triangle_kernel, edges.iter().copied()).unwrap();

            verify_triangulation(&triangle_kernel, &triangles);
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
            let mut kernel = BasicKernelF32 {
                points: comb.to_vec(),
            };

            let mut triangle_kernel = TriangleKernelF32::new();

            let triangles =
                triangulate(&mut kernel, &mut triangle_kernel, edges.iter().copied()).unwrap();

            verify_triangulation(&triangle_kernel, &triangles);
            assert_eq!(triangles.len(), 3);
        }
    }

    #[test]
    fn test_square_with_holes() {
        // Square with three triangular holes (not touching edges)
        let mut kernel = BasicKernelF32 {
            points: vec![
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
            ],
        };

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

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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
        let mut kernel = BasicKernelF32 {
            points: vec![
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
            ],
        };

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

        let mut triangle_kernel = TriangleKernelF32::new();

        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();

        verify_triangulation(&triangle_kernel, &triangles);
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

        let mut kernel = BasicKernelF32 {
            points: vec![
                [1.0, 0.0], // 0
                [3.0, 0.0], // 1
                [4.0, 1.0], // 2
                [4.0, 3.0], // 3
                [3.0, 4.0], // 4
                [1.0, 4.0], // 5
                [0.0, 3.0], // 6
                [0.0, 1.0], // 7
                [2.0, 2.0], // 8
            ],
        };

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

        let mut triangle_kernel = TriangleKernelF32::new();
        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();
        verify_triangulation(&triangle_kernel, &triangles);
        assert_eq!(triangles.len(), 4);
    }

    #[test]
    fn test_star_cross_2() {
        // Same as star_cross but rotated 45 degrees
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.0, 1.0],
                [1.0, 0.0],
                [3.0, 0.0],
                [4.0, 1.0],
                [4.0, 3.0],
                [3.0, 4.0],
                [1.0, 4.0],
                [0.0, 3.0],
                [2.0, 2.0],
            ],
        };

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

        let mut triangle_kernel = TriangleKernelF32::new();
        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();
        verify_triangulation(&triangle_kernel, &triangles);
        assert_eq!(triangles.len(), 4);
    }

    #[test]
    fn test_touching_triangles() {
        let mut kernel = BasicKernelF32 {
            points: vec![
                [1.0034493, 0.002921236],
                [2.0, 0.0],
                [2.0, 2.0],
                [1.9927379, 0.9950979],
                [3.0, 1.0],
                [3.0, 3.0],
            ],
        };
        let edges = vec![(0, 1), (1, 3), (2, 0), (3, 2), (3, 4), (4, 5), (5, 3)];

        let mut triangle_kernel = TriangleKernelF32::new();
        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();
        verify_triangulation(&triangle_kernel, &triangles);
        assert_eq!(triangles.len(), 3);

        // Make sure the point 2., 2. appears in the output triangulation
        let v = triangle_kernel
            .points
            .iter()
            .position(|&v| v == [2., 2.])
            .unwrap() as u32;
        assert!(
            triangles
                .iter()
                .filter(|&&[a, b, c]| a == v || b == v || c == v)
                .count()
                > 0
        );
    }

    #[test]
    fn test_touching_triangles_2() {
        let mut kernel = BasicKernelF32 {
            points: vec![
                [0.004632079, -0.001372324],
                [2.0, 0.0],
                [1.0135162, 2.9981308],
                [-1.0050312, 0.9939545],
                [0.22551377, 1.0972531],
                [1.6716268, 0.9979948],
                [0.34009355, 0.9959848],
            ],
        };
        let edges = vec![
            (0, 1),
            (1, 5),
            (2, 6),
            (3, 6),
            (4, 3),
            (5, 2),
            (6, 0),
            (6, 4),
        ];

        let mut triangle_kernel = TriangleKernelF32::new();
        let triangles = triangulate(&mut kernel, &mut triangle_kernel, edges.into_iter()).unwrap();
        verify_triangulation(&triangle_kernel, &triangles);
        assert_eq!(triangles.len(), 4);
    }
}

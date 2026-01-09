use crate::Geometry;
use std::cmp::Ordering;

/// Event type in the sweep line algorithm
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum EventType {
    End, // End comes before Start in sorting
    Start,
}

/// Edge direction relative to vertex
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    Incoming,
    Outgoing,
}

/// An event in the sweep line algorithm
#[derive(Debug, Clone)]
struct Event<G: Geometry> {
    /// The type of event (start or end)
    event_type: EventType,
    /// The vertex index where this event occurs
    vertex: G::Vertex,
    /// The other (non-event) vertex index
    other_vertex: G::Vertex,
    /// The direction of the edge relative to the event vertex
    direction: Direction,
}

/// An entry in the status structure
#[derive(Debug, Clone)]
struct StatusEntry {
    /// Left vertex index (in sweep line order)
    left_vertex: VertexIndex,
    /// Right vertex index (in sweep line order)
    right_vertex: VertexIndex,
    /// Winding number above this edge
    winding_above: i32,
}

/// State for edge cancellation (tracking 2-cycles and redundant edges)
#[derive(Debug, Clone)]
struct CancellationState {
    event_type: EventType,
    event_vertex: VertexIndex,
    other_vertex: VertexIndex,
    weight: i32,
}

/// Main clipping function
///
/// Takes clean, manifold geometry and applies a winding rule to produce output edges.
///
/// # Arguments
/// * `vertices` - The vertex list
/// * `edges` - Iterator over edges specified as (from_index, to_index)
/// * `winding_rule` - Function that takes a winding number and returns true if it's "inside"
///
/// # Returns
/// A list of edges as (from_index, to_index) wound positively around the "inside" area
pub fn clip<G: Geometry>(
    geometry: &mut G,
    edges: impl Iterator<Item = G::Edge>,
    winding_rule: impl Fn(i32) -> bool,
) -> Vec<G::Edge> {
    // Build event queue
    let mut events = build_event_queue(geometry, edges);

    // Sort events
    sort_events(geometry, &mut events);

    // Cancel opposing edges (2-cycles)
    //let events = cancel_opposing_edges(events);

    // Run sweep line algorithm
    sweep_line(vertices, events, winding_rule)
}

/// Build the event queue from edges
fn build_event_queue<V: Vertex>(
    vertices: &[V],
    edges: impl Iterator<Item = (VertexIndex, VertexIndex)>,
) -> Vec<Event> {
    let mut events = Vec::new();

    for (tail, head) in edges {
        let tail_v = &vertices[tail as usize];
        let head_v = &vertices[head as usize];

        // Sort vertices in sweep-line order
        let (start_vertex, end_vertex, start_direction, end_direction) =
            match tail_v.sweep_line_cmp(head_v) {
                Ordering::Less => {
                    // Tail comes first, so edge goes outgoing from tail
                    (tail, head, Direction::Outgoing, Direction::Incoming)
                }
                Ordering::Greater => {
                    // Head comes first, so edge goes incoming to tail
                    (head, tail, Direction::Incoming, Direction::Outgoing)
                }
                Ordering::Equal => {
                    // Reflex edge (both vertices are the same), skip it
                    continue;
                }
            };

        // Create start event
        events.push(Event {
            event_type: EventType::Start,
            vertex: start_vertex,
            other_vertex: end_vertex,
            direction: start_direction,
        });

        // Create end event
        events.push(Event {
            event_type: EventType::End,
            vertex: end_vertex,
            other_vertex: start_vertex,
            direction: end_direction,
        });
    }

    events
}

/// Sort events in sweep-line order
fn sort_events<V: Vertex>(vertices: &[V], events: &mut [Event]) {
    events.sort_by(|a, b| {
        let a_vertex = &vertices[a.vertex as usize];
        let b_vertex = &vertices[b.vertex as usize];

        // First, compare by vertex position
        a_vertex
            .sweep_line_cmp(b_vertex)
            // Then by event type (End before Start, thanks to Ord derivation)
            .then_with(|| a.event_type.cmp(&b.event_type))
            // Then by incidence angle, bottom-to-top
            .then_with(|| {
                let a_other = &vertices[a.other_vertex as usize];
                let b_other = &vertices[b.other_vertex as usize];
                match a.event_type {
                    EventType::End => a_vertex.sin_cmp(a_other, b_other),
                    EventType::Start => a_vertex.sin_cmp(b_other, a_other),
                }
            })
    });
}

/// Cancel opposing edges (2-cycles) from the event queue.
/// When two edges go between the same vertices in opposite directions,
/// they cancel each other out (net weight of zero).
fn cancel_opposing_edges(events: Vec<Event>) -> Vec<Event> {
    let mut result = Vec::new();
    let mut state = CancellationState {
        event_type: EventType::Start,
        event_vertex: 0,
        other_vertex: 0,
        weight: 0,
    };

    /// Emit events based on state weight
    fn emit_state(result: &mut Vec<Event>, state: &CancellationState) {
        if state.weight > 0 {
            // Emit incoming events
            for _ in 0..state.weight {
                result.push(Event {
                    event_type: state.event_type,
                    vertex: state.event_vertex,
                    other_vertex: state.other_vertex,
                    direction: Direction::Incoming,
                });
            }
        } else if state.weight < 0 {
            // Emit outgoing events
            for _ in 0..(-state.weight) {
                result.push(Event {
                    event_type: state.event_type,
                    vertex: state.event_vertex,
                    other_vertex: state.other_vertex,
                    direction: Direction::Outgoing,
                });
            }
        }
    }

    for event in events {
        let event_weight = match event.direction {
            Direction::Incoming => 1,
            Direction::Outgoing => -1,
        };

        // Check if this event matches the current state
        if event.event_type == state.event_type
            && event.vertex == state.event_vertex
            && event.other_vertex == state.other_vertex
        {
            // Same event signature, update weight
            state.weight += event_weight;
        } else {
            // Different event, emit previous state
            emit_state(&mut result, &state);

            // Start new state
            state = CancellationState {
                event_type: event.event_type,
                event_vertex: event.vertex,
                other_vertex: event.other_vertex,
                weight: event_weight,
            };
        }
    }

    // Emit final state
    emit_state(&mut result, &state);
    result
}

/// Run the sweep line algorithm
fn sweep_line<V: Vertex>(
    vertices: &[V],
    events: Vec<Event>,
    winding_rule: impl Fn(i32) -> bool,
) -> Vec<(VertexIndex, VertexIndex)> {
    let mut status: Vec<StatusEntry> = Vec::new();
    let mut output: Vec<(VertexIndex, VertexIndex)> = Vec::new();

    for event in events {
        match event.event_type {
            EventType::End => {
                // Find and remove the edge from status
                let event_vertex = &vertices[event.vertex as usize];

                let pos = status.partition_point(|entry| {
                    let left = &vertices[entry.left_vertex as usize];
                    let right = &vertices[entry.right_vertex as usize];

                    // We are looking for this point
                    // Greater -|- Equal -|- Less
                    //           ^
                    matches!(left.sin_cmp(right, event_vertex), Ordering::Greater)
                });

                // Linear search around pos to find the matching edge
                // Search forward first, then backward
                let mut found_at = None;

                // Search forward from pos
                for i in pos..status.len() {
                    let entry = &status[i];

                    if entry.left_vertex == event.other_vertex && entry.right_vertex == event.vertex
                    {
                        found_at = Some(i);
                        break;
                    }
                }

                // If not found forward, search backward from pos
                if found_at.is_none() {
                    for i in (0..pos).rev() {
                        let entry = &status[i];

                        if entry.left_vertex == event.other_vertex
                            && entry.right_vertex == event.vertex
                        {
                            found_at = Some(i);
                            break;
                        }
                    }
                }

                let found_at = found_at.expect("Previously inserted edge not found in status");

                if found_at != pos {
                    eprintln!("Warning: Edge removal required linear search (dirty input?)");
                }
                status.remove(found_at);
            }
            EventType::Start => {
                // Find insertion point in status
                let event_vertex_ref = &vertices[event.vertex as usize];

                let pos = status.partition_point(|entry| {
                    let left = &vertices[entry.left_vertex as usize];
                    let right = &vertices[entry.right_vertex as usize];

                    // We are looking for this point
                    // Greater -|- Equal -|- Less
                    //                     ^
                    matches!(
                        left.sin_cmp(right, event_vertex_ref),
                        Ordering::Greater | Ordering::Equal
                    )
                });

                // Calculate winding number
                let winding_below = if pos > 0 {
                    status[pos - 1].winding_above
                } else {
                    0
                };

                let winding_above = match event.direction {
                    Direction::Outgoing => winding_below + 1,
                    Direction::Incoming => winding_below - 1,
                };

                // Determine if we should emit this edge
                let below_inside = winding_rule(winding_below);
                let above_inside = winding_rule(winding_above);

                if below_inside != above_inside {
                    if above_inside {
                        // Bottom edge: emit in forward direction
                        output.push((event.vertex, event.other_vertex));
                    } else {
                        // Top edge: emit in reverse direction
                        output.push((event.other_vertex, event.vertex));
                    }
                }

                // Insert into status
                status.insert(
                    pos,
                    StatusEntry {
                        left_vertex: event.vertex,
                        right_vertex: event.other_vertex,
                        winding_above,
                    },
                );
            }
        }
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clean;

    #[test]
    fn test_simple_square() {
        // Simple square
        let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        // Positive winding rule
        let result = clip(&vertices, edges.iter().copied(), |w| w > 0);

        // All edges should be in the output
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_empty_input() {
        // No edges
        let vertices: Vec<[f32; 2]> = vec![[0.0, 0.0], [1.0, 0.0]];
        let edges: Vec<(u32, u32)> = vec![];

        let result = clip(&vertices, edges.iter().copied(), |w| w > 0);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_overlapping_squares() {
        // Two overlapping squares that should be merged
        let mut vertices = vec![
            // First square
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            // Second square (offset by 0.5)
            [0.5, 0.5],
            [1.5, 0.5],
            [1.5, 1.5],
            [0.5, 1.5],
        ];

        let edges = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
        ];

        // Clean the edges first to handle intersections
        let cleaned_edges = clean(&mut vertices, edges.iter().copied());
        assert!(cleaned_edges.len() == 12);

        let union = clip(&vertices, cleaned_edges.iter().copied(), |w| w > 0);
        // 4 edges lie on the interior and should have been clipped
        assert!(union.len() == 8);

        let intersection = clip(&vertices, cleaned_edges.iter().copied(), |w| w > 1);
        // Intersection is a square
        assert!(intersection.len() == 4);
    }

    #[test]
    fn test_square_with_hole() {
        // Outer square and inner square (hole)
        let vertices = vec![
            // Outer square
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            // Inner square (reversed winding)
            [0.5, 0.5],
            [0.5, 1.5],
            [1.5, 1.5],
            [1.5, 0.5],
        ];

        let edges = vec![
            // Outer square (counter-clockwise)
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            // Inner square (clockwise - opposite winding)
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
        ];

        // Positive winding rule shouldn't change geometry
        let result = clip(&vertices, edges.iter().copied(), |w| w > 0);
        assert_eq!(result.len(), 8);

        // Non-positive winding rule also shouldn't change geometry
        let result = clip(&vertices, edges.iter().copied(), |w| w <= 0);
        assert_eq!(result.len(), 8);
    }

    #[test]
    fn test_empty_output() {
        // Square with winding number that will be filtered out
        let vertices = vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        // Rule that filters out everything
        let result = clip(&vertices, edges.iter().copied(), |w| w > 1);
        assert_eq!(result.len(), 0);
    }
}

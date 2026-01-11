use std::{
    cmp::Ordering,
    collections::{BTreeMap, btree_map::Entry},
};

use crate::kernel::{Edge, Kernel, SweepLineChain, SweepLineEvent, SweepLineEventType};

/// Edge direction
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    Forward,
    Reverse,
}

/// An entry in the status structure
#[derive(Debug, Clone)]
struct StatusEntry<K: Kernel> {
    /// The active edge
    edge: K::Edge,
    /// The active segment
    segment: K::SweepLineEdgeSegment,
    /// Whether this is a top or bottom edge
    chain: SweepLineChain,
    /// The winding number above this segment
    winding_above: i32,
}

/// Main clipping function
///
/// Takes clean, manifold geometry and applies a winding rule to produce output edges.
///
/// # Arguments
/// * `geometry` - The geometric context
/// * `edges` - Iterator over edges
/// * `winding_rule` - Function that takes a winding number and returns true if it's "inside"
///
/// # Returns
/// A list of edges wound positively around the "inside" area
pub fn clip<K: Kernel>(
    geometry: &mut K,
    edges: impl Iterator<Item = K::Edge>,
    winding_rule: impl Fn(i32) -> bool,
) -> Vec<K::Edge> {
    // Build sorted event queue
    let events = build_event_queue(geometry, edges);

    // Run sweep line algorithm
    sweep_line(geometry, events, winding_rule)
}

/// Build the event queue from edges
fn build_event_queue<K: Kernel>(
    geometry: &mut K,
    edges: impl Iterator<Item = K::Edge>,
) -> Vec<SweepLineEvent<K>> {
    let mut events = Vec::new();

    for e in edges {
        for event in geometry.sweep_line_events_for_edge(e) {
            events.push(event);
        }
    }

    // Sort events
    events.sort_by(|a, b| geometry.sweep_line_event_cmp_bottom_up(a, b));
    events
}

/// Run the sweep line algorithm
fn sweep_line<K: Kernel>(
    geometry: &mut K,
    events: Vec<SweepLineEvent<K>>,
    winding_rule: impl Fn(i32) -> bool,
) -> Vec<K::Edge> {
    let mut status: Vec<StatusEntry<K>> = Vec::new();
    let mut output: BTreeMap<K::Edge, Option<Direction>> = BTreeMap::new();

    for event in events.iter() {
        match event.event_type {
            SweepLineEventType::End => {
                // Find and remove the edge from status
                let event_point = geometry.sweep_line_event_point(event);
                let pos = status.partition_point(|status_entry| {
                    let ord = geometry.sweep_line_segment_cmp(
                        status_entry.edge,
                        status_entry.segment,
                        status_entry.chain,
                        event_point,
                    );
                    // We are looking for this point
                    // Less -|- Equal -|- Greater
                    //        ^
                    matches!(ord, Ordering::Less)
                });

                // Linear search around pos to find the matching edge
                let mut found_at = None;

                // Search forward from pos
                for i in pos..status.len() {
                    let entry = &status[i];

                    if entry.edge == event.edge && entry.segment == event.segment {
                        found_at = Some(i);
                        break;
                    }
                }

                // If not found forward, search backward from pos
                if found_at.is_none() {
                    for i in (0..pos).rev() {
                        let entry = &status[i];

                        if entry.edge == event.edge && entry.segment == event.segment {
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
            SweepLineEventType::Start => {
                // Find insertion point in status
                let event_point = geometry.sweep_line_event_point(event);

                let pos = status.partition_point(|status_entry| {
                    let ord = geometry.sweep_line_segment_cmp(
                        status_entry.edge,
                        status_entry.segment,
                        status_entry.chain,
                        event_point,
                    );
                    // We are looking for this point
                    // Less -|- Equal -|- Greater
                    //                    ^
                    matches!(ord, Ordering::Less | Ordering::Equal)
                });

                // Calculate winding number
                let winding_below = if pos > 0 {
                    status[pos - 1].winding_above
                } else {
                    0
                };

                let winding_above = winding_below
                    + match event.chain {
                        SweepLineChain::Bottom => 1,
                        SweepLineChain::Top => -1,
                    };

                // Insert into status
                status.insert(
                    pos,
                    StatusEntry {
                        edge: event.edge,
                        segment: event.segment,
                        chain: event.chain,
                        winding_above,
                    },
                );

                // Determine if we should emit this edge
                // and in which direction
                if let Entry::Vacant(e) = output.entry(event.edge) {
                    let below_inside = winding_rule(winding_below);
                    let above_inside = winding_rule(winding_above);

                    if below_inside != above_inside {
                        if above_inside {
                            // Bottom edge: emit in forward direction
                            e.insert(Some(Direction::Forward));
                        } else {
                            // Top edge: emit in reverse direction
                            e.insert(Some(Direction::Reverse));
                        }
                    } else {
                        e.insert(None);
                    }
                }
            }
        }
    }

    // Output relevant edges, reversing them if necessary
    output
        .iter()
        .filter_map(|(&edge, dir)| match dir {
            Some(Direction::Forward) => Some(edge),
            Some(Direction::Reverse) => Some(edge.reversed()),
            None => None,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clean;
    use crate::kernel::polyline::F32 as Kernel;

    #[test]
    fn test_simple_square() {
        // Simple square
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        // Positive winding rule
        let result = clip(&mut kernel, edges.iter().copied(), |w| w > 0);

        // All edges should be in the output
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_empty_input() {
        // No edges
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [1.0, 0.0]]);
        let edges: Vec<(u32, u32)> = vec![];

        let result = clip(&mut kernel, edges.iter().copied(), |w| w > 0);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_overlapping_squares() {
        // Two overlapping squares that should be merged
        let mut kernel = Kernel::new(vec![
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
        ]);

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
        let cleaned_edges = clean(&mut kernel, edges.iter().copied());
        assert!(cleaned_edges.len() == 12);

        let union = clip(&mut kernel, cleaned_edges.iter().copied(), |w| w > 0);
        // 4 edges lie on the interior and should have been clipped
        assert!(union.len() == 8);

        let intersection = clip(&mut kernel, cleaned_edges.iter().copied(), |w| w > 1);
        // Intersection is a square
        assert!(intersection.len() == 4);
    }

    #[test]
    fn test_square_with_hole() {
        // Outer square and inner square (hole)
        let mut kernel = Kernel::new(vec![
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
        ]);

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

        // Positive winding rule shouldn't change kernel
        let result = clip(&mut kernel, edges.iter().copied(), |w| w > 0);
        assert_eq!(result.len(), 8);

        // Non-positive winding rule also shouldn't change kernel
        let result = clip(&mut kernel, edges.iter().copied(), |w| w <= 0);
        assert_eq!(result.len(), 8);
    }

    #[test]
    fn test_empty_output() {
        // Square with winding number that will be filtered out
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        // Rule that filters out everything
        let result = clip(&mut kernel, edges.iter().copied(), |w| w > 1);
        assert_eq!(result.len(), 0);
    }
}

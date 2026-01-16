use std::collections::{BTreeMap, btree_map::Entry};

use crate::kernel::{Edge, Kernel};
use crate::sweep_line::{
    SweepLineChain, SweepLineEvent, SweepLineEventType, SweepLineStatus, SweepLineStatusEntry,
};

/// Error type for clipping operations
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClippingError {
    /// Invalid topology encountered during clipping
    /// For example:
    /// * segment not found in sweep line status when expected
    /// * inconsistent winding numbers
    InvalidTopology,
}

impl std::fmt::Display for ClippingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ClippingError::InvalidTopology => {
                write!(f, "Invalid topology encountered during clipping")
            }
        }
    }
}

impl std::error::Error for ClippingError {}

/// Edge direction
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    Forward,
    Reverse,
}

/// Algorithm-specific data for status entries in clipping
#[derive(Debug, Clone)]
struct StatusData {
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
/// A list of edges wound positively around the "inside" area, or a ClippingError if
/// invalid topology is encountered.
pub fn clip<K: Kernel>(
    geometry: &mut K,
    edges: impl Iterator<Item = K::Edge>,
    winding_rule: impl Fn(i32) -> bool,
) -> Result<Vec<K::Edge>, ClippingError> {
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
    events.sort_by(|a, b| geometry.sweep_line_event_cmp(a, b));
    events
}

/// Run the sweep line algorithm
fn sweep_line<K: Kernel>(
    geometry: &mut K,
    events: Vec<SweepLineEvent<K>>,
    winding_rule: impl Fn(i32) -> bool,
) -> Result<Vec<K::Edge>, ClippingError> {
    let mut status = SweepLineStatus::<K, StatusData>::new();
    let mut output: BTreeMap<K::Edge, Direction> = BTreeMap::new();

    for event in events.iter() {
        match event.event_type {
            SweepLineEventType::End => {
                // Find and remove the edge from status
                let event_point = geometry.sweep_line_event_point(event);
                status
                    .remove(geometry, event_point, &event.segment)
                    .ok_or(ClippingError::InvalidTopology)?;
            }
            SweepLineEventType::Start => {
                // Find insertion point in status
                let event_point = geometry.sweep_line_event_point(event);

                // Calculate winding number
                let winding_below = status
                    .get_below(geometry, event_point)
                    .map(|entry| entry.data.winding_above)
                    .unwrap_or(0);

                let winding_above = winding_below
                    + match event.segment.chain {
                        SweepLineChain::Bottom => 1,
                        SweepLineChain::Top => -1,
                    };

                // Insert into status
                status.insert(
                    geometry,
                    event_point,
                    SweepLineStatusEntry::new(event.segment, StatusData { winding_above }),
                );

                // Determine if we should emit this edge
                // and in which direction
                let below_inside = winding_rule(winding_below);
                let above_inside = winding_rule(winding_above);

                if below_inside != above_inside {
                    let dir = if above_inside {
                        // Bottom edge
                        match event.segment.chain {
                            SweepLineChain::Bottom => Direction::Forward,
                            SweepLineChain::Top => Direction::Reverse,
                        }
                    } else {
                        // Top edge
                        match event.segment.chain {
                            SweepLineChain::Top => Direction::Forward,
                            SweepLineChain::Bottom => Direction::Reverse,
                        }
                    };
                    // TODO: avoid this panic possibility
                    // by handling each edge only once,
                    // and using its multiplicity in the winding rule check.
                    match output.entry(event.segment.edge) {
                        Entry::Vacant(e) => {
                            e.insert(dir);
                        }
                        Entry::Occupied(e) => match (e.get(), dir) {
                            (Direction::Forward, Direction::Reverse)
                            | (Direction::Reverse, Direction::Forward) => {
                                e.remove();
                            }
                            _ => {
                                panic!("Invalid topology encountered during clipping");
                            }
                        },
                    }
                }
            }
        }
    }

    // Output relevant edges, reversing them if necessary
    Ok(output
        .iter()
        .map(|(&edge, dir)| match dir {
            Direction::Forward => edge,
            Direction::Reverse => edge.reversed(),
        })
        .collect())
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
        let result = clip(&mut kernel, edges.iter().copied(), |w| w > 0).unwrap();

        // All edges should be in the output
        assert_eq!(result.len(), 4);
    }

    #[test]
    fn test_empty_input() {
        // No edges
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [1.0, 0.0]]);
        let edges: Vec<(u32, u32)> = vec![];

        let result = clip(&mut kernel, edges.iter().copied(), |w| w > 0).unwrap();
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

        let union = clip(&mut kernel, cleaned_edges.iter().copied(), |w| w > 0).unwrap();
        // 4 edges lie on the interior and should have been clipped
        assert!(union.len() == 8);

        let intersection = clip(&mut kernel, cleaned_edges.iter().copied(), |w| w > 1).unwrap();
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
        let result = clip(&mut kernel, edges.iter().copied(), |w| w > 0).unwrap();
        assert_eq!(result.len(), 8);

        // Non-positive winding rule also shouldn't change kernel
        let result = clip(&mut kernel, edges.iter().copied(), |w| w <= 0).unwrap();
        assert_eq!(result.len(), 8);
    }

    #[test]
    fn test_empty_output() {
        // Square with winding number that will be filtered out
        let mut kernel = Kernel::new(vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        // Rule that filters out everything
        let result = clip(&mut kernel, edges.iter().copied(), |w| w > 1).unwrap();
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_multiplicity() {
        let mut kernel = Kernel::new(vec![
            [-0.28425846, -0.37691897],
            [2.0, 0.0],
            [2.0, 2.0],
            [-0.020360839, 0.06198269],
            [2.939571, 0.18600994],
            [2.4747195, 1.9946594],
            [2.0, 1.626944],
            [2.0, 0.14663999],
            [0.5973753, 0.540478],
            [0.14415829, 0.068876386],
        ]);
        let edges = vec![
            (0, 1),
            (1, 7),
            (2, 8),
            (4, 5),
            (5, 6),
            (6, 2),
            (6, 8),
            (7, 4),
            (7, 6),
            (8, 9),
            (8, 9),
            (9, 0),
            (9, 7),
        ];
        let result = clip(&mut kernel, edges.iter().copied(), |w| w > 0).unwrap();

        // Result should have one copy of edge 8->9
        assert_eq!(result.iter().filter(|&&e| e == (8, 9)).count(), 1);
    }
}

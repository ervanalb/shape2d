use std::cmp::Ordering;

use crate::kernel::Kernel;

/// Chain designation for a segment in the sweep line algorithm
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SweepLineChain {
    Bottom, // Bottom edges go from left-to-right
    Top,    // Top edges go from right-to-left
}

/// Event type in the sweep line algorithm
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SweepLineEventType {
    End, // End comes before Start in sorting
    Start,
}

impl SweepLineEventType {
    pub fn other(self) -> Self {
        match self {
            SweepLineEventType::End => SweepLineEventType::Start,
            SweepLineEventType::Start => SweepLineEventType::End,
        }
    }
}

/// A segment in the sweep line algorithm
#[derive(Debug)]
pub struct SweepLineSegment<G: Kernel> {
    /// The edge containing this segment
    pub edge: G::Edge,
    /// The specific segment of the edge
    pub portion: G::SweepLineEdgePortion,
    /// Whether this segment is the top or bottom of an enclosed area
    pub chain: SweepLineChain,
}

impl<G: Kernel> Clone for SweepLineSegment<G> {
    fn clone(&self) -> Self {
        Self {
            edge: self.edge,
            portion: self.portion,
            chain: self.chain,
        }
    }
}

impl<G: Kernel> Copy for SweepLineSegment<G> {}

impl<G: Kernel> SweepLineSegment<G> {
    pub fn new(edge: G::Edge, portion: G::SweepLineEdgePortion, chain: SweepLineChain) -> Self {
        Self {
            edge,
            portion,
            chain,
        }
    }
}

/// An event in the sweep line algorithm
#[derive(Debug, Clone)]
pub struct SweepLineEvent<G: Kernel> {
    /// The type of event (start or end)
    pub event_type: SweepLineEventType,
    /// The segment this event is for
    pub segment: SweepLineSegment<G>,
}

/// An entry in the sweep-line status structure
///
/// Bundles a SweepLineSegment with algorithm-specific data
#[derive(Debug, Clone)]
pub struct SweepLineStatusEntry<G: Kernel, T> {
    /// The sweep-line segment (edge, segment, chain)
    pub segment: SweepLineSegment<G>,
    /// Algorithm-specific data
    pub data: T,
}

impl<G: Kernel, T> SweepLineStatusEntry<G, T> {
    pub fn new(segment: SweepLineSegment<G>, data: T) -> Self {
        Self { segment, data }
    }
}

/// A sweep-line status structure that maintains a sorted list of active segments
///
/// This wraps a Vec and provides methods for finding insertion points,
/// finding specific entries, and finding entries below a point.
pub struct SweepLineStatus<G: Kernel, T> {
    entries: Vec<SweepLineStatusEntry<G, T>>,
}

impl<G: Kernel, T> SweepLineStatus<G, T> {
    /// Create a new empty status structure
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Find the insertion point for a new segment at the given event point
    ///
    /// Returns the index where a new entry should be inserted to maintain sorted order.
    //TODO(Claude): return a StatusInsertionPoint struct which wraps the returned usize here.
    //This struct should borrow the Status mutably, and provide methods below() (which returns an
    //immutable reference to the status entry below the insertion point) and insert() (which consumes the
    //InsertionPoint and runs insert() on the vec.)
    pub fn find_insertion_point(&self, geometry: &G, event_point: G::SweepLineEventPoint) -> usize {
        self.entries.partition_point(|entry| {
            let ord = geometry.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less | Ordering::Equal)
        })
    }

    /// Find a specific status entry by matching segment
    ///
    /// Uses binary search followed by linear search to handle degenerate cases.
    /// Panics if the entry is not found.
    pub fn find_entry(
        &self,
        geometry: &G,
        event_point: G::SweepLineEventPoint,
        target_segment: &SweepLineSegment<G>,
    ) -> usize {
        let pos = self.entries.partition_point(|entry| {
            let ord = geometry.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less)
        });

        // Linear search around pos to find matching edge and portion
        for i in pos..self.entries.len() {
            if self.entries[i].segment.edge == target_segment.edge
                && self.entries[i].segment.portion == target_segment.portion
            {
                return i;
            }
        }
        for i in (0..pos).rev() {
            if self.entries[i].segment.edge == target_segment.edge
                && self.entries[i].segment.portion == target_segment.portion
            {
                return i;
            }
        }

        panic!("Status entry not found");
    }

    /// Find the status entry directly below a vertex
    ///
    /// Returns the index of the entry immediately below the event point.
    /// Panics if there is no segment below the vertex.
    //TODO(Claude): This method should no longer be necessary. Instead, callers should use
    //find_insertion_point().below().
    pub fn find_below(&self, geometry: &G, event_point: G::SweepLineEventPoint) -> usize {
        let pos = self.entries.partition_point(|entry| {
            let ord = geometry.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less)
        });

        if pos == 0 {
            panic!("No segment below vertex");
        }
        pos - 1
    }

    /// Insert an entry at the specified position
    //TODO(Claude): Remove this method since we will use find_insertion_point().insert()
    pub fn insert(&mut self, pos: usize, entry: SweepLineStatusEntry<G, T>) {
        self.entries.insert(pos, entry);
    }

    /// Remove and return the entry at the specified position
    //TODO(Claude): Combine this with find_entry() since we always want to remove an entry after
    //finding it.
    pub fn remove(&mut self, pos: usize) -> SweepLineStatusEntry<G, T> {
        self.entries.remove(pos)
    }
}

impl<G: Kernel, T> Default for SweepLineStatus<G, T> {
    fn default() -> Self {
        Self::new()
    }
}

// Implement indexing for convenient access to the data
impl<G: Kernel, T> std::ops::Index<usize> for SweepLineStatus<G, T> {
    type Output = SweepLineStatusEntry<G, T>;

    fn index(&self, index: usize) -> &SweepLineStatusEntry<G, T> {
        &self.entries[index]
    }
}

impl<G: Kernel, T> std::ops::IndexMut<usize> for SweepLineStatus<G, T> {
    fn index_mut(&mut self, index: usize) -> &mut SweepLineStatusEntry<G, T> {
        &mut self.entries[index]
    }
}

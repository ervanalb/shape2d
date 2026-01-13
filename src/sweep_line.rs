use std::cmp::Ordering;

use crate::kernel::Kernel;

/// Chain designation for a segment in the sweep line algorithm
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
pub struct SweepLineSegment<G: Kernel> {
    /// The edge containing this segment
    pub edge: G::Edge,
    /// The specific segment of the edge
    pub portion: G::SweepLineEdgePortion,
    /// Whether this segment is the top or bottom of an enclosed area
    pub chain: SweepLineChain,
}

impl<G: Kernel> std::fmt::Debug for SweepLineSegment<G> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SweepLineSegment")
            .field("edge", &self.edge)
            .field("portion", &self.portion)
            .field("chain", &self.chain)
            .finish()
    }
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

impl<G: Kernel> PartialEq for SweepLineSegment<G> {
    fn eq(&self, other: &Self) -> bool {
        self.edge == other.edge && self.portion == other.portion && self.chain == other.chain
    }
}

impl<G: Kernel> Eq for SweepLineSegment<G> {}

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
#[derive(Clone)]
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

impl<G: Kernel, T: std::fmt::Debug> std::fmt::Debug for SweepLineStatusEntry<G, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SweepLineStatus")
            .field("segment", &self.segment)
            .field("data", &self.data)
            .finish()
    }
}

/// A sweep-line status structure that maintains a sorted list of active segments
///
/// This wraps a Vec and provides methods for finding insertion points,
/// finding specific entries, and finding entries below a point.
#[derive(Clone)]
pub struct SweepLineStatus<G: Kernel, T> {
    entries: Vec<SweepLineStatusEntry<G, T>>,
}

impl<G: Kernel, T: std::fmt::Debug> std::fmt::Debug for SweepLineStatus<G, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SweepLineStatus")
            .field("entries", &self.entries)
            .finish()
    }
}

impl<G: Kernel, T> SweepLineStatus<G, T> {
    /// Create a new empty status structure
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    fn search(
        &self,
        geometry: &G,
        event_point: G::SweepLineEventPoint,
        target_segment: &SweepLineSegment<G>,
    ) -> Option<usize> {
        let start_i = self.entries.partition_point(|entry| {
            let ord = geometry.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less)
        });
        let end_i = self.entries.partition_point(|entry| {
            let ord = geometry.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less | Ordering::Equal)
        });

        // Linear search over all equal segments
        for i in start_i..end_i {
            if &self.entries[i].segment == target_segment {
                return Some(i);
            }
        }
        None
    }

    /// Remove a segment from the status and return it
    ///
    /// Panics if the segment is not found.
    pub fn remove(
        &mut self,
        geometry: &G,
        event_point: G::SweepLineEventPoint,
        target_segment: &SweepLineSegment<G>,
    ) -> SweepLineStatusEntry<G, T> {
        let i = self
            .search(geometry, event_point, target_segment)
            .expect("Segment not found in status");
        return self.entries.remove(i);
    }

    /// Get a reference to the segment equal or below the given event point
    /// Returns None if there is no segment below.
    pub fn get_equal_or_below(
        &self,
        geometry: &G,
        event_point: G::SweepLineEventPoint,
    ) -> Option<&SweepLineStatusEntry<G, T>> {
        let pos = self.entries.partition_point(|entry| {
            let ord = geometry.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less | Ordering::Equal)
        });

        if pos == 0 {
            None
        } else {
            Some(&self.entries[pos - 1])
        }
    }

    /// Get a mutable reference to the segment strictly below the given event point
    ///
    /// Uses strict ordering (Ordering::Less only) to find the segment below.
    /// Returns None if there is no segment below.
    pub fn get_below_mut(
        &mut self,
        geometry: &G,
        event_point: G::SweepLineEventPoint,
    ) -> Option<&mut SweepLineStatusEntry<G, T>> {
        let pos = self.entries.partition_point(|entry| {
            let ord = geometry.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less)
        });

        if pos == 0 {
            None
        } else {
            Some(&mut self.entries[pos - 1])
        }
    }

    /// Insert a new entry at the appropriate position for the given event point
    /// behind any entries that are equal
    pub fn insert_back(
        &mut self,
        geometry: &G,
        event_point: G::SweepLineEventPoint,
        entry: SweepLineStatusEntry<G, T>,
    ) {
        let pos = self.entries.partition_point(|e| {
            let ord = geometry.sweep_line_segment_cmp(&e.segment, event_point);
            matches!(ord, Ordering::Less | Ordering::Equal)
        });
        self.entries.insert(pos, entry);
    }

    /// Insert a new entry at the appropriate position for the given event point
    /// in front of any entries that are equal
    pub fn insert_front(
        &mut self,
        geometry: &G,
        event_point: G::SweepLineEventPoint,
        entry: SweepLineStatusEntry<G, T>,
    ) {
        let pos = self.entries.partition_point(|e| {
            let ord = geometry.sweep_line_segment_cmp(&e.segment, event_point);
            matches!(ord, Ordering::Less)
        });
        self.entries.insert(pos, entry);
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

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
pub struct SweepLineSegment<K: Kernel> {
    /// The edge containing this segment
    pub edge: K::Edge,
    /// The specific segment of the edge
    pub portion: K::SweepLineEdgePortion,
    /// Whether this segment is the top or bottom of an enclosed area
    pub chain: SweepLineChain,
}

impl<K: Kernel> std::fmt::Debug for SweepLineSegment<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SweepLineSegment")
            .field("edge", &self.edge)
            .field("portion", &self.portion)
            .field("chain", &self.chain)
            .finish()
    }
}

impl<K: Kernel> Clone for SweepLineSegment<K> {
    fn clone(&self) -> Self {
        Self {
            edge: self.edge,
            portion: self.portion,
            chain: self.chain,
        }
    }
}
impl<K: Kernel> Copy for SweepLineSegment<K> {}

impl<K: Kernel> PartialEq for SweepLineSegment<K> {
    fn eq(&self, other: &Self) -> bool {
        self.edge == other.edge && self.portion == other.portion && self.chain == other.chain
    }
}

impl<K: Kernel> Eq for SweepLineSegment<K> {}

impl<K: Kernel> SweepLineSegment<K> {
    pub fn new(edge: K::Edge, portion: K::SweepLineEdgePortion, chain: SweepLineChain) -> Self {
        Self {
            edge,
            portion,
            chain,
        }
    }
}

/// An event in the sweep line algorithm
#[derive(Debug, Clone)]
pub struct SweepLineEvent<K: Kernel> {
    /// The type of event (start or end)
    pub event_type: SweepLineEventType,
    /// The segment this event is for
    pub segment: SweepLineSegment<K>,
}

/// An entry in the sweep-line status structure
///
/// Bundles a SweepLineSegment with algorithm-specific data
#[derive(Clone)]
pub struct SweepLineStatusEntry<K: Kernel, T> {
    /// The sweep-line segment (edge, segment, chain)
    pub segment: SweepLineSegment<K>,
    /// Algorithm-specific data
    pub data: T,
}

impl<K: Kernel, T> SweepLineStatusEntry<K, T> {
    pub fn new(segment: SweepLineSegment<K>, data: T) -> Self {
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
pub struct SweepLineStatus<K: Kernel, T> {
    entries: Vec<SweepLineStatusEntry<K, T>>,
}

impl<G: Kernel, T: std::fmt::Debug> std::fmt::Debug for SweepLineStatus<G, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SweepLineStatus")
            .field("entries", &self.entries)
            .finish()
    }
}

impl<G: Kernel, T: std::fmt::Debug> SweepLineStatus<G, T> {
    /// Create a new empty status structure
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Remove a segment from the status and return it
    ///
    /// Returns None if the segment is not found or doesn't match the target.
    pub fn remove(
        &mut self,
        kernel: &G,
        event_point: G::SweepLineEventPoint,
        target_segment: &SweepLineSegment<G>,
    ) -> Option<SweepLineStatusEntry<G, T>> {
        let i = self.entries.partition_point(|entry| {
            let ord = kernel.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less)
        });

        // Check if the index is valid and the segment matches
        if i < self.entries.len() && &self.entries[i].segment == target_segment {
            Some(self.entries.remove(i))
        } else {
            None
        }
    }

    /// Get a reference to the segment equal or below the given event point
    /// Returns None if there is no segment below.
    pub fn get_below(
        &self,
        kernel: &G,
        event_point: G::SweepLineEventPoint,
    ) -> Option<&SweepLineStatusEntry<G, T>> {
        let pos = self.entries.partition_point(|entry| {
            let ord = kernel.sweep_line_segment_cmp(&entry.segment, event_point);
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
        kernel: &G,
        event_point: G::SweepLineEventPoint,
    ) -> Option<&mut SweepLineStatusEntry<G, T>> {
        let pos = self.entries.partition_point(|entry| {
            let ord = kernel.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less | Ordering::Equal)
        });

        if pos == 0 {
            None
        } else {
            Some(&mut self.entries[pos - 1])
        }
    }

    /// Remove a segment from the status and return it,
    /// along a reference to the segment below it
    ///
    /// Returns None if the segment is not found or doesn't match the target.
    pub fn get_below_mut_and_remove(
        &mut self,
        kernel: &G,
        event_point: G::SweepLineEventPoint,
        target_segment: &SweepLineSegment<G>,
    ) -> Option<(
        Option<&mut SweepLineStatusEntry<G, T>>,
        SweepLineStatusEntry<G, T>,
    )> {
        let i = self.entries.partition_point(|entry| {
            let ord = kernel.sweep_line_segment_cmp(&entry.segment, event_point);
            matches!(ord, Ordering::Less)
        });

        // Check if the index is valid and the segment matches
        if i < self.entries.len() && &self.entries[i].segment == target_segment {
            let removed = self.entries.remove(i);
            let below = if i > 0 {
                Some(&mut self.entries[i - 1])
            } else {
                None
            };
            Some((below, removed))
        } else {
            None
        }
    }

    /// Insert a new entry at the appropriate position for the given event point
    /// behind any entries that are equal.
    pub fn insert(
        &mut self,
        kernel: &G,
        event_point: G::SweepLineEventPoint,
        entry: SweepLineStatusEntry<G, T>,
    ) {
        let pos = self.entries.partition_point(|e| {
            let ord = kernel.sweep_line_segment_cmp(&e.segment, event_point);
            matches!(ord, Ordering::Less | Ordering::Equal)
        });
        self.entries.insert(pos, entry);
    }
}

impl<G: Kernel, T: std::fmt::Debug> Default for SweepLineStatus<G, T> {
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

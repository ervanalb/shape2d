/// Axis-aligned bounding rectangle with u16 coordinates
/// Inclusive on all bounds
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rect {
    pub min: [u16; 2],
    pub max: [u16; 2],
}

impl Default for Rect {
    /// This default value is an invalid rect
    /// that is designed so that when combined with any valid rect,
    /// it will yield that rect back.
    fn default() -> Self {
        Self {
            min: [u16::MAX, u16::MAX],
            max: [0, 0],
        }
    }
}

impl Rect {
    /// Combine this rectangle with another, returning the smallest rectangle
    /// that bounds both
    pub fn combine(&self, other: &Self) -> Self {
        Self {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
            ],
        }
    }

    /// Check if two rectangles overlap (including touching on any side or corner)
    /// Returns false if either rect is invalid (max < min)
    pub fn overlaps(&self, other: &Self) -> bool {
        self.min[0] <= self.max[0]
            && self.min[1] <= self.max[1]
            && other.min[0] <= other.max[0]
            && other.min[1] <= other.max[1]
            && self.min[0] <= other.max[0]
            && self.max[0] >= other.min[0]
            && self.min[1] <= other.max[1]
            && self.max[1] >= other.min[1]
    }

    /// Compute the Hilbert curve value for the center point of this rectangle
    /// u16his is slightly lossy as it maps (u32, u32) -> u32
    pub fn hilbert_value(&self) -> u32 {
        // Calculate center point, being careful about overflow
        let center_x = ((self.min[0] as u32 + self.max[0] as u32) / 2) as u16;
        let center_y = ((self.min[1] as u32 + self.max[1] as u32) / 2) as u16;
        crate::hilbert::xy_to_hilbert(center_x, center_y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_combine() {
        let r1 = Rect { min: [0, 0], max: [10, 10] };
        let r2 = Rect { min: [5, 5], max: [15, 15] };
        let combined = r1.combine(&r2);

        assert_eq!(combined.min, [0, 0]);
        assert_eq!(combined.max, [15, 15]);
    }

    #[test]
    fn test_combine_with_default() {
        let r1 = Rect::default();
        let r2 = Rect { min: [5, 5], max: [15, 15] };
        let combined = r1.combine(&r2);

        assert_eq!(combined.min, [5, 5]);
        assert_eq!(combined.max, [15, 15]);
    }

    #[test]
    fn test_is_overlapping() {
        let r1 = Rect { min: [0, 0], max: [10, 10] };
        let r2 = Rect { min: [5, 5], max: [15, 15] };
        let r3 = Rect { min: [20, 20], max: [30, 30] };

        assert!(r1.overlaps(&r2));
        assert!(r2.overlaps(&r1));
        assert!(!r1.overlaps(&r3));
        assert!(!r3.overlaps(&r1));
    }

    #[test]
    fn test_overlapping_touching() {
        let r1 = Rect { min: [0, 0], max: [10, 10] };
        let r2 = Rect { min: [10, 0], max: [20, 10] }; // Touching on right edge

        assert!(r1.overlaps(&r2));
        assert!(r2.overlaps(&r1));
    }

    #[test]
    fn test_overlapping_corner() {
        let r1 = Rect { min: [0, 0], max: [10, 10] };
        let r2 = Rect { min: [10, 10], max: [20, 20] }; // Touching at corner

        assert!(r1.overlaps(&r2));
        assert!(r2.overlaps(&r1));
    }
}

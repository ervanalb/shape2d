/// Axis-aligned bounding rectangle with u16 coordinates
/// Inclusive on all bounds
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rect {
    pub min_x: u16,
    pub max_x: u16,
    pub min_y: u16,
    pub max_y: u16,
}

impl Default for Rect {
    fn default() -> Self {
        Self {
            min_x: 0,
            max_x: u16::MAX,
            min_y: 0,
            max_y: u16::MAX,
        }
    }
}

impl Rect {
    /// Create a new rectangle
    pub fn new(min_x: u16, max_x: u16, min_y: u16, max_y: u16) -> Self {
        Self {
            min_x,
            max_x,
            min_y,
            max_y,
        }
    }

    /// Combine this rectangle with another, returning the smallest rectangle
    /// that bounds both
    pub fn combine(&self, other: &Self) -> Self {
        Self {
            min_x: if self.min_x < other.min_x {
                self.min_x
            } else {
                other.min_x
            },
            max_x: if self.max_x > other.max_x {
                self.max_x
            } else {
                other.max_x
            },
            min_y: if self.min_y < other.min_y {
                self.min_y
            } else {
                other.min_y
            },
            max_y: if self.max_y > other.max_y {
                self.max_y
            } else {
                other.max_y
            },
        }
    }

    /// Check if two rectangles overlap (including touching on any side or corner)
    pub fn overlaps(&self, other: &Self) -> bool {
        self.min_x <= other.max_x
            && self.max_x >= other.min_x
            && self.min_y <= other.max_y
            && self.max_y >= other.min_y
    }

    /// Compute the Hilbert curve value for the center point of this rectangle
    /// u16his is slightly lossy as it maps (u32, u32) -> u32
    pub fn hilbert_value(&self) -> u32 {
        // Calculate center point, being careful about overflow
        let center_x = ((self.min_x as u32 + self.max_x as u32) / 2) as u16;
        let center_y = ((self.min_y as u32 + self.max_y as u32) / 2) as u16;
        crate::hilbert::xy_to_hilbert(center_x, center_y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_combine() {
        let r1 = Rect::new(0, 10, 0, 10);
        let r2 = Rect::new(5, 15, 5, 15);
        let combined = r1.combine(&r2);

        assert_eq!(combined.min_x, 0);
        assert_eq!(combined.max_x, 15);
        assert_eq!(combined.min_y, 0);
        assert_eq!(combined.max_y, 15);
    }

    #[test]
    fn test_is_overlapping() {
        let r1 = Rect::new(0, 10, 0, 10);
        let r2 = Rect::new(5, 15, 5, 15);
        let r3 = Rect::new(20, 30, 20, 30);

        assert!(r1.overlaps(&r2));
        assert!(r2.overlaps(&r1));
        assert!(!r1.overlaps(&r3));
        assert!(!r3.overlaps(&r1));
    }

    #[test]
    fn test_overlapping_touching() {
        let r1 = Rect::new(0, 10, 0, 10);
        let r2 = Rect::new(10, 20, 0, 10); // u16ouching on right edge

        assert!(r1.overlaps(&r2));
        assert!(r2.overlaps(&r1));
    }

    #[test]
    fn test_overlapping_corner() {
        let r1 = Rect::new(0, 10, 0, 10);
        let r2 = Rect::new(10, 20, 10, 20); // u16ouching at corner

        assert!(r1.overlaps(&r2));
        assert!(r2.overlaps(&r1));
    }
}

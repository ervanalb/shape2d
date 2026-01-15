use std::collections::BTreeMap;

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
            min: [self.min[0].min(other.min[0]), self.min[1].min(other.min[1])],
            max: [self.max[0].max(other.max[0]), self.max[1].max(other.max[1])],
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
        xy_to_hilbert(center_x, center_y)
    }
}

/// A node in the packed R-tree array
#[derive(Debug, Clone, Default)]
pub struct RTreeNode {
    pub bbox: Rect,
    pub hilbert_value: u32,
}

/// Packed array Hilbert R-tree for spatial indexing
#[derive(Debug, Clone, Default)]
pub struct RTree {
    nodes: Vec<RTreeNode>,
    first_leaf_index: usize,
}

impl RTree {
    /// Build a new R-tree from a mapping of hilbert values to bounding rectangles
    pub fn build(hilbert_to_rect: BTreeMap<u32, Rect>) -> Self {
        let num_leaves = hilbert_to_rect.len();

        if num_leaves == 0 {
            return Default::default();
        }

        // Calculate first leaf index: L rounded up to next power of 2, minus 1
        let leaf_capacity = num_leaves.next_power_of_two();
        let first_leaf_index = leaf_capacity - 1;

        // Allocate space for all nodes, including a full leaf level
        // to avoid ever having to do bounds checks.
        // The "empty" nodes have hilbert value 0 and invalid bounding boxes.
        // These values should combine correctly during the tree building step.
        let mut nodes = vec![Default::default(); first_leaf_index + leaf_capacity];

        // Populate leaf nodes in hilbert-value order
        let mut leaf_idx = first_leaf_index;
        for (hilbert_value, bbox) in hilbert_to_rect {
            nodes[leaf_idx] = RTreeNode {
                bbox,
                hilbert_value,
            };
            leaf_idx += 1;
        }

        // Build parent nodes bottom-up
        if first_leaf_index > 0 {
            for i in (0..first_leaf_index).rev() {
                let left_child = &nodes[Self::left_child_idx(i)];
                let right_child = &nodes[Self::right_child_idx(i)];

                nodes[i] = RTreeNode {
                    bbox: left_child.bbox.combine(&right_child.bbox),
                    hilbert_value: left_child.hilbert_value.max(right_child.hilbert_value),
                };
            }
        }

        Self {
            nodes,
            first_leaf_index,
        }
    }

    fn left_child_idx(i: usize) -> usize {
        2 * i + 1
    }

    fn right_child_idx(i: usize) -> usize {
        2 * i + 2
    }

    fn parent_idx(i: usize) -> usize {
        (i - 1) / 2
    }

    /// Find the leaf node index for a given hilbert value
    /// Panics if the hilbert value does not exist in the R-tree
    fn find_leaf(&self, hilbert_value: u32) -> usize {
        if self.nodes.is_empty() {
            panic!("R-tree is empty");
        }

        let mut i = 0;

        // Walk down the tree
        while i < self.first_leaf_index {
            let left_child_idx = Self::left_child_idx(i);
            let right_child_idx = Self::right_child_idx(i);

            let left_hilbert = self.nodes[left_child_idx].hilbert_value;

            if hilbert_value <= left_hilbert {
                i = left_child_idx;
            } else {
                i = right_child_idx;
            }
        }

        // Check if we found the correct leaf
        if self.nodes[i].hilbert_value != hilbert_value {
            panic!("Hilbert value not found in R-tree");
        }
        i
    }

    /// Enlarge the bounding box for a given hilbert value
    /// Panics if the given hilbert value does not exist in the R-tree
    pub fn enlarge(&mut self, hilbert_value: u32, bbox: Rect) {
        let mut idx = self.find_leaf(hilbert_value);

        let mut combined_bbox = bbox;

        loop {
            // Combine with existing bbox
            let bbox = &mut self.nodes[idx].bbox;
            combined_bbox = bbox.combine(&combined_bbox);

            if *bbox == combined_bbox {
                return;
            }

            *bbox = combined_bbox;

            // Stop if we are at the root
            if idx == 0 {
                break;
            }

            // Propagate changes up the tree
            idx = Self::parent_idx(idx);
        }
    }

    /// Search for all leaves whose rectangles overlap with the test rectangle
    pub fn search(&self, test_rect: &Rect) -> RTreeSearchIterator<'_> {
        RTreeSearchIterator::new(self, *test_rect)
    }
}

/// Iterator for R-tree spatial search
pub struct RTreeSearchIterator<'a> {
    tree: &'a RTree,
    test_rect: Rect,
    cur: Option<usize>,
}

impl<'a> RTreeSearchIterator<'a> {
    // Given an overlapping node,
    // returns the next overlapping node in preorder traversal of the tree
    fn advance(tree: &RTree, test_rect: &Rect, mut cur: usize) -> Option<usize> {
        if cur < tree.first_leaf_index {
            let left_child_idx = RTree::left_child_idx(cur);
            if tree.nodes[left_child_idx].bbox.overlaps(test_rect) {
                return Some(left_child_idx);
            }
            let right_child_idx = RTree::right_child_idx(cur);
            if tree.nodes[right_child_idx].bbox.overlaps(test_rect) {
                return Some(right_child_idx);
            }
        }

        // No valid paths downward--go up until we find an unvisited right subtree
        while cur > 0 {
            let parent = RTree::parent_idx(cur);
            let left_sibling = RTree::left_child_idx(parent);
            let right_sibling = RTree::right_child_idx(parent);
            if cur == left_sibling && tree.nodes[right_sibling].bbox.overlaps(test_rect) {
                return Some(right_sibling);
            }
            cur = parent;
        }
        None
    }

    // Given an overlapping node (leaf or otherwise),
    // returns the next overlapping leaf node in preorder traversal of the tree
    fn advance_to_leaf(tree: &RTree, test_rect: &Rect, mut cur: usize) -> Option<usize> {
        while let Some(next) = Self::advance(tree, test_rect, cur) {
            if next >= tree.first_leaf_index {
                return Some(next);
            }
            cur = next;
        }
        None
    }

    fn new(tree: &'a RTree, test_rect: Rect) -> Self {
        let cur = if let Some(root_node) = tree.nodes.get(0) {
            if root_node.bbox.overlaps(&test_rect) {
                if tree.nodes.len() == 1 {
                    // The root is a leaf
                    Some(0)
                } else {
                    Self::advance_to_leaf(tree, &test_rect, 0)
                }
            } else {
                None
            }
        } else {
            None
        };

        Self {
            tree,
            test_rect,
            cur,
        }
    }
}

impl<'a> Iterator for RTreeSearchIterator<'a> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        let cur = self.cur?;
        let result = self.tree.nodes[cur].hilbert_value;
        self.cur = Self::advance_to_leaf(self.tree, &self.test_rect, cur);
        Some(result)
    }
}

/// Convert (x, y) coordinates to a Hilbert curve index
/// Input coordinates are u16 values, output is u32
// TODO(Eric): See if this can be replaced with clever bit fiddling
pub fn xy_to_hilbert(x: u16, y: u16) -> u32 {
    let mut x = x as u32;
    let mut y = y as u32;
    let mut d = 0u32;
    let n = 65536u32; // 2^16 for u16 coordinates

    /// Rotate/flip a quadrant appropriately
    fn rotate(n: u32, x: &mut u32, y: &mut u32, rx: bool, ry: bool) {
        if !ry {
            if rx {
                *x = n.saturating_sub(1).saturating_sub(*x);
                *y = n.saturating_sub(1).saturating_sub(*y);
            }

            // Swap x and y
            std::mem::swap(x, y);
        }
    }

    // Process from most significant bit to least
    let mut s = n / 2;
    while s > 0 {
        let rx = (x & s) != 0;
        let ry = (y & s) != 0;

        d += s * s * ((3 * (rx as u32)) ^ (ry as u32));

        // Rotate coordinates
        rotate(n, &mut x, &mut y, rx, ry);

        s /= 2;
    }

    d
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hilbert_corners() {
        // Test the four corners of the space
        let h00 = xy_to_hilbert(0, 0);
        let h01 = xy_to_hilbert(0, u16::MAX);
        let h10 = xy_to_hilbert(u16::MAX, 0);
        let h11 = xy_to_hilbert(u16::MAX, u16::MAX);

        // All should be different
        assert_ne!(h00, h01);
        assert_ne!(h00, h10);
        assert_ne!(h00, h11);
        assert_ne!(h01, h10);
        assert_ne!(h01, h11);
        assert_ne!(h10, h11);
    }

    #[test]
    fn test_hilbert_locality() {
        // Points that are close in 2D space should have similar Hilbert values
        let h1 = xy_to_hilbert(1000, 1000);
        let h2 = xy_to_hilbert(1001, 1000);
        let h3 = xy_to_hilbert(1000, 1001);
        let h_far = xy_to_hilbert(50000, 50000);

        // Nearby points should be closer in Hilbert space than far points
        let diff_near = h1.abs_diff(h2).min(h1.abs_diff(h3));
        let diff_far = h1.abs_diff(h_far);

        assert!(diff_near < diff_far);
    }

    #[test]
    fn test_combine() {
        let r1 = Rect {
            min: [0, 0],
            max: [10, 10],
        };
        let r2 = Rect {
            min: [5, 5],
            max: [15, 15],
        };
        let combined = r1.combine(&r2);

        assert_eq!(combined.min, [0, 0]);
        assert_eq!(combined.max, [15, 15]);
    }

    #[test]
    fn test_combine_with_default() {
        let r1 = Rect::default();
        let r2 = Rect {
            min: [5, 5],
            max: [15, 15],
        };
        let combined = r1.combine(&r2);

        assert_eq!(combined.min, [5, 5]);
        assert_eq!(combined.max, [15, 15]);
    }

    #[test]
    fn test_is_overlapping() {
        let r1 = Rect {
            min: [0, 0],
            max: [10, 10],
        };
        let r2 = Rect {
            min: [5, 5],
            max: [15, 15],
        };
        let r3 = Rect {
            min: [20, 20],
            max: [30, 30],
        };

        assert!(r1.overlaps(&r2));
        assert!(r2.overlaps(&r1));
        assert!(!r1.overlaps(&r3));
        assert!(!r3.overlaps(&r1));
    }

    #[test]
    fn test_overlapping_touching() {
        let r1 = Rect {
            min: [0, 0],
            max: [10, 10],
        };
        let r2 = Rect {
            min: [10, 0],
            max: [20, 10],
        }; // Touching on right edge

        assert!(r1.overlaps(&r2));
        assert!(r2.overlaps(&r1));
    }

    #[test]
    fn test_overlapping_corner() {
        let r1 = Rect {
            min: [0, 0],
            max: [10, 10],
        };
        let r2 = Rect {
            min: [10, 10],
            max: [20, 20],
        }; // Touching at corner

        assert!(r1.overlaps(&r2));
        assert!(r2.overlaps(&r1));
    }

    #[test]
    fn test_rtree_build() {
        let mut map = BTreeMap::new();
        map.insert(
            10,
            Rect {
                min: [0, 0],
                max: [10, 10],
            },
        );
        map.insert(
            20,
            Rect {
                min: [10, 10],
                max: [20, 20],
            },
        );
        map.insert(
            30,
            Rect {
                min: [20, 20],
                max: [30, 30],
            },
        );

        let tree = RTree::build(map);

        assert!(tree.nodes.len() > 0);
        assert_eq!(tree.first_leaf_index, 3); // 3 leaves -> next power of 2 is 4, minus 1 = 3
    }

    #[test]
    fn test_rtree_search() {
        let mut map = BTreeMap::new();
        map.insert(
            10,
            Rect {
                min: [0, 0],
                max: [10, 10],
            },
        );
        map.insert(
            20,
            Rect {
                min: [10, 10],
                max: [20, 20],
            },
        );
        map.insert(
            30,
            Rect {
                min: [20, 20],
                max: [30, 30],
            },
        );

        let tree = RTree::build(map);

        // Search for rect that overlaps first two
        let test_rect = Rect {
            min: [5, 5],
            max: [15, 15],
        };
        let results: Vec<u32> = tree.search(&test_rect).collect();

        assert!(results.len() == 2);
        assert!(results.contains(&10));
        assert!(results.contains(&20));
    }

    #[test]
    fn test_rtree_enlarge() {
        let mut map = BTreeMap::new();
        map.insert(
            10,
            Rect {
                min: [0, 0],
                max: [10, 10],
            },
        );
        let test_rect_a = Rect {
            min: [5, 5],
            max: [8, 8],
        };
        let test_rect_b = Rect {
            min: [16, 16],
            max: [17, 17],
        };

        let mut tree = RTree::build(map);

        // Check that rect A is found
        let results: Vec<u32> = tree.search(&test_rect_a).collect();
        assert_eq!(results, vec![10]);

        // Check that rect B is NOT found
        let results: Vec<u32> = tree.search(&test_rect_b).collect();
        assert_eq!(results, vec![]);

        // Enlarge the bounding box
        tree.enlarge(
            10,
            Rect {
                min: [15, 15],
                max: [20, 20],
            },
        );

        // Check that rect A is found
        let results: Vec<u32> = tree.search(&test_rect_a).collect();
        assert_eq!(results, vec![10]);

        // Check that rect B is found
        let results: Vec<u32> = tree.search(&test_rect_b).collect();
        assert_eq!(results, vec![10]);
    }

    #[test]
    fn test_rtree_empty() {
        let map = BTreeMap::new();
        let tree = RTree::build(map);

        assert_eq!(tree.nodes.len(), 0);

        let test_rect = Rect {
            min: [0, 0],
            max: [10, 10],
        };
        let results: Vec<u32> = tree.search(&test_rect).collect();
        assert_eq!(results.len(), 0);
    }
}

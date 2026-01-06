use crate::Rect;
use std::collections::BTreeMap;

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

        while idx > 0 {
            // Combine with existing bbox
            let bbox = &mut self.nodes[idx].bbox;
            combined_bbox = bbox.combine(&combined_bbox);

            if *bbox == combined_bbox {
                return;
            }

            *bbox = combined_bbox;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rtree_build() {
        let mut map = BTreeMap::new();
        map.insert(10, Rect::new(0, 10, 0, 10));
        map.insert(20, Rect::new(10, 20, 10, 20));
        map.insert(30, Rect::new(20, 30, 20, 30));

        let tree = RTree::build(map);

        assert!(tree.nodes.len() > 0);
        assert_eq!(tree.first_leaf_index, 3); // 3 leaves -> next power of 2 is 4, minus 1 = 3
    }

    #[test]
    fn test_rtree_search() {
        let mut map = BTreeMap::new();
        map.insert(10, Rect::new(0, 10, 0, 10));
        map.insert(20, Rect::new(10, 20, 10, 20));
        map.insert(30, Rect::new(20, 30, 20, 30));

        let tree = RTree::build(map);

        // Search for rect that overlaps first two
        let test_rect = Rect::new(5, 15, 5, 15);
        let results: Vec<u32> = tree.search(&test_rect).collect();

        assert!(results.contains(&10));
        assert!(results.contains(&20));
        assert!(!results.contains(&30));
    }

    #[test]
    fn test_rtree_enlarge() {
        let mut map = BTreeMap::new();
        map.insert(10, Rect::new(0, 10, 0, 10));

        let mut tree = RTree::build(map);

        let changed = tree.enlarge(10, Rect::new(0, 20, 0, 20));
        //assert!(changed);
        // TODO: since changed has been removed,
        // TODO: instead, we should run a search on the tree
        // TODO: to ensure that the given hilbert value was enlarged.

        // Enlarge with same rect should not change anything
        let changed = tree.enlarge(10, Rect::new(0, 15, 0, 15));
        //assert!(!changed);
        // TODO: same as above
    }

    #[test]
    fn test_rtree_empty() {
        let map = BTreeMap::new();
        let tree = RTree::build(map);

        assert_eq!(tree.nodes.len(), 0);

        let test_rect = Rect::new(0, 10, 0, 10);
        let results: Vec<u32> = tree.search(&test_rect).collect();
        assert_eq!(results.len(), 0);
    }
}

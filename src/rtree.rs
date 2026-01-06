use crate::Rect;
use std::collections::BTreeMap;

/// A node in the packed R-tree array
#[derive(Debug, Clone, Default)]
pub struct RTreeNode {
    pub bbox: Rect,
    pub hilbert_value: u32,
}

/// Packed array Hilbert R-tree for spatial indexing
pub struct RTree {
    nodes: Vec<RTreeNode>,
    first_leaf_index: usize,
}

impl RTree {
    /// Build a new R-tree from a mapping of hilbert values to bounding rectangles
    pub fn build(hilbert_to_rect: BTreeMap<u32, Rect>) -> Self {
        let num_leaves = hilbert_to_rect.len();

        if num_leaves == 0 {
            return Self {
                nodes: vec![],
                first_leaf_index: 0,
            };
        }

        // Calculate first leaf index: L rounded up to next power of 2, minus 1
        let first_leaf_index = num_leaves.next_power_of_two() - 1;
        let total_nodes = first_leaf_index + num_leaves;

        // Allocate space for all nodes
        let mut nodes = vec![Default::default(); total_nodes];

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
                let left_child_idx = 2 * i + 1;
                let right_child_idx = 2 * i + 2;

                if right_child_idx < total_nodes {
                    // Both children exist
                    let left_child = &nodes[left_child_idx];
                    let right_child = &nodes[right_child_idx];

                    nodes[i] = RTreeNode {
                        bbox: left_child.bbox.combine(&right_child.bbox),
                        hilbert_value: left_child.hilbert_value.max(right_child.hilbert_value),
                    };
                } else if left_child_idx < total_nodes {
                    // Only left child exists
                    let left_child = nodes[left_child_idx].clone();
                    nodes[i] = left_child;
                } else {
                    // No children (shouldn't happen but handle it)
                    // Leave the node as default
                }
            }
        }

        Self {
            nodes,
            first_leaf_index,
        }
    }

    /// Find the leaf node index for a given hilbert value
    fn find_leaf(&self, hilbert_value: u32) -> Option<usize> {
        if self.nodes.is_empty() {
            return None;
        }

        let mut idx = 0;

        // Walk down the tree
        while idx < self.first_leaf_index {
            let left_child_idx = 2 * idx + 1;
            let right_child_idx = 2 * idx + 2;

            if right_child_idx >= self.nodes.len() {
                // No right child, must be a leaf
                break;
            }

            let left_hilbert = self.nodes[left_child_idx].hilbert_value;

            if hilbert_value <= left_hilbert {
                idx = left_child_idx;
            } else {
                idx = right_child_idx;
            }
        }

        // Check if we found the right leaf
        if idx >= self.first_leaf_index && self.nodes[idx].hilbert_value == hilbert_value {
            Some(idx)
        } else {
            None
        }
    }

    /// Enlarge the bounding box for a given hilbert value
    /// Returns true if any change was made
    pub fn enlarge(&mut self, hilbert_value: u32, new_bbox: Rect) -> bool {
        if let Some(mut idx) = self.find_leaf(hilbert_value) {
            // Combine with existing bbox
            let combined = self.nodes[idx].bbox.combine(&new_bbox);

            if combined == self.nodes[idx].bbox {
                return false; // No change needed
            }

            self.nodes[idx].bbox = combined;

            // Propagate changes up the tree
            while idx > 0 {
                let parent_idx = (idx - 1) / 2;
                let left_child_idx = 2 * parent_idx + 1;
                let right_child_idx = 2 * parent_idx + 2;

                let new_parent_bbox = if right_child_idx < self.nodes.len() {
                    self.nodes[left_child_idx]
                        .bbox
                        .combine(&self.nodes[right_child_idx].bbox)
                } else {
                    self.nodes[left_child_idx].bbox
                };

                if new_parent_bbox == self.nodes[parent_idx].bbox {
                    break; // No change, stop propagating
                }

                self.nodes[parent_idx].bbox = new_parent_bbox;
                idx = parent_idx;
            }

            true
        } else {
            false
        }
    }

    /// Search for all leaf hilbert values whose rectangles overlap with the test rectangle
    pub fn search(&self, test_rect: &Rect) -> RTreeSearchIterator<'_> {
        RTreeSearchIterator {
            tree: self,
            test_rect: *test_rect,
            stack: if self.nodes.is_empty() {
                vec![]
            } else {
                vec![0]
            },
        }
    }

    #[cfg(test)]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
}

/// Iterator for R-tree spatial search
pub struct RTreeSearchIterator<'a> {
    tree: &'a RTree,
    test_rect: Rect,
    stack: Vec<usize>,
}

impl<'a> Iterator for RTreeSearchIterator<'a> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(idx) = self.stack.pop() {
            if idx >= self.tree.nodes.len() {
                continue;
            }

            let node = &self.tree.nodes[idx];

            // Check if this node overlaps with test rectangle
            if !node.bbox.overlaps(&self.test_rect) {
                continue;
            }

            // If this is a leaf node, return its hilbert value
            if idx >= self.tree.first_leaf_index {
                return Some(node.hilbert_value);
            }

            // Otherwise, push children to stack (right first, then left, so left is processed first)
            let left_child_idx = 2 * idx + 1;
            let right_child_idx = 2 * idx + 2;

            if right_child_idx < self.tree.nodes.len() {
                self.stack.push(right_child_idx);
            }
            if left_child_idx < self.tree.nodes.len() {
                self.stack.push(left_child_idx);
            }
        }

        None
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

        assert!(tree.node_count() > 0);
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
        assert!(changed);

        // Enlarge with same rect should not change anything
        let changed = tree.enlarge(10, Rect::new(0, 15, 0, 15));
        assert!(!changed);
    }

    #[test]
    fn test_rtree_empty() {
        let map = BTreeMap::new();
        let tree = RTree::build(map);

        assert_eq!(tree.node_count(), 0);

        let test_rect = Rect::new(0, 10, 0, 10);
        let results: Vec<u32> = tree.search(&test_rect).collect();
        assert_eq!(results.len(), 0);
    }
}

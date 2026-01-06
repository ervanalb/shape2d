/// Convert (x, y) coordinates to a Hilbert curve index
/// Input coordinates are u16 values, output is u32
// TODO(Eric): See if this can be replaced with clever bit fiddling
pub fn xy_to_hilbert(x: u16, y: u16) -> u32 {
    let mut x = x as u32;
    let mut y = y as u32;
    let mut d = 0u32;
    let n = 65536u32; // 2^16 for u16 coordinates

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
}

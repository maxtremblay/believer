/// An iterator over all edges in `self` ordered by check first.

use super::ParityCheckMatrix;

/// An iterator over all edges in `self` ordered by check first.
pub struct EdgesIter<'a> {
    active_check: usize,
    index: usize,
    check_ranges: &'a [usize],
    bit_indices: &'a [usize],
}

impl<'a> EdgesIter<'a> {
    pub(super) fn from(matrix: &'a ParityCheckMatrix) -> Self {
        Self {
            active_check: 0,
            index: 0,
            check_ranges: &matrix.check_ranges,
            bit_indices: &matrix.bit_indices,
        }
    }

    fn get_active_edge(&self) -> Option<(usize, usize)> {
        self.bit_indices
            .get(self.index)
            .map(|bit| (self.active_check, *bit))
    }

    fn go_to_next_edge(&mut self) {
        self.index += 1;
        if self.has_reached_end_of_a_check() {
            self.go_to_next_check();
        }
    }

    fn has_reached_end_of_a_check(&self) -> bool {
        self.get_end_of_active_check()
            .map(|check_end| self.index >= check_end)
            .unwrap_or(true)
    }

    fn get_end_of_active_check(&self) -> Option<usize> {
        self.check_ranges.get(self.active_check + 1).cloned()
    }

    fn go_to_next_check(&mut self) {
        self.active_check += 1;
    }

}

impl<'a> Iterator for EdgesIter<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let edge = self.get_active_edge();
        self.go_to_next_edge();
        edge
    }
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn edges_iterator() {
        let parity_check = ParityCheckMatrix::with_n_bits(3)
            .with_checks(vec![vec![0, 1], vec![1, 2]]);
        
        let mut iter = parity_check.edges_iter();

        assert_eq!(iter.next(), Some((0, 0)));
        assert_eq!(iter.next(), Some((0, 1)));
        assert_eq!(iter.next(), Some((1, 1)));
        assert_eq!(iter.next(), Some((1, 2)));
        assert_eq!(iter.next(), None);
    }
}

//! An iterator over the checks of a parity check matrix.

use super::{ParityCheckMatrix, CheckView};

/// An iterator over the checks of a parity check matrix.
/// 
/// Returns a `CheckView` at each iteration.
pub struct ChecksIter<'a> {
    matrix: &'a ParityCheckMatrix,
    active_check: usize,
}

impl<'a> ChecksIter<'a> {
    pub(super) fn from(matrix: &'a ParityCheckMatrix) -> Self {
        Self { matrix, active_check: 0 }
    }
}

impl<'a> Iterator for ChecksIter<'a> {
    type Item = CheckView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let slice = self.matrix.get_check(self.active_check);
        self.active_check += 1;
        slice
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn checks_iterator() {
        let checks = vec![vec![0, 1, 2], vec![2, 3], vec![0, 3], vec![1, 2, 3]];
        let matrix = ParityCheckMatrix::with_n_bits(4).with_checks(checks);
        let mut iter = matrix.checks_iter();

        assert_eq!(iter.next(), Some(CheckView::from_slice(&[0, 1, 2])));
        assert_eq!(iter.next(), Some(CheckView::from_slice(&[2, 3])));
        assert_eq!(iter.next(), Some(CheckView::from_slice(&[0, 3])));
        assert_eq!(iter.next(), Some(CheckView::from_slice(&[1, 2, 3])));
        assert_eq!(iter.next(), None);
    }
}
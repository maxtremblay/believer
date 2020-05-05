//! An iterator over the checks of a parity check matrix.

use super::ParityCheckMatrix;

/// An iterator over the checks of a parity check matrix.
///
/// Returns a `CheckView` at each iteration.
pub struct Checks<'a> {
    matrix: &'a ParityCheckMatrix,
    active_check: usize,
}

impl<'a> Checks<'a> {
    pub(super) fn from(matrix: &'a ParityCheckMatrix) -> Self {
        Self {
            matrix,
            active_check: 0,
        }
    }
}

impl<'a> Iterator for Checks<'a> {
    type Item = &'a [usize];

    fn next(&mut self) -> Option<Self::Item> {
        let slice = self.matrix.check(self.active_check);
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
        let matrix = ParityCheckMatrix::with_block_size(4).with_checks(checks.clone());
        let mut iter = Checks::from(&matrix);

        assert_eq!(iter.next(), Some(checks[0].as_slice()));
        assert_eq!(iter.next(), Some(checks[1].as_slice()));
        assert_eq!(iter.next(), Some(checks[2].as_slice()));
        assert_eq!(iter.next(), Some(checks[3].as_slice()));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn checks_iterator_for_empty_matrix() {
        let matrix = ParityCheckMatrix::new();
        let mut iter = Checks::from(&matrix);
        assert_eq!(iter.next(), None);
    }
}

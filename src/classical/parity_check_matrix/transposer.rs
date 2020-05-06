use super::{Check, CheckSlice, ParityCheckMatrix};

// A tool to help transpose a parity check matrix.
pub(super) struct Transposer<'a> {
    matrix: &'a ParityCheckMatrix,
    checks: Vec<Check>,
    active_check: usize,
}

impl<'a> Transposer<'a> {
    pub(super) fn from(matrix: &'a ParityCheckMatrix) -> Self {
        let checks = vec![Vec::new(); matrix.number_of_bits()];
        Self {
            matrix,
            checks,
            active_check: 0,
        }
    }

    pub(super) fn transposed(mut self) -> ParityCheckMatrix {
        self.transpose_checks();
        ParityCheckMatrix::with_block_size(self.matrix.number_of_checks()).with_checks(self.checks)
    }

    fn transpose_checks(&mut self) {
        self.matrix
            .checks()
            .for_each(|check| self.insert_bits_from(check));
    }

    fn insert_bits_from(&mut self, check: CheckSlice) {
        check.iter().for_each(|bit| self.insert(*bit));
        self.active_check += 1;
    }

    fn insert(&mut self, bit: usize) {
        self.checks[bit].push(self.active_check);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn transposition_of_empty_matrix() {
        let transposed = Transposer::from(&ParityCheckMatrix::new()).transposed();
        assert_eq!(transposed, ParityCheckMatrix::new());
    }

    #[test]
    fn transposition_of_general_matrix() {
        let matrix = ParityCheckMatrix::with_block_size(5).with_checks(vec![
            vec![0, 1, 4],
            vec![2, 3],
            vec![1, 3, 4],
            vec![0, 2],
        ]);
        let transposed = Transposer::from(&matrix).transposed();

        let expected = ParityCheckMatrix::with_block_size(4).with_checks(vec![
            vec![0, 3],
            vec![0, 2],
            vec![1, 3],
            vec![1, 2],
            vec![0, 2],
        ]);

        assert_eq!(transposed, expected);
    }

    #[test]
    fn tranposition_with_some_empty_transposed_checks() {
        let matrix = ParityCheckMatrix::with_block_size(5).with_checks(vec![
            vec![0, 1, 4],
            vec![2, 4],
            vec![0, 1, 2],
        ]);
        let transposed = Transposer::from(&matrix).transposed();

        let expected = ParityCheckMatrix::with_block_size(3).with_checks(vec![
            vec![0, 2],
            vec![0, 2],
            vec![1, 2],
            vec![0, 1],
        ]);

        assert_eq!(transposed, expected);
    }
}

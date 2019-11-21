use super::{ParityCheckMatrix, Check, CheckView};
use itertools::EitherOrBoth;
use itertools::Itertools;

pub(super) struct Concatener<'a> {
    left_matrix: &'a ParityCheckMatrix,
    right_matrix: &'a ParityCheckMatrix,
}

impl<'a> Concatener<'a> {
    pub(super) fn from(
        left_matrix: &'a ParityCheckMatrix,
        right_matrix: &'a ParityCheckMatrix,
    ) -> Self {
        Self { left_matrix, right_matrix }
    }

    pub(super) fn concat_horizontally(&self) -> ParityCheckMatrix {
        let n_bits = self.left_matrix.get_n_bits() + self.right_matrix.get_n_bits();
        let checks = self.get_checks_of_horizontal_concatenation();
        ParityCheckMatrix::with_n_bits(n_bits).with_checks(checks)
    }

    fn get_checks_of_horizontal_concatenation(&self) -> Vec<Check> {
        self.left_matrix
            .checks_iter()
            .zip_longest(self.right_matrix.checks_iter())
            .map(|checks| self.concat_horizontally_checks(checks))
            .collect()
    }

    fn concat_horizontally_checks(&self, checks: EitherOrBoth<CheckView, CheckView>) -> Check {
        match checks {
            EitherOrBoth::Both(left_check, right_check) => self.concat(left_check, right_check),
            EitherOrBoth::Left(check) => check.to_vec(),
            EitherOrBoth::Right(check) => self.pad_right_check(check)
        }
    }

    fn concat(&self, left_check: CheckView, right_check: CheckView) -> Check {
        let mut check = left_check.to_vec();
        check.append(&mut self.pad_right_check(right_check));
        check
    }

    fn pad_right_check(&self, check: CheckView) -> Check {
        check.iter().map(|bit| bit + self.left_matrix.get_n_bits()).collect()
    }

    pub(super) fn concat_diagonally(&self) -> ParityCheckMatrix {
        let n_bits = self.left_matrix.get_n_bits() + self.right_matrix.get_n_bits();
        let checks = self.get_checks_of_diagonal_concatenation();
        ParityCheckMatrix::with_n_bits(n_bits).with_checks(checks)
    }

    fn get_checks_of_diagonal_concatenation(&self) -> Vec<Check> {
        let n_checks = self.left_matrix.get_n_checks() + self.right_matrix.get_n_checks();
        let mut checks = Vec::with_capacity(n_checks);
        checks.append(&mut self.get_all_left_checks());
        checks.append(&mut self.get_all_padded_right_checks());
        checks
    }

    fn get_all_left_checks(&self) -> Vec<Check> {
        self.left_matrix.checks_iter().map(|check| check.to_vec()).collect()
    }

    fn get_all_padded_right_checks(&self) -> Vec<Check> {
        self.right_matrix.checks_iter().map(|check| self.pad_right_check(check)).collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn horizontal_concat_with_empty_matrix() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);

        let right_matrix = ParityCheckMatrix::new();

        let concatened = Concatener::from(&left_matrix, &right_matrix).concat_horizontally();

        assert_eq!(concatened, left_matrix);
    }

    #[test]
    fn horizontal_concat_from_empty_matrix() {
        let left_matrix = ParityCheckMatrix::new();

        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);


        let concatened = Concatener::from(&left_matrix, &right_matrix).concat_horizontally();

        assert_eq!(concatened, right_matrix);
    }

    #[test]
    fn horizontal_concat_with_smaller_left_matrix() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);

        let right_matrix = ParityCheckMatrix::with_n_bits(3)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![0, 2]]);

        let concatened = Concatener::from(&left_matrix, &right_matrix).concat_horizontally();
        let expected = ParityCheckMatrix::with_n_bits(7)
            .with_checks(vec![vec![0, 1, 4, 5], vec![1, 2, 3, 5, 6], vec![4, 6]]);

        assert_eq!(concatened, expected);
    }

    #[test]
    fn horizontal_concat_with_smaller_right_matrix() {
        let left_matrix = ParityCheckMatrix::with_n_bits(3)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![0, 2]]);
        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);

        let concatened = Concatener::from(&left_matrix, &right_matrix).concat_horizontally();
        let expected = ParityCheckMatrix::with_n_bits(7)
            .with_checks(vec![vec![0, 1, 3, 4], vec![1, 2, 4, 5, 6], vec![0, 2]]);

        assert_eq!(concatened, expected);
    }

    #[test]
    fn horizontal_concat_with_equal_length_matrices() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![2, 3]]);
        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3], vec![0, 2, 3]]);

        let concatened = Concatener::from(&left_matrix, &right_matrix).concat_horizontally();
        let expected = ParityCheckMatrix::with_n_bits(8)
            .with_checks(vec![vec![0, 1, 4, 5], vec![1, 2, 5, 6 ,7], vec![2, 3, 4, 6, 7]]);

        assert_eq!(concatened, expected);
    }


    #[test]
    fn diagonal_concat_with_empty_matrix() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![2, 3]]);
        let right_matrix = ParityCheckMatrix::new();

        let concatened = Concatener::from(&left_matrix, &right_matrix).concat_diagonally();

        assert_eq!(concatened, left_matrix);
    }


    #[test]
    fn diagonal_concat_from_empty_matrix() {
        let left_matrix = ParityCheckMatrix::new();

        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![2, 3]]);

        let concatened = Concatener::from(&left_matrix, &right_matrix).concat_diagonally();

        assert_eq!(concatened, right_matrix);
    }

    #[test]
    fn diagonal_concat() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![2, 3]]);
        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);

        let concatened = Concatener::from(&left_matrix, &right_matrix).concat_diagonally();
        let expected = ParityCheckMatrix::with_n_bits(8)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![2, 3], vec![4, 5], vec![5, 6, 7]]);

        assert_eq!(concatened, expected);
    }
}

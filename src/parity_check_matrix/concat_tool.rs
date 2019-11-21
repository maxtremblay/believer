use super::ParityCheckMatrix;

pub(super) struct ConcatTool {}

impl ConcatTool {
    pub(super) fn concat_diagonally_with(other: &ParityCheckMatrix) -> ParityCheckMatrix {
        unimplemented!()
    }

    pub(super) fn concat_right_with(other: &ParityCheckMatrix) -> ParityCheckMatrix {
        unimplemented!()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn diagonal_concat_with_empty_matrix() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![2, 3]]);
        let right_matrix = ParityCheckMatrix::new()

        let concatened = ConcatTool::from(left_matrix).concat_diagonally_with(right_matrix);

        assert_eq!(concatened, left_matrix);
    }


    #[test]
    fn diagonal_concat_from_empty_matrix() {
        let left_matrix = ParityCheckMatrix::new()

        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![2, 3]]);

        let concatened = ConcatTool::from(left_matrix).concat_diagonally_with(right_matrix);

        assert_eq!(concatened, right_matrix);
    }

    #[test]
    fn diagonal_concat() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![2, 3]]);
        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);

        let concatened = ConcatTool::from(left_matrix).concat_diagonally_with(right_matrix);
        let expected = ParityCheckMatrix::with_n_bits(8)
            .with_checks(vec![0, 1], vec![1, 2], vec![2, 3], vec![4, 5], vec![5, 6, 7]);

        assert_eq!(concatened, expected);
    }


    #[test]
    fn right_concat_with_empty_matrix() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);

        let right_matrix = ParityCheckMatrix::new();

        let concatened = ConcatTool::from(left_matrix).concat_right_with(right_matrix);

        assert_eq!(concatened, left_matrix);
    }

    #[test]
    fn right_concat_from_empty_matrix() {
        let left_matrix = ParityCheckMatrix::new();

        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);


        let concatened = ConcatTool::from(left_matrix).concat_right_with(right_matrix);

        assert_eq!(concatened, right_matrix);
    }

    #[test]
    fn right_concat_with_smaller_left_matrix() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);

        let right_matrix = ParityCheckMatrix::with_n_bits(3)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![0, 2]]);

        let concatened = ConcatTool::from(left_matrix).concat_right_with(right_matrix);
        let expected = ParityCheckMatrix::with_n_bits(8)
            .with_checks(vec![0, 1, 4, 5], vec![1, 2, 3, 5, 6], vec![4, 6]);

        assert_eq!(concatened, expected);
    }

    #[test]
    fn right_concat_with_smaller_right_matrix() {
        let left_matrix = ParityCheckMatrix::with_n_bits(3)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![0, 2]]);
        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3]]);

        let concatened = ConcatTool::from(left_matrix).concat_right_with(right_matrix);
        let expected = ParityCheckMatrix::with_n_bits(8)
            .with_checks(vec![0, 1, 4, 5], vec![1, 2, 5, 6 ,7], vec![0, 2]);

        assert_eq!(concatened, expected);
    }

    #[test]
    fn right_concat_with_equal_length_matrices() {
        let left_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2], vec![2, 3]]);
        let right_matrix = ParityCheckMatrix::with_n_bits(4)
            .with_checks(vec![vec![0, 1], vec![1, 2, 3], vec![0, 2, 3]]);

        let concatened = ConcatTool::from(left_matrix).concat_right_with(right_matrix);
        let expected = ParityCheckMatrix::with_n_bits(8)
            .with_checks(vec![0, 1, 4, 5], vec![1, 2, 5, 6 ,7], vec![2, 3, 4, 6, 7]);

        assert_eq!(concatened, expected);
    }
}

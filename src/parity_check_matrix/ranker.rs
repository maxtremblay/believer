use super::get_bitwise_sum;
use super::ParityCheckMatrix;

type Row = Vec<usize>;

pub(super) struct Ranker {
    n_columns: usize,
    n_rows: usize,
    rank: usize,
    active_column: usize,
    rows: Vec<Row>,
}

impl Ranker {
    pub(super) fn from_parity_check_matrix(matrix: &ParityCheckMatrix) -> Self {
        Self {
            n_columns: matrix.get_n_bits(),
            n_rows: matrix.get_n_checks(),
            rank: 0,
            active_column: 0,
            rows: matrix.checks_iter().map(|check| check.to_vec()).collect(),
        }
    }

    pub(super) fn get_rank(mut self) -> usize {
        while self.is_not_in_echelon_form() {
            self.pivot_active_column();
            self.go_to_next_column();
        }
        self.rank
    }

    fn is_not_in_echelon_form(&self) -> bool {
        self.active_column < self.n_columns
    }

    fn pivot_active_column(&mut self) {
        if let Some(pivot) = self.find_and_remove_pivot() {
            self.pivot_rows_that_start_in_active_column_with(pivot);
            self.rank += 1;
        }
    }

    fn find_and_remove_pivot(&mut self) -> Option<Row> {
        for row_index in 0..self.n_rows {
            if self.row_at_index_start_at_active_column(row_index) {
                let row = self.get_and_remove_row_at_index(row_index);
                return Some(row);
            }
        }
        None
    }

    fn row_at_index_start_at_active_column(&self, index: usize) -> bool {
        self.rows[index]
            .first()
            .map(|column| *column == self.active_column)
            .unwrap_or(false)
    }

    fn get_and_remove_row_at_index(&mut self, index: usize) -> Row {
        std::mem::replace(&mut self.rows[index], Vec::new())
    }

    fn pivot_rows_that_start_in_active_column_with(&mut self, pivot: Row) {
        for row_index in 0..self.n_rows {
            if self.row_at_index_start_at_active_column(row_index) {
                self.rows[row_index] = get_bitwise_sum(&pivot, &self.rows[row_index])
            }
        }
    }

    fn go_to_next_column(&mut self) {
        self.active_column += 1;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn rank_for_empty_matrix_is_zero() {
        let matrix = ParityCheckMatrix::new();
        let rank = Ranker::from_parity_check_matrix(&matrix).get_rank();
        assert_eq!(rank, 0);
    }

    #[test]
    fn rank_for_the_repetition_code_is_one_less_than_the_number_of_bits() {
        let checks = vec![vec![0, 1], vec![1, 2], vec![2, 3], vec![3, 4], vec![4, 0]];
        let matrix = ParityCheckMatrix::with_n_bits(5).with_checks(checks);
        let rank = Ranker::from_parity_check_matrix(&matrix).get_rank();
        assert_eq!(rank, 4);
    }

    #[test]
    fn rank_for_a_full_rank_matrix_is_the_number_of_checks() {
        let checks = vec![vec![0, 1], vec![1, 2], vec![0, 1, 3]];
        let matrix = ParityCheckMatrix::with_n_bits(4).with_checks(checks);
        let rank = Ranker::from_parity_check_matrix(&matrix).get_rank();
        assert_eq!(rank, 3);
    }

    #[test]
    fn rank_for_a_given_matrix_is_accurate() {
        let checks = vec![
            vec![0, 1, 2],
            vec![1, 2, 3],
            vec![0, 3],
            vec![3, 4, 5],
            vec![0, 4, 6],
            vec![5, 6],
        ];
        let matrix = ParityCheckMatrix::with_n_bits(7).with_checks(checks);
        let rank = Ranker::from_parity_check_matrix(&matrix).get_rank();
        assert_eq!(rank, 4);
    }
}

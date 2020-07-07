use super::ParityCheckMatrix;

pub(super) struct HyperGraphProduct<'a> {
    first_matrix: &'a ParityCheckMatrix,
    second_matrix: &'a ParityCheckMatrix,
}

impl<'a> HyperGraphProduct<'a> {
    pub(super) fn of(
        first_matrix: &'a ParityCheckMatrix,
        second_matrix: &'a ParityCheckMatrix,
    ) -> Self {
        Self {
            first_matrix,
            second_matrix,
        }
    }

    pub(super) fn compute(&self) -> (ParityCheckMatrix, ParityCheckMatrix) {
        (self.x_matrix(), self.z_matrix())
    }

    fn x_matrix(&self) -> ParityCheckMatrix {
        self.x_left_matrix()
            .get_horizontal_concat_with(&self.x_right_matrix())
    }

    fn x_left_matrix(&self) -> ParityCheckMatrix {
        self.first_matrix
            .tensor_product_with(&ParityCheckMatrix::identity_with_n_bits(
                self.second_matrix.get_n_bits(),
            ))
    }

    fn x_right_matrix(&self) -> ParityCheckMatrix {
        ParityCheckMatrix::identity_with_n_bits(self.first_matrix.get_n_checks())
            .tensor_product_with(&self.second_matrix.get_transposed_matrix())
    }

    fn z_matrix(&self) -> ParityCheckMatrix {
        self.z_left_matrix()
            .get_horizontal_concat_with(&self.z_right_matrix())
    }

    fn z_left_matrix(&self) -> ParityCheckMatrix {
        ParityCheckMatrix::identity_with_n_bits(self.first_matrix.get_n_bits())
            .tensor_product_with(&self.second_matrix)
    }

    fn z_right_matrix(&self) -> ParityCheckMatrix {
        self.first_matrix
            .get_transposed_matrix()
            .tensor_product_with(&ParityCheckMatrix::identity_with_n_bits(
                self.second_matrix.get_n_checks(),
            ))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn hyper_graph_product_identity_with_repetition_code() {
        let (x_matrix, z_matrix) =
            ParityCheckMatrix::identity_with_n_bits(2).hyper_graph_product_with(&repetition_code());
        let x_expected = ParityCheckMatrix::with_n_bits(10).with_checks(vec![
            vec![0, 6],
            vec![1, 6, 7],
            vec![2, 7],
            vec![3, 8],
            vec![4, 8, 9],
            vec![5, 9],
        ]);
        let z_expected = ParityCheckMatrix::with_n_bits(10).with_checks(vec![
            vec![0, 1, 6],
            vec![1, 2, 7],
            vec![3, 4, 8],
            vec![4, 5, 9],
        ]);
        println!("{}", x_matrix);
        assert_eq!(x_matrix, x_expected);
        assert_eq!(z_matrix, z_expected);
    }

    #[test]
    fn hyper_graph_product_repetition_code_with_identity() {
        let (x_matrix, z_matrix) =
            repetition_code().hyper_graph_product_with(&ParityCheckMatrix::identity_with_n_bits(2));
        let x_expected = ParityCheckMatrix::with_n_bits(10).with_checks(vec![
            vec![0, 2, 6],
            vec![1, 3, 7],
            vec![2, 4, 8],
            vec![3, 5, 9],
        ]);
        let z_expected = ParityCheckMatrix::with_n_bits(10).with_checks(vec![
            vec![0, 6],
            vec![1, 7],
            vec![2, 6, 8],
            vec![3, 7, 9],
            vec![4, 8],
            vec![5, 9],
        ]);
        println!("{}", x_matrix);
        assert_eq!(x_matrix, x_expected);
        assert_eq!(z_matrix, z_expected);
    }

    #[test]
    fn hyper_graph_product_repetition_code_with_repetition_code() {
        let (x_matrix, z_matrix) = repetition_code().hyper_graph_product_with(&repetition_code());
        let x_expected = ParityCheckMatrix::with_n_bits(13).with_checks(vec![
            vec![0, 3, 9],
            vec![1, 4, 9, 10],
            vec![2, 5, 10],
            vec![3, 6, 11],
            vec![4, 7, 11, 12],
            vec![5, 8, 12],
        ]);
        let z_expected = ParityCheckMatrix::with_n_bits(13).with_checks(vec![
            vec![0, 1, 9],
            vec![1, 2, 10],
            vec![3, 4, 9, 11],
            vec![4, 5, 10, 12],
            vec![6, 7, 11],
            vec![7, 8, 12],
        ]);
        println!("{}", x_matrix);
        assert_eq!(x_matrix, x_expected);
        assert_eq!(z_matrix, z_expected);
    }

    fn repetition_code() -> ParityCheckMatrix {
        ParityCheckMatrix::with_n_bits(3).with_checks(vec![vec![0, 1], vec![1, 2]])
    }
}

use super::{CheckView, ParityCheckMatrix};

pub(super) struct TensorProduct<'a> {
    left_matrix: &'a ParityCheckMatrix,
    right_matrix: &'a ParityCheckMatrix,
}

impl<'a> TensorProduct<'a> {
    pub(super) fn of(
        left_matrix: &'a ParityCheckMatrix,
        right_matrix: &'a ParityCheckMatrix,
    ) -> Self {
        Self {
            left_matrix,
            right_matrix,
        }
    }

    pub(super) fn compute(&self) -> ParityCheckMatrix {
        ParityCheckMatrix::with_n_bits(self.get_n_bits()).with_checks(self.get_checks())
    }

    fn get_n_bits(&self) -> usize {
        self.left_matrix.get_n_bits() * self.right_matrix.get_n_bits()
    }

    fn get_checks(&self) -> Vec<Vec<usize>> {
        self.left_matrix
            .checks_iter()
            .flat_map(|left_check| self.left_tensor_product(left_check.as_ref()))
            .collect()
    }

    fn left_tensor_product(&self, left_check: &[usize]) -> Vec<Vec<usize>> {
        self.right_matrix
            .checks_iter()
            .map(|right_check| self.vec_tensor_product(left_check, right_check.as_ref()))
            .collect()
    }

    fn vec_tensor_product(&self, left_check: &[usize], right_check: &[usize]) -> Vec<usize> {
        left_check
            .iter()
            .flat_map(|left_bit| {
                right_check
                    .iter()
                    .map(move |right_bit| left_bit * self.right_matrix.get_n_bits() + right_bit)
            })
            .collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn tensor_product_identity_with_repetition_code() {
        let result =
            ParityCheckMatrix::identity_with_n_bits(2).tensor_product_with(&repetition_code());
        let expected = ParityCheckMatrix::with_n_bits(6).with_checks(vec![
            vec![0, 1],
            vec![1, 2],
            vec![3, 4],
            vec![4, 5],
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn tensor_product_repetition_code_with_identity() {
        let result =
            repetition_code().tensor_product_with(&ParityCheckMatrix::identity_with_n_bits(2));
        let expected = ParityCheckMatrix::with_n_bits(6).with_checks(vec![
            vec![0, 2],
            vec![1, 3],
            vec![2, 4],
            vec![3, 5],
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn tensor_product_repetition_code_with_repetition_code() {
        let result = repetition_code().tensor_product_with(&repetition_code());
        let expected = ParityCheckMatrix::with_n_bits(9).with_checks(vec![
            vec![0, 1, 3, 4],
            vec![1, 2, 4, 5],
            vec![3, 4, 6, 7],
            vec![4, 5, 7, 8],
        ]);
        assert_eq!(result, expected);
    }

    fn repetition_code() -> ParityCheckMatrix {
        ParityCheckMatrix::with_n_bits(3).with_checks(vec![vec![0, 1], vec![1, 2]])
    }
}

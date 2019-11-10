//! NOTE: Need to check the math to be sure everything work.

use super::{Decoder, DecoderBuilder};
use crate::ErasureResult;
use crate::GF4Stabilizers;
use crate::ParityCheckMatrix;
use rand::{thread_rng, Rng};

/// Decoder for quantum erasure channel.
///
/// # Example
///
/// ```
/// # use believer::*;
/// let x_checks = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
/// let z_checks = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
/// let erasure_prob = 0.25;
/// let stabilizers = GF4Stabilizers::from_parity_check_matrices(x_checks, z_checks);
/// let decoder = QuantumErasureDecoder::new(stabilizers, erasure_prob);
/// decoder.decode(&decoder.random_error());
/// ```
pub struct QuantumErasureDecoder {
    stabilizers: GF4Stabilizers,
    erasure_prob: f64,
}

impl QuantumErasureDecoder {
    /// Creates an quantum erasure decoder.
    ///
    /// The two `ParityCheckMatrix` represent the X and Z components for the GF4 representation of
    /// a quantum code. That is, if `x_k` is the `k` check of the first `ParityCheckMatrix` and
    /// `z_k` the `k` check of the second `ParityCheckMatrix`, then `(x_kl, z_kl)` is the `l`
    /// element of the GF4 representation of the check.
    ///
    /// # Panic
    ///
    /// Panics if `erasure_prob` is not between 0.0 and 1.0.
    pub fn new(stabilizers: GF4Stabilizers, erasure_prob: f64) -> Self {
        if erasure_prob < 0.0 || erasure_prob > 1.0 {
            panic!("invalid probability");
        }
        Self {
            stabilizers,
            erasure_prob,
        }
    }
}

impl Decoder for QuantumErasureDecoder {
    // The error is the positions of the erased qubits.
    type Error = Vec<usize>;

    type Result = ErasureResult;
    type Checks = GF4Stabilizers;

    fn decode(&self, error: &Self::Error) -> ErasureResult {
        let erased_rank = self.stabilizers.keep(error).merge().rank();
        if erased_rank >= error.len() {
            let not_erased_rank = self.stabilizers.without(error).merge().rank();
            let non_commuting_overhead = erased_rank - error.len();
            let rank_sum = not_erased_rank + erased_rank - 2 * non_commuting_overhead;
            let no_error_rank = self.stabilizers.merge().rank();
            if rank_sum == no_error_rank {
                ErasureResult::Success
            } else {
                ErasureResult::Failure
            }
        } else {
            ErasureResult::Failure
        }
    }

    fn random_error(&self) -> Vec<usize> {
        let mut rng = thread_rng();
        (0..self.stabilizers.n_qubits())
            .filter(|_| rng.gen::<f64>() < self.erasure_prob)
            .collect()
    }

    fn take_checks(&mut self) -> Self::Checks {
        std::mem::replace(
            &mut self.stabilizers,
            GF4Stabilizers::from_parity_check_matrices(
                ParityCheckMatrix::empty(),
                ParityCheckMatrix::empty(),
            ),
        )
    }
}

pub struct QuantumErasureDecoderBuilder {
    erasure_prob: f64,
}

impl QuantumErasureDecoderBuilder {
    pub fn new(erasure_prob: f64) -> Self {
        Self { erasure_prob }
    }
}

impl DecoderBuilder for QuantumErasureDecoderBuilder {
    type Checks = GF4Stabilizers;
    type Decoder = QuantumErasureDecoder;

    fn build_from(&self, checks: Self::Checks) -> Self::Decoder {
        QuantumErasureDecoder::new(checks, self.erasure_prob)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::Pauli::{I, X, Z};

    #[test]
    fn steane_code() {
        let stabilizers = GF4Stabilizers::from_sparse_paulis(
            vec![
                vec![(X, 0), (X, 1), (X, 2), (X, 4)],
                vec![(X, 0), (X, 1), (X, 3), (X, 5)],
                vec![(X, 0), (X, 2), (X, 3), (X, 6)],
                vec![(Z, 0), (Z, 1), (Z, 2), (Z, 4)],
                vec![(Z, 0), (Z, 1), (Z, 3), (Z, 5)],
                vec![(Z, 0), (Z, 2), (Z, 3), (Z, 6)],
            ],
            7,
        );

        let decoder = QuantumErasureDecoder::new(stabilizers, 0.25);

        assert_eq!(decoder.decode(&vec![]), ErasureResult::Success);
        for i in 0..=6 {
            assert_eq!(decoder.decode(&vec![i]), ErasureResult::Success);
            for j in (i + 1)..=6 {
                assert_eq!(decoder.decode(&vec![i, j]), ErasureResult::Success);
            }
        }

        assert_eq!(decoder.decode(&vec![0, 1, 2]), ErasureResult::Success);
        assert_eq!(decoder.decode(&vec![2, 4, 5]), ErasureResult::Success);
        assert_eq!(decoder.decode(&vec![0, 1, 4]), ErasureResult::Success);
        assert_eq!(decoder.decode(&vec![3, 4, 5]), ErasureResult::Success);

        assert_eq!(decoder.decode(&vec![2, 4, 6]), ErasureResult::Failure);
        assert_eq!(decoder.decode(&vec![1, 2, 3]), ErasureResult::Failure);
        assert_eq!(decoder.decode(&vec![0, 3, 4]), ErasureResult::Failure);
        assert_eq!(decoder.decode(&vec![0, 2, 5]), ErasureResult::Failure);

        assert_eq!(
            decoder.decode(&vec![0, 1, 2, 3, 4, 5, 6]),
            ErasureResult::Failure
        );
    }

    #[test]
    fn shor_code_from_builder() {
        let builder = QuantumErasureDecoderBuilder::new(0.25);
        let stabilizers = GF4Stabilizers::from_sparse_paulis(
            vec![
                vec![(Z, 0), (Z, 1)],
                vec![(Z, 1), (Z, 2)],
                vec![(Z, 3), (Z, 4)],
                vec![(Z, 4), (Z, 5)],
                vec![(Z, 6), (Z, 7)],
                vec![(Z, 7), (Z, 8)],
                vec![(X, 0), (X, 1), (X, 2), (X, 3), (X, 4), (X, 5)],
                vec![(X, 3), (X, 4), (X, 5), (X, 6), (X, 7), (X, 8)],
            ],
            9,
        );

        let decoder = builder.build_from(stabilizers);

        assert_eq!(decoder.decode(&vec![]), ErasureResult::Success);
        for i in 0..9 {
            assert_eq!(decoder.decode(&vec![i]), ErasureResult::Success);
            for j in (i + 1)..9 {
                assert_eq!(decoder.decode(&vec![i, j]), ErasureResult::Success);
            }
        }
        assert_eq!(decoder.decode(&vec![0, 1, 2]), ErasureResult::Failure);
        assert_eq!(decoder.decode(&vec![0, 3, 6]), ErasureResult::Failure);
        assert_eq!(decoder.decode(&vec![0, 1, 3, 4]), ErasureResult::Success);
    }

    #[test]
    fn five_qubits_code() {
        let stabilizers = GF4Stabilizers::from_dense_paulis(
            vec![
                vec![X, Z, Z, X, I],
                vec![I, X, Z, Z, X],
                vec![X, I, X, Z, Z],
                vec![Z, X, I, X, Z],
            ],
            5,
        );

        let decoder = QuantumErasureDecoder::new(stabilizers, 0.25);

        assert_eq!(decoder.decode(&vec![]), ErasureResult::Success);
        for i in 0..5 {
            assert_eq!(decoder.decode(&vec![i]), ErasureResult::Success);
            for j in (i + 1)..5 {
                assert_eq!(decoder.decode(&vec![i, j]), ErasureResult::Success);
                for k in (j + 1)..5 {
                    assert_eq!(decoder.decode(&vec![i, j, k]), ErasureResult::Failure);
                }
            }
        }
    }
}

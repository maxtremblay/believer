//! A classical erasure decoder.

use super::{Decoder, DecodingResult};
use crate::ParityCheckMatrix;
use rand::{thread_rng, Rng};

pub struct ErasureDecoder<'a> {
    checks: &'a ParityCheckMatrix,
    erasure_prob: f64,
}

impl<'a> ErasureDecoder<'a> {
    pub fn new(checks: &'a ParityCheckMatrix, erasure_prob: f64) -> Self {
        if erasure_prob < 0.0 || erasure_prob > 1.0 {
            panic!("invalid probability");
        }
        Self {
            checks,
            erasure_prob,
        }
    }
}

impl<'a> Decoder for ErasureDecoder<'a> {
    type Error = Vec<usize>;
    type Result = ErasureResult;

    fn decode(&self, error: Self::Error) -> Self::Result {
        let erased_parity_check = self.checks.keep(&error);
        if error.len() - erased_parity_check.rank() == 0 {
            ErasureResult::Succeed
        } else {
            ErasureResult::Failed
        }
    }

    fn random_error(&self) -> Vec<usize> {
        let mut rng = thread_rng();
        (0..self.checks.n_bits())
            .filter(|_| rng.gen::<f64>() < self.erasure_prob)
            .collect()
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum ErasureResult {
    Failed,
    Succeed,
}

impl DecodingResult for ErasureResult {
    fn succeed(&self) -> bool {
        self == &Self::Succeed
    }
}

#[cfg(test)]

mod test {
    use super::*;

    #[test]
    fn repetition_code() {
        let matrix = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]]);
        let decoder = ErasureDecoder::new(&matrix, 0.2);

        assert_eq!(decoder.decode(vec![]), ErasureResult::Succeed);
        for i in 0..=2 {
            assert_eq!(decoder.decode(vec![i]), ErasureResult::Succeed);
            for j in (i + 1)..=2 {
                assert_eq!(decoder.decode(vec![i, j]), ErasureResult::Succeed);
            }
        }
        assert_eq!(decoder.decode(vec![0, 1, 2]), ErasureResult::Failed);
    }

    #[test]
    fn hamming_code() {
        let matrix =
            ParityCheckMatrix::new(vec![vec![0, 1, 2, 4], vec![0, 1, 3, 5], vec![0, 2, 3, 6]]);
        let decoder = ErasureDecoder::new(&matrix, 0.2);

        assert_eq!(decoder.decode(vec![]), ErasureResult::Succeed);
        for i in 0..=6 {
            assert_eq!(decoder.decode(vec![i]), ErasureResult::Succeed);
            for j in (i + 1)..=6 {
                assert_eq!(decoder.decode(vec![i, j]), ErasureResult::Succeed);
            }
        }
        assert_eq!(decoder.decode(vec![0, 1, 2]), ErasureResult::Succeed);
        assert_eq!(decoder.decode(vec![2, 4, 5]), ErasureResult::Succeed);
        assert_eq!(decoder.decode(vec![0, 1, 4]), ErasureResult::Succeed);
        assert_eq!(decoder.decode(vec![3, 4, 5]), ErasureResult::Succeed);

        assert_eq!(decoder.decode(vec![2, 4, 6]), ErasureResult::Failed);
        assert_eq!(decoder.decode(vec![1, 2, 3]), ErasureResult::Failed);
        assert_eq!(decoder.decode(vec![0, 3, 4]), ErasureResult::Failed);
        assert_eq!(decoder.decode(vec![0, 2, 5]), ErasureResult::Failed);

        assert_eq!(
            decoder.decode(vec![0, 1, 2, 3, 4, 5, 6]),
            ErasureResult::Failed
        );
    }
}

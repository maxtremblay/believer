//! A classical erasure decoder.

use super::{Decoder, DecodingResult};
use crate::ParityCheckMatrix;
use rand::Rng;

/// Decoder for classical erasure channel.
///
/// # Example
///
/// ```
/// # use believer::*;
/// let code = ParityCheckMatrix::with_n_bits(3).with_checks(vec![vec![0, 1], vec![1, 2]]);
/// let decoder = ErasureDecoder::with_prob(0.25).with_code(code);
/// decoder.decode(&decoder.random_error());
/// ```
pub struct ErasureDecoder {
    code: ParityCheckMatrix,
    erasure_prob: f64,
}

impl ErasureDecoder {
    /// Creates an erasure decoder.
    ///
    /// # Panic
    ///
    /// Panics if `erasure_prob` is not between 0.0 and 1.0.
    pub fn with_prob(erasure_prob: f64) -> Self {
        if erasure_prob < 0.0 || erasure_prob > 1.0 {
            panic!("invalid probability");
        }
        Self {
            code: ParityCheckMatrix::new(),
            erasure_prob,
        }
    }

    /// Set the code of self to `code` consuming it.
    pub fn with_code(mut self, code: ParityCheckMatrix) -> Self {
        self.code = code;
        self
    }
}

impl Decoder for ErasureDecoder {
    // The error is the positions of the erased bits.
    type Error = Vec<usize>;
    type Result = ErasureResult;
    type Code = ParityCheckMatrix;

    fn set_code(&mut self, code: Self::Code) {
        self.code = code;
    }

    fn take_code(&mut self) -> Self::Code {
        std::mem::replace(&mut self.code, ParityCheckMatrix::new())
    }

    // An erasure error can be corrected if there is no information in the erased submatrix. That
    // is, the number of erased bits is equal to the rank of the parity check matrix restricted to
    // the erased bit columns.
    fn decode(&self, error: &Self::Error) -> Self::Result {
        let erased_parity_check = self.code.keep(error);
        if error.len() - erased_parity_check.get_rank() == 0 {
            ErasureResult::Success
        } else {
            ErasureResult::Failure
        }
    }

    // Erase random bits with given probability.
    fn get_random_error_with_rng<R: Rng>(&self, rng: &mut R) -> Vec<usize> {
        (0..self.code.get_n_bits())
            .filter(|_| rng.gen::<f64>() < self.erasure_prob)
            .collect()
    }
}

/// An erasure decoder can either result in a `Success` when no logical bits are erased or in a
/// `Failure` when some logical bits are erased.
#[derive(Debug, PartialEq, Eq)]
pub enum ErasureResult {
    Failure,
    Success,
}

impl DecodingResult for ErasureResult {
    fn succeed(&self) -> bool {
        self == &Self::Success
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn repetition_code() {
        let code = ParityCheckMatrix::with_n_bits(3).with_checks(vec![vec![0, 1], vec![1, 2]]);
        let decoder = ErasureDecoder::with_prob(0.2).with_code(code);

        assert_eq!(decoder.decode(&vec![]), ErasureResult::Success);
        for i in 0..=2 {
            assert_eq!(decoder.decode(&vec![i]), ErasureResult::Success);
            for j in (i + 1)..=2 {
                assert_eq!(decoder.decode(&vec![i, j]), ErasureResult::Success);
            }
        }
        assert_eq!(decoder.decode(&vec![0, 1, 2]), ErasureResult::Failure);
    }

    #[test]
    fn hamming_code() {
        let code = ParityCheckMatrix::with_n_bits(7).with_checks(vec![
            vec![0, 1, 2, 4],
            vec![0, 1, 3, 5],
            vec![0, 2, 3, 6],
        ]);
        let decoder = ErasureDecoder::with_prob(0.25).with_code(code);

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
}

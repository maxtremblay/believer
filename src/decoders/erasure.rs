//! A classical erasure decoder.

use super::{Decoder, DecoderBuilder, DecodingResult};
use crate::ParityCheckMatrix;
use rand::{thread_rng, Rng};

/// Decoder for classical erasure channel.
///
/// # Example
///
/// ```
/// # use believer::*;
/// let checks = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
/// let erasure_prob = 0.25;
/// let decoder = ErasureDecoder::new(checks, erasure_prob);
/// decoder.decode(&decoder.random_error());
/// ```
pub struct ErasureDecoder {
    checks: ParityCheckMatrix,
    erasure_prob: f64,
}

impl ErasureDecoder {
    /// Creates an erasure decoder.
    ///
    /// # Panic
    ///
    /// Panics if `erasure_prob` is not between 0.0 and 1.0.
    pub fn new(checks: ParityCheckMatrix, erasure_prob: f64) -> Self {
        if erasure_prob < 0.0 || erasure_prob > 1.0 {
            panic!("invalid probability");
        }
        Self {
            checks,
            erasure_prob,
        }
    }

    /// Returns the number of bit the decoder is acting on.
    pub fn n_bits(&self) -> usize {
        self.checks.get_n_bits()
    }
}

impl Decoder for ErasureDecoder {
    // The error is the positions of the erased bits.
    type Error = Vec<usize>;
    type Result = ErasureResult;
    type Checks = ParityCheckMatrix;

    // An erasure error can be corrected if there is no information in the erased submatrix. That
    // is, the number of erased bits is equal to the rank of the parity check matrix restricted to
    // the erased bit columns.
    fn decode(&self, error: &Self::Error) -> Self::Result {
        let erased_parity_check = self.checks.keep(error);
        if error.len() - erased_parity_check.get_rank() == 0 {
            ErasureResult::Success
        } else {
            ErasureResult::Failure
        }
    }

    // Erase random bits with given probability.
    fn random_error(&self) -> Vec<usize> {
        let mut rng = thread_rng();
        (0..self.checks.get_n_bits())
            .filter(|_| rng.gen::<f64>() < self.erasure_prob)
            .collect()
    }

    fn take_checks(&mut self) -> Self::Checks {
        std::mem::replace(&mut self.checks, ParityCheckMatrix::new())
    }
}

/// Builder for ErasureDecoder
///
/// # Example
///
/// ```
/// # use believer::*;
/// let erasure_prob = 0.25;
/// let builder = ErasureDecoderBuilder::new(erasure_prob);
/// let checks = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
/// let decoder = builder.build_from(checks);
///
/// let other_checks = ParityCheckMatrix::new(vec![vec![0, 1], vec![0, 2]], 3);
/// let other_decoder = builder.build_from(other_checks);
/// ```
pub struct ErasureDecoderBuilder {
    erasure_prob: f64,
}

impl ErasureDecoderBuilder {
    /// Creates an erasure decoder builder.
    ///
    /// # Panic
    ///
    /// Panics if `erasure_prob` is not between 0.0 and 1.0.
    pub fn new(erasure_prob: f64) -> Self {
        if erasure_prob < 0.0 || erasure_prob > 1.0 {
            panic!("invalid probability");
        }
        Self { erasure_prob }
    }
}

impl DecoderBuilder for ErasureDecoderBuilder {
    type Checks = ParityCheckMatrix;
    type Decoder = ErasureDecoder;

    fn build_from(&self, checks: Self::Checks) -> Self::Decoder {
        ErasureDecoder::new(checks, self.erasure_prob)
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
        let matrix = ParityCheckMatrix::with_n_bits(3).with_checks(vec![vec![0, 1], vec![1, 2]]);
        let decoder = ErasureDecoder::new(matrix, 0.2);

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
    fn hamming_code_with_builder() {
        let decoder_builder = ErasureDecoderBuilder::new(0.25);
        let code = ParityCheckMatrix::with_n_bits(7)
            .with_checks(vec![vec![0, 1, 2, 4], vec![0, 1, 3, 5], vec![0, 2, 3, 6]]);
        let decoder = decoder_builder.build_from(code);

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

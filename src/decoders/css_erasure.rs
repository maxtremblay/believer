use super::{Decoder, DecoderBuilder, DecodingResult};
use crate::ErasureResult;
use crate::ParityCheckMatrix;
use rand::{thread_rng, Rng};

pub struct CSSErasureDecoder {
    x_checks: ParityCheckMatrix,
    z_checks: ParityCheckMatrix,
    erasure_prob: f64,
}

impl CSSErasureDecoder {
    pub fn new(checks: (ParityCheckMatrix, ParityCheckMatrix), erasure_prob: f64) -> Self {
        if erasure_prob < 0.0 || erasure_prob > 1.0 {
            panic!("invalid probability");
        }
        if checks.0.n_bits() != checks.1.n_bits() {
            panic!("different number of bits in the 2 parity check matrices")
        }
        Self {
            x_checks: checks.0,
            z_checks: checks.1,
            erasure_prob,
        }
    }
}

impl Decoder for CSSErasureDecoder {
    type Error = Vec<usize>;
    type Result = ErasureResult;
    type Checks = (ParityCheckMatrix, ParityCheckMatrix);

    fn decode(&self, error: &Self::Error) -> ErasureResult {
        let erased_x = self.x_checks.keep(error);
        let erased_z = self.x_checks.keep(error);
        let erased_as_gf4 = erased_x.right_concat(&erased_z);

        if error.len() - erased_as_gf4.rank() == 0 {
            ErasureResult::Success
        } else {
            ErasureResult::Failure
        }
    }

    fn random_error(&self) -> Vec<usize> {
        let mut rng = thread_rng();
        (0..self.x_checks.n_bits())
            .filter(|_| rng.gen::<f64>() < self.erasure_prob)
            .collect()
    }

    fn take_checks(&mut self) -> Self::Checks {
        let x_checks = std::mem::replace(&mut self.x_checks, ParityCheckMatrix::empty());
        let z_checks = std::mem::replace(&mut self.z_checks, ParityCheckMatrix::empty());
        (x_checks, z_checks)
    }
}

pub struct CSSErasureDecoderBuilder {
    erasure_prob: f64,
}

impl CSSErasureDecoderBuilder {
    pub fn new(erasure_prob: f64) -> Self {
        Self { erasure_prob }
    }
}

impl DecoderBuilder for CSSErasureDecoderBuilder {
    type Checks = (ParityCheckMatrix, ParityCheckMatrix);
    type Decoder = CSSErasureDecoder;

    fn build_from(&self, checks: Self::Checks) -> Self::Decoder {
        CSSErasureDecoder::new(checks, self.erasure_prob)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn steane_code() {
        let checks_x = ParityCheckMatrix::new(
            vec![vec![0, 1, 2, 4], vec![0, 1, 3, 5], vec![0, 2, 3, 6]],
            7,
        );
        let checks_z = checks_x.clone();

        let decoder = CSSErasureDecoder::new((checks_x, checks_z), 0.25);

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

        assert_eq!(
            decoder.decode(&vec![2, 4, 6]),
            ErasureResult::Failure
        );
        assert_eq!(
            decoder.decode(&vec![1, 2, 3]),
            ErasureResult::Failure
        );
        assert_eq!(
            decoder.decode(&vec![0, 3, 4]),
            ErasureResult::Failure
        );
        assert_eq!(
            decoder.decode(&vec![0, 2, 5]),
            ErasureResult::Failure
        );

        assert_eq!(
            decoder.decode(&vec![0, 1, 2, 3, 4, 5, 6]),
            ErasureResult::Failure
        );
    }

    #[test]
    fn shor_code_from_builder() {
        let builder = CSSErasureDecoderBuilder::new(0.25);
        let checks_x = ParityCheckMatrix::new(
            vec![
                vec![0, 1],
                vec![1, 2],
                vec![3, 4],
                vec![4, 5],
                vec![6, 7],
                vec![7, 8],
            ],
            9,
        );
        let checks_z =
            ParityCheckMatrix::new(vec![vec![0, 1, 2, 3, 4, 5], vec![3, 4, 5, 6, 7, 8]], 9);

        let decoder = builder.build_from((checks_x, checks_z));

        assert_eq!(decoder.decode(&vec![]), ErasureResult::Success);
        for i in 0..9 {
            assert_eq!(decoder.decode(&vec![i]), ErasureResult::Success);
            for j in (i + 1)..9 {
                println!("{} : {}", i, j);
                assert_eq!(decoder.decode(&vec![i, j]), ErasureResult::Success);
            }
        }
    }
}

//! A classical erasure decoder.

use crate::ParityCheckMatrix;

pub struct ErasureDecoder<'a> {
    parity_check: &'a ParityCheckMatrix,
}

impl<'a> ErasureDecoder<'a> {
    pub fn decode(&self, locations: &[usize]) -> ErasureResult {
        let erased_parity_check = self.parity_check.keep(locations);
        if locations.len() - erased_parity_check.rank() == 0  {
            ErasureResult::Success
        } else {
            ErasureResult::Failure
        }
    }

    pub fn new(parity_check: &'a ParityCheckMatrix) -> Self {
        Self { parity_check }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum ErasureResult {
    Failure,
    Success,
}

#[cfg(test)]

mod test {
    use super::*;

    #[test]
    fn repetition_code() {
        let matrix = ParityCheckMatrix::new(vec![
            vec![0, 1],
            vec![1, 2],
        ]);
        let decoder = ErasureDecoder::new(&matrix);

        assert_eq!(decoder.decode(&[]), ErasureResult::Success);
        for i in 0..=2 {
            assert_eq!(decoder.decode(&[i]), ErasureResult::Success);
            for j in (i + 1)..=2 {
                assert_eq!(decoder.decode(&[i, j]), ErasureResult::Success);
            }
        }
        assert_eq!(decoder.decode(&[0, 1, 2]), ErasureResult::Failure);
    }

    #[test]
    fn hamming_code() {
        let matrix = ParityCheckMatrix::new(vec![
            vec![0, 1, 2, 4],
            vec![0, 1, 3, 5],
            vec![0, 2, 3, 6],
        ]);
        let decoder = ErasureDecoder::new(&matrix);
        
        assert_eq!(decoder.decode(&[]), ErasureResult::Success);
        for i in 0..=6 {
            assert_eq!(decoder.decode(&[i]), ErasureResult::Success);
            for j in (i + 1)..=6 {
                assert_eq!(decoder.decode(&[i, j]), ErasureResult::Success);
            }
        }
        assert_eq!(decoder.decode(&[0, 1, 2]), ErasureResult::Success);
        assert_eq!(decoder.decode(&[2, 4, 5]), ErasureResult::Success);
        assert_eq!(decoder.decode(&[0, 1, 4]), ErasureResult::Success);
        assert_eq!(decoder.decode(&[3, 4, 5]), ErasureResult::Success);

        assert_eq!(decoder.decode(&[2, 4, 6]), ErasureResult::Failure);
        assert_eq!(decoder.decode(&[1, 2, 3]), ErasureResult::Failure);
        assert_eq!(decoder.decode(&[0, 3, 4]), ErasureResult::Failure);
        assert_eq!(decoder.decode(&[0, 2, 5]), ErasureResult::Failure);


        assert_eq!(decoder.decode(&[0, 1, 2, 3, 4, 5, 6]), ErasureResult::Failure);
    }
}
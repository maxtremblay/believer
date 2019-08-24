//! A sparse implementation of the belief propagation decoder for a binary channel.
//! 
//! The implementation is based on "Error Correction Coding: Mathematical Methods
//! and Algorithms (Chapter 15), Todd K. Moon, 2005, Wiley".

use crate::channel::BinaryChannel;
use crate::sparse_matrix::{SparseMatrix, Transposer};
use crate::ParityCheckMatrix;
use crate::GF2;

/// A `Decoder` can decode a message received from a given `channel` using a given
/// `parity_check` matrix.
///
/// # Example
///
/// ```
/// # use believer::*;
/// // A bsc with error prob 0.2.
/// let channel = channel::BinarySymmetricChannel::new(0.2);
///
/// // A 5 bits repetition code.
/// let parity_check = ParityCheckMatrix::new(vec![
///     vec![0, 1],
///     vec![1, 2],
///     vec![2, 3],
///     vec![3, 4],
/// ]);
///
/// let decoder = Decoder::new(&channel, &parity_check);
///
/// // A message with 2 errors.
/// let received_message = vec![GF2::B1, GF2::B1, GF2::B0, GF2::B0, GF2::B0];
///
/// // Should be decoded to the all zero codeword.
/// let decoded_message = decoder.decode(&received_message, 10);
/// assert_eq!(decoded_message, Some(vec![GF2::B0; 5]));
/// ```
pub struct Decoder<'a, C>
where
    C: BinaryChannel,
{
    channel: &'a C,
    parity_check: &'a ParityCheckMatrix,
    transposer: Transposer,
}

impl<'a, C: BinaryChannel> Decoder<'a, C> {
    /// Decodes a given `message`.
    ///
    /// # Panic
    ///
    /// Panics if `message` lenght doesn't correspond to `self.n_bits()`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// // A bsc with error prob 0.2.
    /// let channel = channel::BinarySymmetricChannel::new(0.2);
    ///
    /// // A 5 bits repetition code.
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     vec![0, 1],
    ///     vec![1, 2],
    ///     vec![2, 3],
    ///     vec![3, 4],
    /// ]);
    ///
    /// let decoder = Decoder::new(&channel, &parity_check);
    ///
    /// // A message with 2 errors.
    /// let received_message = vec![GF2::B1, GF2::B1, GF2::B0, GF2::B0, GF2::B0];
    ///
    /// // Should be decoded to the all zero codeword.
    /// let decoded_message = decoder.decode(&received_message, 10);
    /// assert_eq!(decoded_message, Some(vec![GF2::B0; 5]));
    /// ```
    pub fn decode(&self, message: &[GF2], max_iters: usize) -> Option<Vec<GF2>> {
        if message.len() != self.n_bits() {
            panic!("message doesn't have the right length")
        }

        let mut extrinsec_likelyhoods =
            SparseMatrix::from_parity_check(self.parity_check, vec![0.0; self.parity_check.len()]);
        let intrinsec_likelyhoods = self.channel.message_likelyhood(&message);
        let mut total_likelyhoods = intrinsec_likelyhoods.clone();
        let mut iter = 0;
        let mut decoded_message = message.to_vec();

        while !self.is_codeword(&decoded_message) && iter < max_iters {
            iter += 1;

            extrinsec_likelyhoods =
                self.check_node_update(&extrinsec_likelyhoods, &total_likelyhoods);

            total_likelyhoods =
                self.bit_node_update(&intrinsec_likelyhoods, &extrinsec_likelyhoods);

            decoded_message = total_likelyhoods
                .iter()
                .map(|l| if l > &0.0 { GF2::B1 } else { GF2::B0 })
                .collect();
        }

        if self.is_codeword(&decoded_message) {
            Some(decoded_message)
        } else {
            None
        }
    }

    pub fn new(channel: &'a C, parity_check: &'a ParityCheckMatrix) -> Self {
        Self {
            channel,
            parity_check,
            transposer: Transposer::new(parity_check),
        }
    }

    pub fn n_bits(&self) -> usize {
        self.parity_check.n_bits()
    }

    pub fn n_checks(&self) -> usize {
        self.parity_check.n_checks()
    }

    fn bit_node_update(
        &self,
        intrinsec_likelyhoods: &[f64],
        extrinsec_likelyhoods: &SparseMatrix,
    ) -> Vec<f64> {
        self.transposer
            .transpose(extrinsec_likelyhoods)
            .rows_iter()
            .map(|row| row.map(|(val, _)| val).sum())
            .zip(intrinsec_likelyhoods)
            .map(|(ext, int): (f64, &f64)| ext + *int)
            .collect()
    }

    fn check_node_update(
        &self,
        extrinsec_likelyhoods: &SparseMatrix,
        total_likelyhoods: &[f64],
    ) -> SparseMatrix<'a> {
        let likelyhood_diff_products: Vec<f64> = extrinsec_likelyhoods
            .rows_iter()
            .map(|row| {
                row.map(|(val, pos): (&f64, &usize)| {
                    ((*val - total_likelyhoods[*pos]) / 2.0).tanh()
                })
                .product::<f64>()
            })
            .collect();

        let updated_values = self
            .parity_check
            .positions_iter()
            .map(|(row, col)| {
                extrinsec_likelyhoods
                    .get(row, col)
                    .map(|val| ((val - total_likelyhoods[col]) / 2.0).tanh())
                    .map(|denominator| -2.0 * (likelyhood_diff_products[row] / denominator).atanh())
                    .unwrap_or(0.0)
            })
            .collect();

        SparseMatrix::from_parity_check(self.parity_check, updated_values)
    }

    // NOTE : Maybe implement parity_check.is_codeword(...).
    // It could be faster if the first 1 is near the beggining of message.
    fn is_codeword(&self, message: &[GF2]) -> bool {
        self.parity_check.dot(message).iter().all(|x| x == &GF2::B0)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::channel::BinarySymmetricChannel;

    #[test]
    fn general_usage() {
        // Tests with a 3 bits repetition code over a
        // bsc with error probability of 0.2.s
        let channel = BinarySymmetricChannel::new(0.2);
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]]);
        let decoder = Decoder::new(&channel, &parity_check);

        let decoded_message = decoder.decode(&vec![GF2::B0, GF2::B0, GF2::B1], 10);
        assert_eq!(decoded_message, Some(vec![GF2::B0; 3]));
    }
}

//! A sparse implementation of the belief propagation decoder for a binary channel.
//!
//! The implementation is based on "Error Correction Coding: Mathematical Methods
//! and Algorithms (Chapter 15), Todd K. Moon, 2005, Wiley".

use super::{Decoder, DecodingResult};
use crate::channel::BinaryChannel;
use crate::sparse_matrix::{SparseMatrix, Transposer};
use crate::ParityCheckMatrix;
use crate::GF2;

/// A `BPDecoder` can decode a message received from a given `channel` using a given
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
/// let parity_check = ParityCheckMatrix::new(
///     vec![
///         vec![0, 1],
///         vec![1, 2],
///         vec![2, 3],
///         vec![3, 4],
///     ],
///     5
/// );
///
/// let decoder = BPDecoder::new(&channel, &parity_check, 10);
///
/// // A message with 2 errors.
/// let received_message = vec![GF2::B1, GF2::B1, GF2::B0, GF2::B0, GF2::B0];
///
/// // Should be decoded to the all zero codeword.
/// let decoded_message = decoder.decode(&received_message);
/// assert_eq!(decoded_message, BPResult::Codeword(vec![GF2::B0; 5]));
/// ```
pub struct BPDecoder<'a, C>
where
    C: BinaryChannel,
{
    channel: &'a C,
    parity_check: &'a ParityCheckMatrix,
    transposer: Transposer,
    max_iters: usize,
}

impl<'a, C: BinaryChannel> BPDecoder<'a, C> {
    // *
    // Public methods
    // *

    /// Returns the number of bits in `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let channel = channel::BinarySymmetricChannel::new(0.2);
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///         vec![2, 3],
    ///         vec![3, 4],
    ///     ],
    ///     5
    /// );
    /// let decoder = BPDecoder::new(&channel, &parity_check, 3);
    ///
    /// assert_eq!(decoder.n_bits(), 5);
    /// ```
    pub fn n_bits(&self) -> usize {
        self.parity_check.n_bits()
    }

    /// Returns the number of checks in `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let channel = channel::BinarySymmetricChannel::new(0.2);
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///         vec![2, 3],
    ///         vec![3, 4],
    ///     ],
    ///     5
    /// );
    /// let decoder = BPDecoder::new(&channel, &parity_check, 3);
    ///
    /// assert_eq!(decoder.n_checks(), 4);
    /// ```
    pub fn n_checks(&self) -> usize {
        self.parity_check.n_checks()
    }

    /// Creates a new decoder from references to a `channel` and a `parity_check` matrix.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// // A bsc with error prob 0.2.
    /// let channel = channel::BinarySymmetricChannel::new(0.2);
    ///
    /// // A 5 bits repetition code.
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///         vec![2, 3],
    ///         vec![3, 4],
    ///     ],
    ///     5
    /// );
    /// let max_iters = 10;
    ///
    /// // The decoder
    /// let decoder = BPDecoder::new(&channel, &parity_check, max_iters);
    /// ```
    pub fn new(channel: &'a C, parity_check: &'a ParityCheckMatrix, max_iters: usize) -> Self {
        Self {
            channel,
            parity_check,
            transposer: Transposer::new(parity_check),
            max_iters,
        }
    }

    // *
    // Private methods
    // *

    // Inits all the likelyhood (total, intrinsec and extrinsec).
    fn init_likelyhoods(&self, message: &[GF2]) -> Likelyhoods<C> {
        let intrinsec = self.channel.message_likelyhood(message);
        let total = intrinsec.clone();
        let extrinsec =
            SparseMatrix::from_parity_check(self.parity_check, vec![0.0; self.parity_check.len()]);
        Likelyhoods {
            decoder: &self,
            total,
            intrinsec,
            extrinsec,
        }
    }
}
/*
impl<C> Decoder for BPDecoder<C>
where
    C: BinaryChannel,
{
    type Error = Vec<GF2>;
    type Result = BPResult;

    /// Decodes a given `message` doing at most `max_iters` iterations. Returns
    /// `Some(codeword)` if the decoder converge to a solution of `None` if
    /// it didn't find a solution in the given the number of iterations.
    ///
    /// # Panic
    ///
    /// Panics if `message` lenght doesn't correspond to `self.n_bits()`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let channel = channel::BinarySymmetricChannel::new(0.2);
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///         vec![2, 3],
    ///     ],
    ///     4
    /// );
    /// let decoder = BPDecoder::new(&channel, &parity_check, 10);
    ///
    /// // Should be able to decode this message to the 1 codeword.
    /// let easy_message = vec![GF2::B0, GF2::B1, GF2::B1, GF2::B1];
    /// let easy_decoded = decoder.decode(&easy_message);
    /// assert_eq!(easy_decoded, BPResult::Codeword(vec![GF2::B1; 4]));
    ///
    /// // Should get stuck while decoding this message.
    /// let impossible_message = vec![GF2::B0, GF2::B0, GF2::B1, GF2::B1];
    /// let impossible_decoded = decoder.decode(&impossible_message);
    /// assert_eq!(impossible_decoded, BPResult::GotStuck);
    /// ```
    fn decode(&self, message: &Self::Error) -> Self::Result {
        if message.len() != self.n_bits() {
            panic!("message doesn't have the right length")
        }

        let mut likelyhoods = self.init_likelyhoods(message);
        let mut iter = 0;
        let mut result = None;

        while result.is_none() {
            iter += 1;

            likelyhoods.check_node_update();
            likelyhoods.bit_node_update();

            if iter == self.max_iters {
                result = Some(BPResult::ReachedMaxIter);
            } else if likelyhoods.is_stuck() {
                result = Some(BPResult::GotStuck);
            } else {
                let m: Vec<GF2> = likelyhoods.message();
                if self.parity_check.has_codeword(&m) {
                    result = Some(BPResult::Codeword(m));
                }
            }
        }

        result.unwrap()
    }

    fn random_error(&self) -> Self::Error {
        self.channel.sample_uniform(GF2::B0, self.n_bits())
    }
}
*/
// ****************************
// Public utilitary constructs
// ****************************

/// The possible results of the decoding process.
///
/// `Codeword(c)`: The algorithm converge to the codeword `c`.
/// `GotStuck`: Some codewords are equally likely and the algorithm can't converge.
/// `ReachedMaxIter`: Was not able to find a valid codeword in the given number of iterations.
#[derive(Debug, PartialEq, Eq)]
pub enum BPResult {
    Codeword(Vec<GF2>),
    GotStuck,
    ReachedMaxIter,
}

impl DecodingResult for BPResult {
    fn succeed(&self) -> bool {
        if let Self::Codeword(c) = self {
            c.iter().all(|b| b == &GF2::B0)
        } else {
            false
        }
    }
}

// ****************************
// Private utilitary constructs
// ****************************

struct Likelyhoods<'a, C>
where
    C: BinaryChannel,
{
    decoder: &'a BPDecoder<'a, C>,
    total: Vec<f64>,
    intrinsec: Vec<f64>,
    extrinsec: SparseMatrix<'a>,
}

impl<'a, C: BinaryChannel> Likelyhoods<'a, C> {
    fn bit_node_update(&mut self) {
        self.total = self
            .decoder
            .transposer
            .transpose(&self.extrinsec)
            .rows_iter()
            .map(|row| row.map(|(val, _)| val).sum())
            .zip(&self.intrinsec)
            .map(|(ext, int): (f64, &f64)| ext + int)
            .collect()
    }

    fn check_node_update(&mut self) {
        let updated_values = self
            .decoder
            .parity_check
            .positions_iter()
            .map(|(row, col)| {
                self.extrinsec
                    .row_slice(row)
                    .map(|slice| {
                        -2.0 * slice
                            .filter(|(_, c)| c != &&col)
                            .map(|(val, c)| ((val - self.total[*c]) / 2.0).tanh())
                            .product::<f64>()
                            .atanh()
                    })
                    .unwrap_or(0.0)
            })
            .collect();

        self.extrinsec = SparseMatrix::from_parity_check(self.decoder.parity_check, updated_values);
    }

    fn is_stuck(&self) -> bool {
        self.total.iter().all(|l| l.abs() < 1e-12)
    }

    fn message(&self) -> Vec<GF2> {
        self.total
            .iter()
            .map(|l| if l > &0.0 { GF2::B1 } else { GF2::B0 })
            .collect()
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
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
        let decoder = BPDecoder::new(&channel, &parity_check, 5);

        // Should decode 1 error.
        let decoded_message = decoder.decode(&vec![GF2::B0, GF2::B0, GF2::B1]);
        assert_eq!(decoded_message, BPResult::Codeword(vec![GF2::B0; 3]));

        let decoded_message = decoder.decode(&vec![GF2::B1, GF2::B0, GF2::B1]);
        assert_eq!(decoded_message, BPResult::Codeword(vec![GF2::B1; 3]));
    }

    #[test]
    fn get_stuck_in_cycle() {
        let channel = BinarySymmetricChannel::new(0.2);
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1]], 2);
        let decoder = BPDecoder::new(&channel, &parity_check, 10);

        assert_eq!(decoder.decode(&vec![GF2::B0, GF2::B1]), BPResult::GotStuck);
    }
}

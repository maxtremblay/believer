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
/// assert_eq!(decoded_message, DecodingResult::Codeword(vec![GF2::B0; 5]));
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
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     vec![0, 1],
    ///     vec![1, 2],
    ///     vec![2, 3],
    /// ]);
    /// let decoder = Decoder::new(&channel, &parity_check);
    ///
    /// // Should be able to decode this message to the 1 codeword.
    /// let easy_message = vec![GF2::B0, GF2::B1, GF2::B1, GF2::B1];
    /// let easy_decoded = decoder.decode(&easy_message, 10);
    /// assert_eq!(easy_decoded, DecodingResult::Codeword(vec![GF2::B1; 4]));
    ///
    /// // Should get stuck while decoding this message.
    /// let impossible_message = vec![GF2::B0, GF2::B0, GF2::B1, GF2::B1];
    /// let impossible_decoded = decoder.decode(&impossible_message, 4);
    /// assert_eq!(impossible_decoded, DecodingResult::GotStuck);
    /// ```
    pub fn decode(&self, message: &[GF2], max_iters: usize) -> DecodingResult {
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

            if iter == max_iters {
                result = Some(DecodingResult::ReachedMaxIter);
            } else if likelyhoods.is_stuck() {
                result = Some(DecodingResult::GotStuck);
            } else {
                let m: Vec<GF2> = likelyhoods.message();
                if self.is_codeword(&m) {
                    result = Some(DecodingResult::Codeword(m));
                }
            }
        }

        result.unwrap()
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
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     vec![0, 1],
    ///     vec![1, 2],
    ///     vec![2, 3],
    ///     vec![3, 4],
    /// ]);
    ///
    /// // The decoder
    /// let decoder = Decoder::new(&channel, &parity_check);
    /// ```
    pub fn new(channel: &'a C, parity_check: &'a ParityCheckMatrix) -> Self {
        Self {
            channel,
            parity_check,
            transposer: Transposer::new(parity_check),
        }
    }

    /// Returns the number of bits in `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let channel = channel::BinarySymmetricChannel::new(0.2);
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     vec![0, 1],
    ///     vec![1, 2],
    ///     vec![2, 3],
    ///     vec![3, 4],
    /// ]);
    /// let decoder = Decoder::new(&channel, &parity_check);
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
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     vec![0, 1],
    ///     vec![1, 2],
    ///     vec![2, 3],
    ///     vec![3, 4],
    /// ]);
    /// let decoder = Decoder::new(&channel, &parity_check);
    ///
    /// assert_eq!(decoder.n_checks(), 4);
    /// ```
    pub fn n_checks(&self) -> usize {
        self.parity_check.n_checks()
    }

    // TO DO
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

    // NOTE : Maybe implement parity_check.is_codeword(...).
    // It could be faster if the first 1 is near the beggining of message.

    /// Checks if a `message` is a valid codeword.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let channel = channel::BinarySymmetricChannel::new(0.2);
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     vec![0, 1],
    ///     vec![1, 2],
    /// ]);
    /// let decoder = Decoder::new(&channel, &parity_check);
    ///
    /// assert_eq!(decoder.is_codeword(&vec![GF2::B0; 3]), true);
    /// assert_eq!(decoder.is_codeword(&vec![GF2::B0, GF2::B1, GF2::B0]), false);
    /// ```
    pub fn is_codeword(&self, message: &[GF2]) -> bool {
        self.parity_check.dot(message).iter().all(|x| x == &GF2::B0)
    }
}

/// The possible results of the decoding process.
///
/// `Codeword(c)`: The algorithm converge to the codeword `c`.
/// `GotStuck`: Some codewords are equally likely and the algorithm can't converge.
/// `ReachedMaxIter`: Was not able to find a valid codeword in the given number of iterations.
#[derive(Debug, PartialEq, Eq)]
pub enum DecodingResult {
    Codeword(Vec<GF2>),
    GotStuck,
    ReachedMaxIter,
}

struct Likelyhoods<'a, C>
where
    C: BinaryChannel,
{
    decoder: &'a Decoder<'a, C>,
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
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]]);
        let decoder = Decoder::new(&channel, &parity_check);

        // Should decode 1 error.
        let decoded_message = decoder.decode(&vec![GF2::B0, GF2::B0, GF2::B1], 5);
        assert_eq!(decoded_message, DecodingResult::Codeword(vec![GF2::B0; 3]));

        let decoded_message = decoder.decode(&vec![GF2::B1, GF2::B0, GF2::B1], 5);
        assert_eq!(decoded_message, DecodingResult::Codeword(vec![GF2::B1; 3]));
    }

    #[test]
    fn get_stuck_in_cycle() {
        let channel = BinarySymmetricChannel::new(0.2);
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1]]);
        let decoder = Decoder::new(&channel, &parity_check);

        assert_eq!(
            decoder.decode(&vec![GF2::B0, GF2::B1], 10),
            DecodingResult::GotStuck
        );
    }
}

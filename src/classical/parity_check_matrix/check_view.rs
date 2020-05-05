//! A view over a check of a parity check matrix.
//!
//! # Example
//!
//! ```
//! # use believer::ParityCheckMatrix;
//! let all_checks = vec![
//!     vec![0, 1],
//!     vec![0, 3],
//!     vec![1, 2, 3],
//! ];
//! let matrix = ParityCheckMatrix::with_n_bits(4).with_checks(all_checks);
//!
//! let second_check = matrix.get_check(2);
//! ```

/// A view over a check of a parity check matrix.
#[derive(Debug, PartialEq)]
pub struct CheckView<'a> {
    bits: &'a [usize],
}

impl<'a> CheckView<'a> {
    // Useful to create a check from a parity check matrix.
    pub(crate) fn from_slice(slice: &'a [usize]) -> Self {
        Self { bits: slice }
    }

    /// An iterator over the bits in `self`.
    pub fn iter(&self) -> std::slice::Iter<usize> {
        self.bits.iter()
    }

    /*
        /// Returns the syndrome of a given `message`. That is, returns the dot product between
        /// `self` and `message`.
        ///
        /// # Example
        ///
        /// ```
        /// # use believer::{GF2, ParityCheckMatrix};
        /// let parity_check = ParityCheckMatrix::with_n_bits(7)
        ///     .with_checks(vec![
        ///         vec![0, 1, 2, 4],
        ///         vec![0, 1, 3, 5],
        ///         vec![0, 2, 3, 6],
        ///     ]);
        /// let other = vec![GF2::B0, GF2::B1, GF2::B0, GF2::B1, GF2::B0, GF2::B1, GF2::B0];
        ///
        /// assert_eq!(parity_check.get_check(0).unwrap().compute_syndrome(&other), GF2::B1);
        /// assert_eq!(parity_check.get_check(1).unwrap().compute_syndrome(&other), GF2::B1);
        /// assert_eq!(parity_check.get_check(2).unwrap().compute_syndrome(&other), GF2::B1);
        /// ```
        pub fn compute_syndrome(&self, message: &[GF2]) -> GF2 {
            self.iter().fold(GF2::B0, |acc, &bit| {
                message.get(bit).map(|&value| acc + value).unwrap_or(acc)
            })
        }

    /// Returns `true` if the syndrome of `message` is `GF2::B1`.
    pub fn has_non_zero_syndrome(&self, message: &[GF2]) -> bool {
        self.compute_syndrome(message) == GF2::B1
    }

    /// Returns `true` if the syndrome of `message` is `GF2::B0`.
    pub fn has_zero_syndrome(&self, message: &[GF2]) -> bool {
        self.compute_syndrome(message) == GF2::B0
    }
    */

    /// Returns the number of bits connected to `self`.
    pub fn degree(&self) -> usize {
        self.bits.len()
    }

    /*
    /// Returns the spread of `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::ParityCheckMatrix;
    /// let parity_check = ParityCheckMatrix::with_n_bits(6).with_checks(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2, 5],
    ///         vec![3],
    ///         vec![0, 4],
    ///     ]
    /// );
    ///
    /// assert_eq!(parity_check.get_check(0).unwrap().spread(), 2);
    /// assert_eq!(parity_check.get_check(1).unwrap().spread(), 5);
    /// assert_eq!(parity_check.get_check(2).unwrap().spread(), 1);
    /// assert_eq!(parity_check.get_check(3).unwrap().spread(), 5);
    /// ```
    pub fn spread(&self) -> usize {
        self.max() - self.min() + 1 // Well define (>= 0), because max >= min.
    }

    /// Returns the bit connected to `self` with maximal index.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::with_n_bits(6).with_checks(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2, 5],
    ///         vec![3],
    ///         vec![0, 4],
    ///     ]
    /// );
    ///
    /// assert_eq!(parity_check.get_check(0).unwrap().max(), 1);
    /// assert_eq!(parity_check.get_check(1).unwrap().max(), 5);
    /// assert_eq!(parity_check.get_check(2).unwrap().max(), 3);
    /// assert_eq!(parity_check.get_check(3).unwrap().max(), 4);
    /// ```
    pub fn max(&self) -> usize {
        *self.bits.last().unwrap() // Check is always sorted and non empty.
    }

    /// Returns the bit connected to `self` with minimal index if `self`
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::with_n_bits(6).with_checks(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2, 5],
    ///         vec![3],
    ///         vec![0, 4],
    ///     ]
    /// );
    ///
    /// assert_eq!(parity_check.get_check(0).unwrap().min(), 0);
    /// assert_eq!(parity_check.get_check(1).unwrap().min(), 1);
    /// assert_eq!(parity_check.get_check(2).unwrap().min(), 3);
    /// assert_eq!(parity_check.get_check(3).unwrap().min(), 0);
    /// ```
    pub fn min(&self) -> usize {
        *self.bits.first().unwrap() // Check is always sorted and non empty.
    }
    */
    /// Copies `self` into a new `Vec`.
    pub fn to_vec(&self) -> Vec<usize> {
        self.as_ref().to_vec()
    }
}

impl<'a> AsRef<[usize]> for CheckView<'a> {
    /// Returns the indices of the bits connected to `self` as a slice.
    fn as_ref(&self) -> &[usize] {
        self.bits
    }
}

use crate::GF2;

/// A view over a check of a parity check matrix.
#[derive(Debug, PartialEq)]
pub struct Check<'a> {
    bits: &'a [usize],
}

impl<'a> Check<'a> {
    pub(crate) fn from_slice(slice: &'a [usize]) -> Self {
        Self { bits: slice }
    }

    /// An iterator over the bits in `self`.
    pub fn iter(&self) -> std::slice::Iter<usize> {
        self.bits.iter()
    }

    /// Returns the syndrome of a given `message`. That is, returns the dot product between
    /// `self` and `message`.
    /// 
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1, 2, 4],
    ///         vec![0, 1, 3, 5],
    ///         vec![0, 2, 3, 6],
    ///     ],
    ///     7
    /// );
    /// let other = vec![GF2::B0, GF2::B1, GF2::B0, GF2::B1, GF2::B0, GF2::B1, GF2::B0];
    ///
    /// assert_eq!(parity_check.check(0).unwrap().compute_syndrome(&other), GF2::B1);
    /// assert_eq!(parity_check.check(1).unwrap().compute_syndrome(&other), GF2::B1);
    /// assert_eq!(parity_check.check(2).unwrap().compute_syndrome(&other), GF2::B1);
    /// ```
    pub fn compute_syndrome(&self, message: &[GF2]) -> GF2 {
        self.iter().fold(GF2::B0, |acc, &bit| {
            message
                .get(bit)
                .map(|&value| acc + value)
                .unwrap_or(acc)
        })
    }

    pub fn has_non_zero_syndrome(&self, message: &[GF2]) -> bool {
        self.compute_syndrome(message) == GF2::B1
    }

    pub fn has_zero_syndrome(&self, message: &[GF2]) -> bool {
        self.compute_syndrome(message) == GF2::B0
    }

    /// Returns the number of bits connected to `self`.
    pub fn get_n_bits(&self) -> usize {
        self.bits.len()
    }

    /// Returns the spread of `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2, 5],
    ///         vec![3],
    ///         vec![0, 4],
    ///         vec![],
    ///     ],
    ///     6
    /// );
    ///
    /// assert_eq!(parity_check.check(0).unwrap().spread(), 2);
    /// assert_eq!(parity_check.check(1).unwrap().spread(), 5);
    /// assert_eq!(parity_check.check(2).unwrap().spread(), 1);
    /// assert_eq!(parity_check.check(3).unwrap().spread(), 5);
    /// ```
    pub fn spread(&self) -> usize {
        self.min()
            .and_then(|min| self.max().map(|max| max - min + 1))
            .unwrap_or(0)
    }

    pub fn max(&self) -> Option<&usize> {
        self.bits.last() // Check is always sorted
    }

    pub fn min(&self) -> Option<&usize> {
        self.bits.first() // Check is always sorted
    }
}


impl<'a> AsRef<[usize]> for Check<'a> {
    fn as_ref(&self) -> &[usize] {
        self.bits
    } 
}

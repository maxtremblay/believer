//! A sparse implementation of a parity check matrix.

use crate::GF2;

// *******************
// Parity Check Matrix
// *******************

/// A sparse parity check matrix.
#[derive(Debug, PartialEq)]
pub struct ParityCheckMatrix {
    check_ranges: Vec<usize>,
    bit_indices: Vec<usize>,
    n_bits: usize,
}

impl ParityCheckMatrix {
    // *
    // Public methods
    // *

    /// Returns an iterator that yields a slice for each check of `self`.
    ///
    /// # Example
    /// ```
    /// # use::believer::*;
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///     ],
    ///     3
    /// );
    /// let mut iter = parity_check.checks_iter();
    ///
    /// assert_eq!(iter.next(), parity_check.check(0));
    /// assert_eq!(iter.next(), parity_check.check(1));
    /// assert_eq!(iter.next(), None);
    ///
    /// ```
    pub fn checks_iter(&self) -> ChecksIter {
        ChecksIter {
            matrix: &self,
            active_check: 0,
        }
    }

    /// Returns `Some` slice of the given `check` in `self`. Returns `None` if
    /// `check` is out of bound.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///     ],
    ///     3
    /// );
    /// let slice = parity_check.check(0).unwrap();
    /// let vector = vec![GF2::B1, GF2::B1, GF2::B0];
    ///
    /// assert_eq!(slice.dot(&vector), GF2::B0);
    /// ```
    pub fn check(&self, check: usize) -> Option<Check> {
        self.check_ranges.get(check).and_then(|&check_start| {
            self.check_ranges.get(check + 1).map(|&check_end| Check {
                positions: &self.bit_indices[check_start..check_end],
            })
        })
    }

    /// Computes the dot product between `self` and a binary vector.
    ///
    /// # Example
    ///
    /// ```
    /// # use::believer::*;
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///     ],
    ///     3
    /// );
    /// let vector = vec![GF2::B0, GF2::B1, GF2::B1];
    ///
    /// assert_eq!(parity_check.dot(&vector), vec![GF2::B1, GF2::B0]);
    /// ```
    pub fn dot(&self, vector: &[GF2]) -> Vec<GF2> {
        self.checks_iter().map(|check| check.dot(vector)).collect()
    }

    /// Checks if a given `message` is a codeword of `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///     ],
    ///     3
    /// );
    /// let message = vec![GF2::B0, GF2::B1, GF2::B1];
    /// let codeword = vec![GF2::B0; 3];
    ///
    /// assert_eq!(parity_check.has_codeword(&message), false);
    /// assert_eq!(parity_check.has_codeword(&codeword), true);
    /// ```
    pub fn has_codeword(&self, message: &[GF2]) -> bool {
        self.checks_iter()
            .all(|check| check.dot(message) == GF2::B0)
    }

    /// Returns the extra spread per check of `self`. If all checks are locals (acts on
    /// consecutive bits), then the extra spread is 0. The extra spread of a check is 
    /// the diffrence between its spread and length.
    ///  
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
    /// assert_eq!(parity_check.extra_spread_per_check(), vec![0, 2, 0, 3, 0]);
    /// ```
    pub fn extra_spread_per_check(&self) -> Vec<usize> {
        self
            .checks_iter()
            .map(|check| check.spread() - check.len())
            .collect()
    }

    /// Returns `true` is there is no value in `self`.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of 1 in `self`.
    pub fn len(&self) -> usize {
        self.bit_indices.len()
    }

    /// Returns the number of bits in `self`.
    pub fn n_bits(&self) -> usize {
        self.n_bits
    }

    /// Returns the number of checks in `self`.
    pub fn n_checks(&self) -> usize {
        self.check_ranges().len() - 1
    }

    /// Creates a new `ParityCheckMatrix` from a list of `checks` where
    /// each check is a list of the bits connected to that check and the
    /// number of bits.
    ///
    /// # Example
    ///
    /// ```
    /// # use::believer::ParityCheckMatrix;
    /// // The parity check matrix of a 3 bits repetition code.
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///     ],
    ///     3
    /// );
    /// ```
    pub fn new(mut checks: Vec<Vec<usize>>, n_bits: usize) -> Self {
        let mut bit_indices = Vec::new();
        let mut check_ranges = Vec::with_capacity(checks.len() + 1);
        check_ranges.push(0);

        let mut n_elements = 0;
        for check in checks.iter_mut() {
            n_elements += check.len();
            check_ranges.push(n_elements);
            check.sort();
            bit_indices.append(check);
        }

        Self {
            check_ranges,
            bit_indices,
            n_bits,
        }
    }

    /// Returns an iterators over all positions in `self` where the value
    /// is 1. Positions are ordered by check first.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///     ],
    ///     3
    /// );
    /// let mut iter = parity_check.positions_iter();
    ///
    /// assert_eq!(iter.next(), Some((0, 0)));
    /// assert_eq!(iter.next(), Some((0, 1)));
    /// assert_eq!(iter.next(), Some((1, 1)));
    /// assert_eq!(iter.next(), Some((1, 2)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn positions_iter(&self) -> PositionsIter {
        PositionsIter {
            active_check: 0,
            index: 0,
            check_ranges: &self.check_ranges,
            bit_indices: &self.bit_indices,
        }
    }

    /// Computes the rank of `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///         vec![0, 2],
    ///     ],
    ///     3
    /// );
    /// assert_eq!(parity_check.rank(), 2);
    /// ```
    pub fn rank(&self) -> usize {
        Ranker::new(self).rank()
    }

    // *
    // Private methods
    // *

    // Returns a reference to `self.check_ranges`.
    pub(crate) fn check_ranges(&self) -> &[usize] {
        &self.check_ranges
    }

    // Returns a reference to `self.bit_indices
    pub(crate) fn bit_indices(&self) -> &[usize] {
        &self.bit_indices
    }

    pub(crate) fn keep(&self, bits: &[usize]) -> Self {
        let checks = self
            .checks_iter()
            .map(|check| {
                check
                    .positions()
                    .iter()
                    .filter(|&bit| bits.iter().any(|b| b == bit))
                    .cloned()
                    .collect()
            })
            .collect();
        Self::new(checks, self.len())
    }
}

impl std::fmt::Display for ParityCheckMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for check in self.checks_iter() {
            write!(f, "{:?}\n", check.positions())?;
        }
        Ok(())
    }
}


// ************************
// Public utilitary structs
// ************************

/// Iterator over the checks of a parity check matrix.
pub struct ChecksIter<'a> {
    matrix: &'a ParityCheckMatrix,
    active_check: usize,
}

impl<'a> Iterator for ChecksIter<'a> {
    type Item = Check<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let slice = self.matrix.check(self.active_check);
        self.active_check += 1;
        slice
    }
}

/// A wrapper over a check. It contains the positions of the bits in the check.
#[derive(Debug, PartialEq)]
pub struct Check<'a> {
    positions: &'a [usize],
}

impl<'a> Check<'a> {
    /// Returns the dot product between `self` and `other`. In that case, it corresponds to the
    /// parity of the overlap.
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
    /// assert_eq!(parity_check.check(0).unwrap().dot(&other), GF2::B1);
    /// assert_eq!(parity_check.check(1).unwrap().dot(&other), GF2::B1);
    /// assert_eq!(parity_check.check(2).unwrap().dot(&other), GF2::B1);
    /// ```
    pub fn dot(&self, other: &[GF2]) -> GF2 {
        let mut total = GF2::B0;
        self.positions.iter().for_each(|&pos| {
            if let Some(&value) = other.get(pos) {
                total = total + value;
            }
        });
        total
    }

    /// Returns the number of positions in `self`.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Returns the positions of the bits in `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![1, 2],
    ///     ],
    ///     3
    /// );
    ///
    /// assert_eq!(parity_check.check(0).unwrap().positions(), &[0, 1]);
    /// assert_eq!(parity_check.check(1).unwrap().positions(), &[1, 2]);
    /// ```
    pub fn positions(&self) -> &[usize] {
        self.positions
    }

    /// Returns the spread of `self`, that is `max(self) - min(self) + 1'. 
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
    /// assert_eq!(parity_check.check(4).unwrap().spread(), 0);
    /// ```
    pub fn spread(&self) -> usize {
        self.positions
            .last()
            .and_then(|last| 
                self.positions
                    .first()
                    .map(|first| last - first + 1)
            )
            .unwrap_or(0)
    }
}

/// An iterator over the position where a parity check matrix is 1.
/// It is ordered by row and then by col.
pub struct PositionsIter<'a> {
    active_check: usize,
    index: usize,
    check_ranges: &'a [usize],
    bit_indices: &'a [usize],
}

impl<'a> Iterator for PositionsIter<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        self.check_ranges
            .get(self.active_check + 1)
            .and_then(|&check_end| {
                if self.index >= check_end {
                    self.active_check += 1;
                }
                let position = self
                    .bit_indices
                    .get(self.index)
                    .map(|&col| (self.active_check, col));
                self.index += 1;
                position
            })
    }
}

// *************************
// Private utilitary structs
// *************************

// Helper struct to compute the rank of a parity check matrix.
struct Ranker {
    checks: Vec<Vec<Vec<usize>>>,
    n_cols: usize,
}

impl Ranker {
    fn new(matrix: &ParityCheckMatrix) -> Self {
        let mut checks = vec![Vec::new(); matrix.n_bits()];
        matrix.checks_iter().for_each(|check| {
            if let Some(col) = check.positions().first() {
                checks[*col].push(check.positions().to_vec());
            }
        });
        Self {
            checks,
            n_cols: matrix.n_bits(),
        }
    }

    fn rank(&mut self) -> usize {
        self.echelon_form();
        let mut r = 0;
        for col in 0..self.n_cols {
            if !self.checks[col].is_empty() {
                r += 1;
            }
        }
        r
    }

    // Private methodes

    fn eliminate_col(&mut self, col: usize) {
        if let Some(pivot_check) = self.checks[col].pop() {
            while let Some(other_check) = self.checks[col].pop() {
                self.insert(add_checks(&pivot_check, &other_check));
            }
            self.checks[col].push(pivot_check)
        }
    }

    fn insert(&mut self, check: Vec<usize>) {
        if let Some(&x) = check.first() {
            self.checks[x].push(check);
        }
    }

    fn echelon_form(&mut self) {
        for col in 0..self.n_cols {
            self.eliminate_col(col);
        }
    }
}

fn add_checks(check_0: &[usize], check_1: &[usize]) -> Vec<usize> {
    let mut sum = Vec::with_capacity(check_0.len() + check_1.len());
    let mut iter_0 = check_0.iter().peekable();
    let mut iter_1 = check_1.iter().peekable();
    loop {
        match (iter_0.peek(), iter_1.peek()) {
            (Some(&&v), None) => {
                sum.push(v);
                iter_0.next();
            }
            (None, Some(&&v)) => {
                sum.push(v);
                iter_1.next();
            }
            (Some(&&v_0), Some(&&v_1)) => {
                if v_0 < v_1 {
                    sum.push(v_0);
                    iter_0.next();
                } else if v_0 > v_1 {
                    sum.push(v_1);
                    iter_1.next();
                } else {
                    iter_0.next();
                    iter_1.next();
                }
            }
            (None, None) => break,
        }
    }
    sum
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn checks_iterator() {
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
        let mut iter = parity_check.checks_iter();

        assert_eq!(iter.next(), parity_check.check(0));
        assert_eq!(iter.next(), parity_check.check(1));
        assert_eq!(iter.next(), None);

        assert_eq!(parity_check.check(0).unwrap().positions(), &[0, 1]);
        assert_eq!(parity_check.check(1).unwrap().positions(), &[1, 2]);
    }

    #[test]
    fn dot_product() {
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
        let bits = vec![GF2::B0, GF2::B1, GF2::B1];

        assert_eq!(parity_check.check(0).unwrap().dot(&bits), GF2::B1);
        assert_eq!(parity_check.check(1).unwrap().dot(&bits), GF2::B0);
        assert_eq!(parity_check.dot(&bits), vec![GF2::B1, GF2::B0]);
    }

    #[test]
    fn positions_iterator() {
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
        let mut iter = parity_check.positions_iter();

        assert_eq!(iter.next(), Some((0, 0)));
        assert_eq!(iter.next(), Some((0, 1)));
        assert_eq!(iter.next(), Some((1, 1)));
        assert_eq!(iter.next(), Some((1, 2)));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn rank() {
        let parity_check_0 =
            ParityCheckMatrix::new(vec![vec![0, 1, 2, 4], vec![0, 1, 3, 5], vec![0, 2, 3, 6]], 7);
        assert_eq!(parity_check_0.rank(), 3);

        let parity_check_1 = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2], vec![0, 2]], 3);
        assert_eq!(parity_check_1.rank(), 2);
    }

    #[test]
    fn size() {
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);

        assert_eq!(parity_check.len(), 4);
        assert_eq!(parity_check.n_bits(), 3);
        assert_eq!(parity_check.n_checks(), 2);
    }
}

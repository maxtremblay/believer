//! A sparse implementation of a parity check matrix.

use crate::GF2;
use itertools::EitherOrBoth::{Both, Left, Right};
use itertools::Itertools;
use std::ops::Add;

pub mod check;
pub use check::{Check, CheckSlice};
use check::get_bitwise_sum;

pub mod check_view;
pub use check_view::CheckView;

pub mod checks_iter;
pub use checks_iter::ChecksIter;

pub mod edges_iter;
pub use edges_iter::EdgesIter;

mod ranker;
use ranker::Ranker;

mod transposer;
use transposer::Transposer;

mod concatener;
use concatener::Concatener;

/// A sparse implementation of a parity check matrix.
#[derive(Debug, PartialEq, Clone)]
pub struct ParityCheckMatrix {
    check_ranges: Vec<usize>,
    bit_indices: Vec<usize>,
    n_bits: usize,
}

impl ParityCheckMatrix {
    // ***** Construction *****

    /// Creates an empty parity check matrix. That is, a parity check with 0 bit and 0 check.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::ParityCheckMatrix;
    /// let matrix = ParityCheckMatrix::new();
    /// ```
    pub fn new() -> Self {
        Self {
            check_ranges: Vec::new(),
            bit_indices: Vec::new(),
            n_bits: 0,
        }
    }

    /// Creates a new `ParityCheckMatrix` with `n_bits` and no checks.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::ParityCheckMatrix;
    /// let matrix = ParityCheckMatrix::with_n_bits(5);
    /// ```
    pub fn with_n_bits(n_bits: usize) -> Self {
        Self {
            check_ranges: Vec::new(),
            bit_indices: Vec::new(),
            n_bits: n_bits,
        }
    }

    /// Set the checks of `self` consuming `checks`.
    ///
    /// # Panic
    ///
    /// Panics if some checks are out of bounds. That is, if they are connected to a bit that is
    /// greater or equal than `self.get_n_bits()`.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::ParityCheckMatrix;
    /// let checks = vec![vec![0, 1], vec![1, 2]];
    /// let mut matrix = ParityCheckMatrix::with_n_bits(3).with_checks(checks);
    /// ```
    pub fn with_checks(mut self, checks: Vec<Check>) -> Self {
        if self.some_checks_are_out_of_bounds(&checks) {
            panic!("some checks are out of bounds");
        }
        self.init_bit_indices(&checks);
        self.init_check_ranges(&checks);
        self.fill_with(checks);
        self
    }

    fn some_checks_are_out_of_bounds(&self, checks: &[Check]) -> bool {
        checks.iter().any(|check| self.is_out_of_bounds(check))
    }

    fn is_out_of_bounds(&self, check: CheckSlice) -> bool {
        check.iter().any(|bit| *bit >= self.n_bits)
    }

    fn init_bit_indices(&mut self, checks: &[Check]) {
        let capacity = checks.iter().map(|check| check.len()).sum();
        self.bit_indices = Vec::with_capacity(capacity)
    }

    fn init_check_ranges(&mut self, checks: &[Check]) {
        self.check_ranges = Vec::with_capacity(checks.len() + 1);
        self.check_ranges.push(0);
    }

    fn fill_with(&mut self, checks: Vec<Check>) {
        checks.into_iter().for_each(|check| if check.len() > 0 { self.add_check(check) });
    }

    fn add_check(&mut self, check: Check) {
        self.add_check_range(&check);
        self.add_bit_indices(check);
    }

    fn add_check_range(&mut self, check: CheckSlice) {
        let n_elements_before = self.bit_indices.len();
        self.check_ranges.push(n_elements_before + check.len());
    }

    fn add_bit_indices(&mut self, mut check: Check) {
        check.sort();
        self.bit_indices.append(&mut check);
    }

    /// Creates the `n_bits` identity matrix.
    /// 
    /// # Example 
    /// 
    /// ```
    /// use believer::ParityCheckMatrix;
    /// 
    /// let matrix = ParityCheckMatrix::identity_with_n_bits(3);
    /// 
    /// let identity_checks = vec![vec![0], vec![1], vec![2]];
    /// let identity_matrix = ParityCheckMatrix::with_n_bits(3)
    ///     .with_checks(identity_checks);
    /// 
    /// assert_eq!(matrix, identity_matrix);
    /// ```
    pub fn identity_with_n_bits(n_bits: usize) -> ParityCheckMatrix {
        Self {
            bit_indices: (0..n_bits).collect(),
            check_ranges: (0..n_bits + 1).collect(),
            n_bits,
        }
    }

    // ***** Getters *****

    /// Returns the number of bits in `self`.
    pub fn get_n_bits(&self) -> usize {
        self.n_bits
    }

    /// Returns the number of checks in `self`.
    pub fn get_n_checks(&self) -> usize {
        if self.check_ranges().len() > 0 {
            self.check_ranges().len() - 1
        } else {
            0
        }
    }

    /// Returns the number of edges in `self`.
    pub fn get_n_edges(&self) -> usize {
        self.bit_indices.len()
    }

    /// Returns the degree of each bit in `self`.
    ///
    /// # Example
    /// 
    /// ```
    /// use believer::ParityCheckMatrix;
    /// 
    /// let checks = vec![vec![0, 1, 2, 5], vec![1, 3, 4], vec![2, 4, 5], vec![0, 5]];
    /// let matrix = ParityCheckMatrix::with_n_bits(7).with_checks(checks);
    /// assert_eq!(matrix.get_bit_degrees(), vec![2, 2, 2, 1, 2, 3, 0]);
    /// ```
    pub fn get_bit_degrees(&self) -> Vec<usize> {
        let mut degrees = vec![0; self.n_bits];
        self.checks_iter().for_each(|check| {
            check.iter().for_each(|bit| degrees[*bit] += 1)
        });
        degrees
    }

    /// Returns the degree of each check in `self`.
    ///
    /// # Example
    /// 
    /// ```
    /// use believer::ParityCheckMatrix;
    /// 
    /// let checks = vec![vec![0, 1, 2, 5], vec![1, 3, 4], vec![2, 4, 5], vec![0, 5]];
    /// let matrix = ParityCheckMatrix::with_n_bits(7).with_checks(checks);
    /// assert_eq!(checks.get_check_degrees(), vec![4, 3, 3, 2]);
    /// ```
    pub fn get_check_degrees(&self) -> Vec<usize> {
        self.checks_iter().map(|check| check.get_n_bits()).collect()
    }

    /// Returns `Some` view over the given `check` in `self`. Returns `None` if
    /// `check` is out of bound.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::ParityCheckMatrix;
    ///
    /// let checks = vec![vec![0, 1], vec![1, 2]];
    /// 
    /// let parity_check = ParityCheckMatrix::with_n_bits(3).with_checks(checks);
    ///
    /// let check = parity_check.get_check(0).unwrap();
    /// assert_eq!(check.as_ref(), &[0, 1]);
    ///
    /// let check = parity_check.get_check(1).unwrap();
    /// assert_eq!(check.as_ref(), &[1, 2]);
    ///
    /// assert!(parity_check.get_check(2).is_none());
    /// ```
    pub fn get_check(&self, check: usize) -> Option<CheckView> {
        self.check_ranges.get(check).and_then(|&check_start| {
            self.check_ranges.get(check + 1).map(|&check_end| {
                CheckView::from_slice(&self.bit_indices[check_start..check_end])
            })
        })
    }

    /// Computes the syndrome of a given `message`.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::{GF2, ParityCheckMatrix};
    ///
    /// let checks = vec![vec![0, 1], vec![1, 2]];
    /// let parity_check = ParityCheckMatrix::with_n_bits(3).with_checks(checks);
    /// 
    /// let message = vec![GF2::B0, GF2::B1, GF2::B1];
    ///
    /// assert_eq!(parity_check.get_syndrome_of(&message), vec![GF2::B1, GF2::B0]);
    /// ```
    pub fn get_syndrome_of(&self, message: &[GF2]) -> Vec<GF2> {
        self.checks_iter().map(|check| check.compute_syndrome(message)).collect()
    }

    /// Computes the rank of `self`.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::ParityCheckMatrix;
    /// 
    /// let checks = vec![vec![0, 1], vec![1, 2], vec![0, 2]];
    /// let parity_check = ParityCheckMatrix::with_n_bits(3).with_checks(checks);
    /// 
    /// assert_eq!(parity_check.get_rank(), 2);
    /// ```
    pub fn get_rank(&self) -> usize {
        Ranker::from_parity_check_matrix(self).get_rank()
    }

    /// Gets the transposed version of `self` by swapping the bits with the checks.
    /// 
    /// # Example
    /// 
    /// ```
    /// use believer::ParityCheckMatrix;
    /// 
    /// let checks = vec![vec![0, 1, 2], vec![1, 3], vec![0, 2, 3]];
    /// let matrix = ParityCheckMatrix::with_n_bits(4).with_checks(checks);
    /// 
    /// let transposed_matrix = matrix.get_transposed_matrix();
    ///
    /// let expected_checks = vec![vec![0, 2], vec![0, 1], vec![0, 2], vec![1, 2]];
    /// let expected_matrix = ParityCheckMatrix::with_n_bits(3).with_checks(checks);
    /// 
    /// assert_eq!(transposed_matrix, expected_matrix);
    /// ```
    pub fn get_transposed_matrix(&self) -> Self {
        let checks = self.get_transposed_checks();
        Self::with_n_bits(self.get_n_checks()).with_checks(checks)
    }

    fn get_transposed_checks(&self) -> Vec<Check> {
        let mut transposer = Transposer::new();
        self.checks_iter().for_each(|check| transposer.insert_bits_from(check));
        transposer.get_checks()
    }

    // ***** Iterators *****

    /// Returns an iterator that yields a slice for each check of `self`.
    ///
    /// # Example
    /// 
    /// ```
    /// use believer::ParityCheckMatrix;
    /// 
    /// let parity_check = ParityCheckMatrix::with_n_bits(3)
    ///     .with_checks(vec![vec![0, 1], vec![1, 2]]);
    ///
    /// let mut iter = parity_check.checks_iter();
    ///
    /// assert_eq!(iter.next(), parity_check.get_check(0));
    /// assert_eq!(iter.next(), parity_check.get_check(1));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn checks_iter(&self) -> ChecksIter {
        ChecksIter::from(self)
    }

    /// An iterators over all edges in `self` ordered by check first.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let parity_check = ParityCheckMatrix::with_n_bits(3)
    ///     .with_checks(vec![vec![0, 1], vec![1, 2]]);
    /// 
    /// let mut iter = parity_check.edges_iter();
    ///
    /// assert_eq!(iter.next(), Some((0, 0)));
    /// assert_eq!(iter.next(), Some((0, 1)));
    /// assert_eq!(iter.next(), Some((1, 1)));
    /// assert_eq!(iter.next(), Some((1, 2)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn edges_iter(&self) -> EdgesIter {
        EdgesIter::from(self)
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
            .all(|check| check.compute_syndrome(message) == GF2::B0)
    }

    /// Returns a truncated parity check matrix with only the column of the given `bits`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let checks = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1, 2],
    ///         vec![2, 3, 4],
    ///         vec![0, 2, 4],
    ///         vec![1, 3],
    ///     ],
    ///     5  
    /// );
    /// let truncated_checks = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1],
    ///         vec![4],
    ///         vec![0, 4],
    ///         vec![1],
    ///     ],
    ///     5  
    /// );
    ///
    /// assert_eq!(checks.keep(&[0, 1, 4]), truncated_checks);
    /// ```
    pub fn keep(&self, bits: &[usize]) -> Self {
        let checks = self
            .checks_iter()
            .map(|check| {
                check
                    .iter()
                    .filter(|&bit| bits.iter().any(|b| b == bit))
                    .cloned()
                    .collect()
            })
            .collect();
        Self::with_n_bits(self.get_n_bits()).with_checks(checks)
    }

    /// Concatenates a parity check matrix with an other.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let x = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2], vec![2, 3]], 4);
    /// let z = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2, 3]], 4);
    ///
    /// let concat_x_z = x.right_concat(&z);
    /// let expected_concat_x_z = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1, 4, 5],
    ///         vec![1, 2, 5, 6, 7],
    ///         vec![2, 3],
    ///     ],
    ///     8
    /// );
    /// assert_eq!(concat_x_z, expected_concat_x_z);
    ///
    /// let concat_z_x = z.right_concat(&x);
    /// let expected_concat_z_x = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1, 4, 5],
    ///         vec![1, 2, 3, 5, 6],
    ///         vec![6, 7],
    ///     ],
    ///     8
    /// );
    /// assert_eq!(concat_z_x, expected_concat_z_x);
    /// ```
    pub fn right_concat(&self, other: &ParityCheckMatrix) -> ParityCheckMatrix {
        let new_checks = self
            .checks_iter()
            .zip_longest(other.checks_iter())
            .map(|x| match x {
                Both(a, b) => {
                    let mut check: Vec<usize> = a.iter().cloned().collect();
                    check.extend(b.iter().map(|p| p + self.get_n_bits()));
                    check
                }
                Left(a) => a.iter().cloned().collect(),
                Right(b) => b.iter().map(|p| p + self.get_n_bits()).collect(),
            })
            .collect();

        Self::with_n_bits(self.get_n_bits() + other.get_n_bits()).with_checks(new_checks)
    }

    pub fn diag_concat(&self, hz: &ParityCheckMatrix) -> ParityCheckMatrix {
        let new_nb_bits = self.n_bits + hz.n_bits;

        let r1 = &self.check_ranges;
        let r2 = &hz.check_ranges;

        let mut new_ranges: Vec<usize> = Vec::with_capacity(r1.len() + r2.len());
        let mut new_indices: Vec<usize> =
            Vec::with_capacity(self.bit_indices.len() + hz.bit_indices.len());

        for r in r1 {
            new_ranges.push(*r);
        }

        for r in r2 {
            //when adding the new ranges, we just need to shift these by the last range of the first matrix, since they all arrive after
            new_ranges.push(*r + r1.last().cloned().unwrap_or(0));
        }

        for i in &self.bit_indices {
            new_indices.push(*i);
        }
        for i in &hz.bit_indices {
            // similar idea here
            new_indices.push(*i + self.n_bits);
        }

        Self {
            bit_indices: new_indices,
            check_ranges: new_ranges,
            n_bits: new_nb_bits,
        }
    }

    pub fn gbc(&self, b: &ParityCheckMatrix) -> ParityCheckMatrix {
        // should check that A and B commute and that Hx*Hz^T = 0

        let hx = self.right_concat(&b);
        let hz = b.get_transposed_matrix().right_concat(&self.get_transposed_matrix());

        hx.diag_concat(&hz)
    }

    pub fn permu_matrix(l: usize) -> ParityCheckMatrix {
        let ranges: Vec<usize> = (0..l + 1).collect(); // ranges = [0,1,...,l] because we have one entries per row
        let mut indices: Vec<usize> = Vec::with_capacity(l);

        indices.push(l - 1);

        for i in 0..l - 1 {
            // indices = [l-1,0,1,...,l-2]
            indices.push(i);
        }

        Self {
            bit_indices: indices,
            check_ranges: ranges,
            n_bits: l,
        }
    }

    pub fn empty_matrix(l: usize) -> ParityCheckMatrix {
        let n_bits = l;
        let check_ranges: Vec<usize> = vec![0; l + 1];
        let bit_indices: Vec<usize> = Vec::with_capacity(0);

        Self {
            check_ranges,
            bit_indices,
            n_bits,
        }
    }

    pub fn circulant_down(indices: &Vec<usize>, l: usize) -> ParityCheckMatrix {
        let w = indices.len();
        let mut checks: Vec<Vec<usize>> = Vec::with_capacity(l);

        for i in 0..l {
            let mut new_row: Vec<usize> = Vec::with_capacity(w);
            for j in indices {
                new_row.push((l - j + i) % l);
            }
            checks.push(new_row);
        }

        ParityCheckMatrix::with_n_bits(l).with_checks(checks)
    }

    pub fn circulant_right(indices: &Vec<usize>, l: usize) -> ParityCheckMatrix {
        let w = indices.len();
        let mut checks: Vec<Vec<usize>> = Vec::with_capacity(l);

        for i in 0..l {
            let mut new_row: Vec<usize> = Vec::with_capacity(w);
            for j in indices {
                new_row.push((j + i) % l);
            }
            checks.push(new_row);
        }

        ParityCheckMatrix::with_n_bits(l).with_checks(checks)
    }

    /// Returns a truncated parity check matrix where the column of the given `bits` are remove.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let checks = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1, 2],
    ///         vec![2, 3, 4],
    ///         vec![0, 2, 4],
    ///         vec![1, 3],
    ///     ],
    ///     5  
    /// );
    /// let truncated_checks = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![1],
    ///         vec![3, 4],
    ///         vec![4],
    ///         vec![1, 3],
    ///     ],
    ///     5  
    /// );
    ///
    /// assert_eq!(checks.without(&[0, 2]), truncated_checks);
    /// ```
    pub fn without(&self, bits: &[usize]) -> Self {
        let to_keep: Vec<usize> = (0..9).filter(|x| !bits.contains(x)).collect();
        self.keep(&to_keep)
    }

    // Returns a reference to `self.check_ranges`.
    pub(crate) fn check_ranges(&self) -> &[usize] {
        &self.check_ranges
    }

    // Returns a reference to `self.bit_indices
    pub(crate) fn bit_indices(&self) -> &[usize] {
        &self.bit_indices
    }
}

impl std::fmt::Display for ParityCheckMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for check in self.checks_iter() {
            write!(f, "[ ")?;
            for bit in check.iter() {
                write!(f, "{} ", bit)?;
            }
            write!(f, "]")?;
        }
        Ok(())
    }
}

// ************************
// Public utilitary structs
// ************************



impl Add for &ParityCheckMatrix {
    type Output = ParityCheckMatrix;

    fn add(self, b: &ParityCheckMatrix) -> ParityCheckMatrix {
        // should check same dimensions

        let nb_bits = self.n_bits;
        let nb_rows = self.check_ranges.len() - 1;
        let mut new_rows: Vec<Vec<usize>> = Vec::with_capacity(nb_rows);

        for nb_row in 0..nb_rows {
            let a_row: Vec<usize> = self.get_check(nb_row).unwrap().iter().cloned().collect();
            let b_row: Vec<usize> = b.get_check(nb_row).unwrap().iter().cloned().collect();

            let mut accum: Vec<usize> = Vec::with_capacity(a_row.len() + b_row.len()); // in the worst case the 1s are in different places and don't cancel out

            let mut i = 0;
            let mut j = 0;

            while (i < a_row.len()) && (j < b_row.len()) {
                if a_row[i] == b_row[j] {
                    // 1+1, we do nothing
                    i += 1;
                    j += 1;
                } else if a_row[i] > b_row[j] {
                    accum.push(b_row[j]);
                    j += 1;
                } else {
                    accum.push(a_row[i]);
                    i += 1;
                }
            }

            // because a and b can have different length we take care of adding potential leftovers
            if j >= b_row.len() {
                for k in i..a_row.len() {
                    accum.push(a_row[k]);
                }
            } else if i >= a_row.len() {
                for k in j..b_row.len() {
                    accum.push(b_row[k]);
                }
            }

            new_rows.push(accum);
        }

        ParityCheckMatrix::with_n_bits(nb_bits).with_checks(new_rows)
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn checks_are_sorted_on_construction() {
        let checks = vec![vec![1, 0], vec![0, 2, 1], vec![1, 2, 3]];
        let matrix = ParityCheckMatrix::with_n_bits(4).with_checks(checks);

        assert_eq!(matrix.get_check(0).unwrap().as_ref(), &[0, 1]);
        assert_eq!(matrix.get_check(1).unwrap().as_ref(), &[0, 1, 2]);
        assert_eq!(matrix.get_check(2).unwrap().as_ref(), &[1, 2, 3]);
    }

    #[test]
    fn empty_checks_are_removed_on_construction() {
        let checks = vec![vec![], vec![0, 1], vec![], vec![1, 2]];
        let matrix = ParityCheckMatrix::with_n_bits(3).with_checks(checks);

        assert_eq!(matrix.get_check(0).unwrap().as_ref(), &[0, 1]);
        assert_eq!(matrix.get_check(1).unwrap().as_ref(), &[1, 2]);
        assert_eq!(matrix.get_n_checks(), 2);
    }

    #[test]
    #[should_panic]
    fn panics_on_construction_if_checks_are_out_of_bound() {
        let checks = vec![vec![0, 1, 5], vec![2, 3, 4]];
        ParityCheckMatrix::with_n_bits(5).with_checks(checks);
    }

    #[test]
    fn syndrome() {
        let parity_check = ParityCheckMatrix::with_n_bits(3)
            .with_checks(vec![vec![0, 1], vec![1, 2]]);
        let bits = vec![GF2::B0, GF2::B1, GF2::B1];

        assert_eq!(
            parity_check.get_check(0).unwrap().compute_syndrome(&bits),
            GF2::B1
        );
        assert_eq!(
            parity_check.get_check(1).unwrap().compute_syndrome(&bits),
            GF2::B0
        );
        assert_eq!(parity_check.get_syndrome_of(&bits), vec![GF2::B1, GF2::B0]);
    }
}

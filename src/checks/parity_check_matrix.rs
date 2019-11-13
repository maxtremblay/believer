//! A sparse implementation of a parity check matrix.

use crate::GF2;
use itertools::EitherOrBoth::{Both, Left, Right};
use itertools::Itertools;
use std::cmp::Ordering;
use std::cmp::{max, min};
use std::ops::{Add, Mul};
use super::Check;

// *******************
// Parity Check Matrix
// *******************

/// A sparse parity check matrix.
#[derive(Debug, PartialEq, Clone)]
pub struct ParityCheckMatrix {
    check_ranges: Vec<usize>,
    bit_indices: Vec<usize>,
    n_bits: usize,
}

impl ParityCheckMatrix {
    // ************
    // Construction
    // ************

    /// Creates a new `ParityCheckMatrix` with `n_bits` and no checks.
    /// 
    /// # Example 
    /// 
    /// ```
    /// # use believer::*;
    /// let matrix = ParityCheckMatrix::with_n_bits(5);
    /// assert_eq!(matrix.n_bits(), 5);
    /// assert_eq!(matrix.n_checks(), 0);
    /// ```
    pub fn with_n_bits(n_bits: usize) -> Self {
        Self {
            check_ranges: vec![],
            bit_indices: vec![],
            n_bits
        }
    }

    /// Set the checks of `self` consuming `checks`. 
    /// 
    /// # Panic 
    /// 
    /// Panics if some checks are out of bounds. That is, if they are connected to a bit that is
    /// greater or equal than `self.n_bits()`.
    /// 
    /// # Example 
    /// 
    /// ```
    /// # use believer::*;
    /// let checks = vec![vec![0, 1], vec![1, 2]];
    /// let mut matrix = ParityCheckMatrix::with_n_bits(3).with_checks(checks);
    /// assert_eq!(matrix.n_checks(), 2);
    /// ```
    pub fn with_checks(mut self, checks: Vec<Vec<usize>>) -> Self {
        if self.some_checks_are_out_of_bounds(&checks) {
            panic!("some checks are out of bounds");
        }
        self.init_bit_indices(&checks);
        self.init_check_ranges(&checks);
        self.fill_with(checks);
        self
    }

    fn some_checks_are_out_of_bounds(&self, checks: &[Vec<usize>]) -> bool {
        checks.iter().any(|check| self.is_out_of_bounds(check))
    }

    fn is_out_of_bounds(&self, check: &[usize]) -> bool {
        check.iter().max().map(|max| *max >= self.n_bits).unwrap_or(false)
    }

    fn init_bit_indices(&mut self, checks: &[Vec<usize>]) {
        let capacity = checks.iter().fold(0, |acc, check| acc + check.len());
        self.bit_indices = Vec::with_capacity(capacity)
    }

    fn init_check_ranges(&mut self, checks: &[Vec<usize>]) {
        self.check_ranges = Vec::with_capacity(checks.len() + 1);
        self.check_ranges.push(0);
    }

    fn fill_with(&mut self, checks: Vec<Vec<usize>>) {
        checks.into_iter().for_each(|check| {
            if check.len() > 0 {
                self.add_check(check);
            }
        })
    }

    fn add_check(&mut self, check: Vec<usize>) {
        self.add_check_range(&check);
        self.add_bit_indices(check);
        
    }

    fn add_check_range(&mut self, check: &[usize]) {
        let n_elements_before = self.bit_indices.len();
        self.check_ranges.push(n_elements_before + check.len());
    }

    fn add_bit_indices(&mut self, mut check: Vec<usize>) {
        check.sort();
        self.bit_indices.append(&mut check);
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
    /// let checks = vec![vec![0, 1], vec![1, 2]];
    /// let n_bits = 3;
    /// let parity_check = ParityCheckMatrix::new(checks, n_bits);
    /// ```
    pub fn new(checks: Vec<Vec<usize>>, n_bits: usize) -> Self {
        Self::with_n_bits(n_bits).with_checks(checks)
    }

    // *******
    // Getters
    // *******

    /// Returns the degree of each bit in `self`.
    ///
    /// # Example
    /// ```
    /// # use believer::*;
    /// let checks = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1, 2, 5],
    ///         vec![1, 3, 4],
    ///         vec![2, 4, 5],
    ///         vec![0, 5],
    ///     ],
    ///     7,
    /// );
    /// assert_eq!(checks.get_bit_degrees(), vec![2, 2, 2, 1, 2, 3, 0]);
    /// ```
    pub fn get_bit_degrees(&self) -> Vec<usize> {
        let mut degrees = vec![0; self.n_bits];
        self.checks_iter().for_each(|check| {
            check.iter().for_each(|bit| {
                degrees[*bit] += 1;
            })
        });
        degrees
    }

    /// Returns the degree of each check in `self`.
    ///
    /// # Example
    /// ```
    /// # use believer::*;
    /// let checks = ParityCheckMatrix::new(
    ///     vec![
    ///         vec![0, 1, 2, 5],
    ///         vec![1, 3, 4],
    ///         vec![2, 4, 5],
    ///         vec![0, 5],
    ///     ],
    ///     7,
    /// );
    /// assert_eq!(checks.get_check_degrees(), vec![4, 3, 3, 2]);
    /// ```
    pub fn get_check_degrees(&self) -> Vec<usize> {
        self.checks_iter().map(|check| check.get_n_bits()).collect()
    }

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
    /// assert_eq!(slice.compute_syndrome(&vector), GF2::B0);
    /// ```
    pub fn check(&self, check: usize) -> Option<Check> {
        self.check_ranges.get(check).and_then(|&check_start| {
            self.check_ranges.get(check + 1).map(|&check_end| {
                Check::from_slice(&self.bit_indices[check_start..check_end])
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
        self.checks_iter().map(|check| check.compute_syndrome(vector)).collect()
    }

    /// Creates an empty parity check matrix. That is, a parity check with 0 bit and 0 check.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let empty_matrix = ParityCheckMatrix::empty();
    /// assert_eq!(empty_matrix.n_bits(), 0);
    /// assert_eq!(empty_matrix.n_checks(), 0);
    /// ```
    pub fn empty() -> Self {
        Self::new(Vec::new(), 0)
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

    /// Returns `true` is there is no value in `self`.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of 1 in `self`.
    pub fn len(&self) -> usize {
        self.bit_indices.len()
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
        Self::new(checks, self.n_bits())
    }

    /// Returns the number of bits in `self`.
    pub fn n_bits(&self) -> usize {
        self.n_bits
    }

    /// Returns the number of checks in `self`.
    pub fn n_checks(&self) -> usize {
        if self.check_ranges().len() > 0 {
            self.check_ranges().len() - 1
        } else {
            0
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

    pub fn transpose(&self) -> Self {
        let mut indices = Vec::with_capacity(self.n_bits());
        let mut column_indices = Vec::with_capacity(self.n_bits());
        let mut row_ranges = Vec::new();
        row_ranges.push(0);

        let mut active_col = 0;
        let mut row_lenght = 0;

        self.positions_iter()
            .enumerate()
            .sorted_by(|(_, (r_0, c_0)), (_, (r_1, c_1))| match c_0.cmp(c_1) {
                Ordering::Equal => r_0.cmp(r_1),
                otherwise => otherwise,
            })
            .for_each(|(idx, (row, col))| {
                if col == active_col {
                    row_lenght += 1;
                } else {
                    while active_col < col {
                        active_col += 1;
                        row_ranges.push(*row_ranges.last().unwrap_or(&0) + row_lenght);
                    }
                    row_lenght = 1;
                }
                column_indices.push(row);
                indices.push(idx);
            });

        row_ranges.push(*row_ranges.last().unwrap_or(&0) + row_lenght);

        let n_bits = self.check_ranges().len() - 1;

        Self {
            check_ranges: row_ranges,
            bit_indices: column_indices,
            n_bits,
        }
    }

    pub fn sparse_dot(v1: &[usize], v2: &[usize]) -> GF2 {
        let mut i = 0;
        let mut j = 0;
        let mut accum = GF2::B0;

        while (i < v1.len()) && (j < v2.len()) {
            if v1[i] == v2[j] {
                accum = accum + GF2::B1;
                i += 1;
                j += 1;
            } else if v1[i] > v2[j] {
                j += 1;
            } else {
                i += 1;
            }
        }

        accum
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
                    check.extend(b.iter().map(|p| p + self.n_bits()));
                    check
                }
                Left(a) => a.iter().cloned().collect(),
                Right(b) => b.iter().map(|p| p + self.n_bits()).collect(),
            })
            .collect();

        Self::new(new_checks, self.n_bits() + other.n_bits())
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
        let hz = b.transpose().right_concat(&self.transpose());

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

    pub fn ident(l: usize) -> ParityCheckMatrix {
        //returns the identity matrix of dimension l
        let n_bits = l;
        let check_ranges: Vec<usize> = (0..l + 1).collect();
        let bit_indices: Vec<usize> = (0..l).collect();

        Self {
            bit_indices,
            check_ranges,
            n_bits,
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

        ParityCheckMatrix::new(checks, l)
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

        ParityCheckMatrix::new(checks, l)
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

impl Add for &ParityCheckMatrix {
    type Output = ParityCheckMatrix;

    fn add(self, b: &ParityCheckMatrix) -> ParityCheckMatrix {
        // should check same dimensions

        let nb_bits = self.n_bits;
        let nb_rows = self.check_ranges.len() - 1;
        let mut new_rows: Vec<Vec<usize>> = Vec::with_capacity(nb_rows);

        for nb_row in 0..nb_rows {
            let a_row: Vec<usize> = self.check(nb_row).unwrap().iter().cloned().collect();
            let b_row: Vec<usize> = b.check(nb_row).unwrap().iter().cloned().collect();

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

        ParityCheckMatrix::new(new_rows, nb_bits)
    }
}

// NOTE: IS IT USEFUL??? 

// impl Mul for &ParityCheckMatrix {
//     type Output = ParityCheckMatrix;

//     fn mul(self, b: &ParityCheckMatrix) -> ParityCheckMatrix {
//         // should check for dimensions, needs some optimization

//         //let mut new_row = vec![GF2::B0; B.n_bits]; // initializes vector that will temporarily store a computed row of the end matrix
//         let nb_cols = b.n_bits;
//         let nb_rows = self.check_ranges.len() - 1;
//         let mut accum_new_rows: Vec<Vec<usize>> = Vec::with_capacity(nb_rows);

//         for v1 in self.checks_iter() {
//             let mut sparse_row = Vec::with_capacity(nb_cols); // in the worst case we have entries for every column

//             for (i, v2) in b.transpose().checks_iter().enumerate() {
//                 let dot_res = ParityCheckMatrix::sparse_dot(
//                     &v1.iter().cloned().collect(), 
//                     &v2.iter().cloned().collect()
//                 );

//                 if let GF2::B1 = dot_res {
//                     //only if the dot product returns a non null value do we mention the column in the sparse row
//                     sparse_row.push(i);
//                 }
//             }

//             accum_new_rows.push(sparse_row);
//         }

//         ParityCheckMatrix::new(accum_new_rows, nb_cols)
//     }
// }

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
            checks[check.min()].push(check.iter().cloned().collect());
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
    fn checks_are_sorted_on_construction() {
        let checks = vec![vec![1, 0], vec![0, 2, 1], vec![1, 2, 3]];
        let matrix = ParityCheckMatrix::with_n_bits(4).with_checks(checks);

        assert_eq!(matrix.check(0).unwrap().as_ref(), &[0, 1]);
        assert_eq!(matrix.check(1).unwrap().as_ref(), &[0, 1, 2]);
        assert_eq!(matrix.check(2).unwrap().as_ref(), &[1, 2, 3]);
    }

    #[test]
    fn empty_checks_are_removed_on_construction() {
        let checks = vec![vec![], vec![0, 1], vec![], vec![1, 2]];
        let matrix = ParityCheckMatrix::with_n_bits(3).with_checks(checks);

        assert_eq!(matrix.check(0).unwrap().as_ref(), &[0, 1]);
        assert_eq!(matrix.check(1).unwrap().as_ref(), &[1, 2]);
        assert_eq!(matrix.n_checks(), 2);
    }

    #[test]
    #[should_panic]
    fn panics_on_construction_if_checks_are_out_of_bound() {
        let checks = vec![vec![0, 1, 5], vec![2, 3, 4]];
        ParityCheckMatrix::with_n_bits(5).with_checks(checks);
    }

    #[test]
    fn checks_iterator() {
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
        let mut iter = parity_check.checks_iter();

        assert_eq!(iter.next(), parity_check.check(0));
        assert_eq!(iter.next(), parity_check.check(1));
        assert_eq!(iter.next(), None);

        assert_eq!(parity_check.check(0).unwrap().as_ref(), &[0, 1]);
        assert_eq!(parity_check.check(1).unwrap().as_ref(), &[1, 2]);
    }

    #[test]
    fn dot_product() {
        let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]], 3);
        let bits = vec![GF2::B0, GF2::B1, GF2::B1];

        assert_eq!(parity_check.check(0).unwrap().compute_syndrome(&bits), GF2::B1);
        assert_eq!(parity_check.check(1).unwrap().compute_syndrome(&bits), GF2::B0);
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
        let parity_check_0 = ParityCheckMatrix::new(
            vec![vec![0, 1, 2, 4], vec![0, 1, 3, 5], vec![0, 2, 3, 6]],
            7,
        );
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

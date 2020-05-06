//! A sparse implementation of a parity check matrix.

pub mod check;
pub use check::{Check, CheckSlice};

pub mod checks_iter;
pub use checks_iter::Checks;

pub mod edges_iter;
pub use edges_iter::Edges;

mod ranker;
use ranker::Ranker;

mod transposer;
use transposer::Transposer;

mod concatener;
use concatener::Concatener;

use qecsim::code::Code;
use qecsim::decoding::syndrome_decoding::{CheckState, Syndrome};
use qecsim::noise::binary::BinaryError;
use qecsim::noise::ErrorBlock;

/// A sparse implementation of a parity check matrix.
#[derive(Debug, PartialEq, Clone)]
pub struct ParityCheckMatrix {
    check_ranges: Vec<usize>,
    bit_indices: Vec<usize>,
    block_size: usize,
}

impl Code for ParityCheckMatrix {
    type Error = BinaryError;

    fn block_size(&self) -> usize {
        self.block_size
    }

    fn dimension(&self) -> usize {
        self.block_size - self.rank()
    }

    fn erased_at<P>(&self, positions: P) -> Self
    where
        P: IntoIterator<Item = usize>,
    {
        let mut checks: Vec<Check> = self.checks().map(|check| check.to_vec()).collect();
        for bit in positions.into_iter() {
            for check in checks.iter_mut() {
                if let Some(position) = check.iter().position(|b| *b == bit).clone() {
                    check.swap_remove(position);
                }
            }
        }
        ParityCheckMatrix::with_block_size(self.block_size).with_checks(checks)
    }

    fn is_equivalent_to_identity(&self, error_block: &ErrorBlock<BinaryError>) -> bool {
        error_block.weight() == 0
    }

    fn syndrome_of(&self, error_block: &ErrorBlock<BinaryError>) -> Syndrome {
        self.checks()
            .map(|check| {
                if error_block.satisfy_check(check.to_vec()) {
                    CheckState::Satisfied
                } else {
                    CheckState::Unsatisfied
                }
            })
            .collect()
    }
}

impl ParityCheckMatrix {
    /// Creates an empty parity check matrix. That is, a parity check with 0 bit and 0 check.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    /// let matrix = ParityCheckMatrix::new();
    /// ```
    pub fn new() -> Self {
        Self {
            check_ranges: Vec::new(),
            bit_indices: Vec::new(),
            block_size: 0,
        }
    }

    /// Creates a new `ParityCheckMatrix` with `block_size` and no checks.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    /// let matrix = ParityCheckMatrix::with_block_size(5);
    /// ```
    pub fn with_block_size(block_size: usize) -> Self {
        Self {
            check_ranges: Vec::new(),
            bit_indices: Vec::new(),
            block_size,
        }
    }

    /// Set the checks of `self` consuming `checks`.
    ///
    /// It also sort the checks before inserting them.
    ///
    /// # Panic
    ///
    /// Panics if some checks are out of bounds. That is, if they are connected to a bit that is
    /// greater or equal than `self.block_size()`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    /// let checks = vec![vec![0, 1], vec![1, 2]];
    /// let mut matrix = ParityCheckMatrix::with_block_size(3).with_checks(checks);
    /// ```
    pub fn with_checks(mut self, checks: Vec<Check>) -> Self {
        if !checks.is_empty() {
            if self.some_checks_are_out_of_bounds(&checks) {
                panic!("some checks are out of bounds");
            }
            self.init_bit_indices(&checks);
            self.init_check_ranges(&checks);
            self.fill_with(checks);
        }
        self
    }

    fn some_checks_are_out_of_bounds(&self, checks: &[Check]) -> bool {
        checks.iter().any(|check| self.is_out_of_bounds(check))
    }

    fn is_out_of_bounds(&self, check: CheckSlice) -> bool {
        check.iter().any(|bit| *bit >= self.block_size)
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
        checks
            .into_iter()
            .filter(|check| check.len() > 0)
            .for_each(|check| self.add_check(check));
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

    /// Creates an identity matrix of the given block size.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    ///
    /// let matrix = ParityCheckMatrix::identity(3);
    ///
    /// let identity_checks = vec![vec![0], vec![1], vec![2]];
    /// let identity_matrix = ParityCheckMatrix::with_block_size(3)
    ///     .with_checks(identity_checks);
    ///
    /// assert_eq!(matrix, identity_matrix);
    /// ```
    pub fn identity(block_size: usize) -> ParityCheckMatrix {
        Self {
            bit_indices: (0..block_size).collect(),
            check_ranges: (0..block_size + 1).collect(),
            block_size,
        }
    }

    /// Returns the number of bits.
    pub fn number_of_bits(&self) -> usize {
        self.block_size
    }

    /// Returns the number of checks.
    pub fn number_of_checks(&self) -> usize {
        match self.check_ranges.len() {
            0 => 0,
            n => n - 1,
        }
    }

    /// Returns the number of edges.
    pub fn number_of_edges(&self) -> usize {
        self.bit_indices.len()
    }

    /// Returns the degree of each bit in `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    ///
    /// let checks = vec![vec![0, 1, 2, 5], vec![1, 3, 4], vec![2, 4, 5], vec![0, 5]];
    /// let matrix = ParityCheckMatrix::with_block_size(7).with_checks(checks);
    /// assert_eq!(matrix.bit_degrees().collect::<Vec<usize>>(), vec![2, 2, 2, 1, 2, 3, 0]);
    /// ```
    pub fn bit_degrees<'a>(&'a self) -> impl Iterator<Item = usize> + 'a {
        (0..self.block_size)
            .map(move |bit| self.checks().filter(|check| check.contains(&bit)).count())
    }

    /// Returns the degree of each check in `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    ///
    /// let checks = vec![vec![0, 1, 2, 5], vec![1, 3, 4], vec![2, 4, 5], vec![0, 5]];
    /// let matrix = ParityCheckMatrix::with_block_size(7).with_checks(checks);
    /// assert_eq!(matrix.check_degrees().collect::<Vec<usize>>(), vec![4, 3, 3, 2]);
    /// ```
    pub fn check_degrees<'a>(&'a self) -> impl Iterator<Item = usize> + 'a {
        self.checks().map(|check| check.len())
    }

    /// Returns `Some` view over the given `check` in `self`. Returns `None` if
    /// `check` is out of bound.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    ///
    /// let checks = vec![vec![0, 1], vec![1, 2]];
    ///
    /// let parity_check = ParityCheckMatrix::with_block_size(3).with_checks(checks);
    ///
    /// let check = parity_check.check(0).unwrap();
    /// assert_eq!(check.as_ref(), &[0, 1]);
    ///
    /// let check = parity_check.check(1).unwrap();
    /// assert_eq!(check.as_ref(), &[1, 2]);
    ///
    /// assert!(parity_check.check(2).is_none());
    /// ```
    pub fn check(&self, check: usize) -> Option<CheckSlice> {
        let check_start = self.check_ranges.get(check)?;
        let check_end = self.check_ranges.get(check + 1)?;
        Some(&self.bit_indices[*check_start..*check_end])
    }

    /// Computes the rank of `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    ///
    /// let checks = vec![vec![0, 1], vec![1, 2], vec![0, 2]];
    /// let parity_check = ParityCheckMatrix::with_block_size(3).with_checks(checks);
    ///
    /// assert_eq!(parity_check.rank(), 2);
    /// ```
    pub fn rank(&self) -> usize {
        Ranker::from_parity_check_matrix(self).rank()
    }

    /// Gets the transposed version of `self` by swapping the bits with the checks.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    ///
    /// let checks = vec![vec![0, 1, 2], vec![1, 3], vec![0, 2, 3]];
    /// let matrix = ParityCheckMatrix::with_block_size(4).with_checks(checks);
    ///
    /// let transposed_matrix = matrix.transposed();
    ///
    /// let expected_checks = vec![vec![0, 2], vec![0, 1], vec![0, 2], vec![1, 2]];
    /// let expected_matrix = ParityCheckMatrix::with_block_size(3).with_checks(expected_checks);
    ///
    /// assert_eq!(transposed_matrix, expected_matrix);
    /// ```
    pub fn transposed(&self) -> Self {
        Transposer::from(self).transposed()
    }

    /// Returns the horizontal concatenation of `self` with `other`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    /// let left_matrix = ParityCheckMatrix::with_block_size(3)
    ///     .with_checks(vec![vec![0, 1], vec![1, 2]]);
    /// let right_matrix = ParityCheckMatrix::with_block_size(4)
    ///     .with_checks(vec![vec![1, 2, 3], vec![0, 1], vec![2, 3]]);
    ///
    /// let concatened = left_matrix.horizontal_concat_with(&right_matrix);
    ///
    /// let expected = ParityCheckMatrix::with_block_size(7)
    ///     .with_checks(vec![vec![0, 1, 4, 5, 6], vec![1, 2, 3, 4], vec![5, 6]]);
    ///
    /// assert_eq!(concatened, expected);
    /// ```
    pub fn horizontal_concat_with(&self, other: &ParityCheckMatrix) -> ParityCheckMatrix {
        Concatener::from(self, other).concat_horizontally()
    }

    /// Returns the diagonal concatenation of `self` with `other`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    /// let left_matrix = ParityCheckMatrix::with_block_size(3)
    ///     .with_checks(vec![vec![0, 1], vec![1, 2]]);
    /// let right_matrix = ParityCheckMatrix::with_block_size(4)
    ///     .with_checks(vec![vec![1, 2, 3], vec![0, 1], vec![2, 3]]);
    ///
    /// let concatened = left_matrix.diagonal_concat_with(&right_matrix);
    ///
    /// let expected = ParityCheckMatrix::with_block_size(7)
    ///     .with_checks(vec![vec![0, 1], vec![1, 2], vec![4, 5, 6], vec![3, 4], vec![5, 6]]);
    ///
    /// assert_eq!(concatened, expected);
    /// ```
    pub fn diagonal_concat_with(&self, other: &ParityCheckMatrix) -> ParityCheckMatrix {
        Concatener::from(self, other).concat_diagonally()
    }

    /// Returns an iterator that yields a slice for each check of `self`.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    ///
    /// let parity_check = ParityCheckMatrix::with_block_size(3)
    ///     .with_checks(vec![vec![0, 1], vec![1, 2]]);
    ///
    /// let mut iter = parity_check.checks();
    ///
    /// assert_eq!(iter.next(), parity_check.check(0));
    /// assert_eq!(iter.next(), parity_check.check(1));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn checks(&self) -> Checks {
        Checks::from(self)
    }

    /// An iterators over all edges in `self` ordered by check first.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::ParityCheckMatrix;
    /// let parity_check = ParityCheckMatrix::with_block_size(3)
    ///     .with_checks(vec![vec![0, 1], vec![1, 2]]);
    ///
    /// let mut iter = parity_check.edges();
    ///
    /// assert_eq!(iter.next(), Some((0, 0)));
    /// assert_eq!(iter.next(), Some((0, 1)));
    /// assert_eq!(iter.next(), Some((1, 1)));
    /// assert_eq!(iter.next(), Some((1, 2)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn edges(&self) -> Edges {
        Edges::from(self)
    }
    /*
    pub fn gbc(&self, b: &ParityCheckMatrix) -> ParityCheckMatrix {
        // should check that A and B commute and that Hx*Hz^T = 0
        let hx = self.horizontal_concat_with(&b);
        let hz = b
            .transposed_matrix()
            .horizontal_concat_with(&self.transposed_matrix());
        hx.diagonal_concat_with(&hz)
    }

    pub fn gbc_from_poly(a: &[usize], b: &[usize], l: usize) -> ParityCheckMatrix {
        //l is block_size_a = block_size_b, and the number of physical qubits is 2*l
        let w_a = a.len();
        let w_b = b.len();
        let w = w_a + w_b;
        let check_ranges: Vec<usize> = (0..2 * l * w + 1).step_by(w).collect();
        let mut bit_indices: Vec<usize> = Vec::with_capacity(2 * l * w);

        for i in 0..l {
            for j in a {
                bit_indices.push((j + i) % l);
            }
            (&mut bit_indices[i * w..i * w + w_a]).sort_unstable();
            for j in b {
                bit_indices.push((j + i) % l + l);
            }
            (&mut bit_indices[i * w + w_a..(i + 1) * w]).sort_unstable();
        }
        for i in l..2 * l {
            for j in b {
                bit_indices.push((l - j + i) % l + 2 * l);
            }
            (&mut bit_indices[i * w..i * w + w_b]).sort_unstable();
            for j in a {
                bit_indices.push((l - j + i) % l + 3 * l);
            }
            (&mut bit_indices[i * w + w_b..(i + 1) * w]).sort_unstable();
        }

        Self {
            block_size: 4 * l,
            check_ranges,
            bit_indices,
        }
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
            block_size: l,
        }
    }

    pub fn circulant_down(indices: &[usize], l: usize) -> ParityCheckMatrix {
        let w = indices.len();
        let mut checks: Vec<Vec<usize>> = Vec::with_capacity(l);

        for i in 0..l {
            let mut new_row: Vec<usize> = Vec::with_capacity(w);
            for j in indices {
                new_row.push((l - j + i) % l);
            }
            checks.push(new_row);
        }

        ParityCheckMatrix::with_block_size(l).with_checks(checks)
    }

    pub fn circulant_right(indices: &[usize], l: usize) -> ParityCheckMatrix {
        let w = indices.len();
        let mut checks: Vec<Vec<usize>> = Vec::with_capacity(l);

        for i in 0..l {
            let mut new_row: Vec<usize> = Vec::with_capacity(w);
            for j in indices {
                new_row.push((j + i) % l);
            }
            checks.push(new_row);
        }

        ParityCheckMatrix::with_block_size(l).with_checks(checks)
    }

    // Returns a reference to `self.check_ranges`.
    pub(crate) fn check_ranges(&self) -> &[usize] {
        &self.check_ranges
    }

    // Returns a reference to `self.bit_indices
    #[allow(dead_code)]
    pub(crate) fn bit_indices(&self) -> &[usize] {
        &self.bit_indices
    }
    */
}
/*
fn transfer_to(v1: &[usize], v2: &mut Vec<usize>) {
    unsafe { v2.set_len(v1.len()) }
    for i in 0..v1.len() {
        v2[i] = v1[i];
    }
}

pub fn add_checks_mut(check_0: &[usize], check_1: &[usize], sum: &mut Vec<usize>) {
    //let mut sum = Vec::with_capacity(check_0.len() + check_1.len());
    sum.clear();
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
}

pub fn binary_search(data: &[usize], target: &usize) -> bool {
    // println!("data:{:?}",data);
    // println!("target:{:?}",target);
    if data.len() > 0 {
        let mut high = data.len() - 1;
        let mut low = 0;

        if target < &data[0] || target > &data[high] {
            return false;
        }

        while low <= high {
            let mid = (low + high) / 2;

            if data[mid] < *target {
                low = mid + 1
            } else if data[mid] > *target {
                high = mid - 1
            } else if data[mid] == *target {
                return true;
            }
        }
        return false;
    } else {
        return false;
    }
}

pub fn add_checks(check_0: &[usize], check_1: &[usize]) -> Vec<usize> {
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
*/
impl std::fmt::Display for ParityCheckMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for check in self.checks() {
            write!(f, "[ ")?;
            for bit in check.iter() {
                write!(f, "{} ", bit)?;
            }
            write!(f, "]")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn checks_are_sorted_on_construction() {
        let checks = vec![vec![1, 0], vec![0, 2, 1], vec![1, 2, 3]];
        let matrix = ParityCheckMatrix::with_block_size(4).with_checks(checks);

        assert_eq!(matrix.check(0).unwrap().as_ref(), &[0, 1]);
        assert_eq!(matrix.check(1).unwrap().as_ref(), &[0, 1, 2]);
        assert_eq!(matrix.check(2).unwrap().as_ref(), &[1, 2, 3]);
    }

    #[test]
    fn empty_checks_are_removed_on_construction() {
        let checks = vec![vec![], vec![0, 1], vec![], vec![1, 2]];
        let matrix = ParityCheckMatrix::with_block_size(3).with_checks(checks);

        assert_eq!(matrix.check(0).unwrap().as_ref(), &[0, 1]);
        assert_eq!(matrix.check(1).unwrap().as_ref(), &[1, 2]);
        assert_eq!(matrix.number_of_checks(), 2);
    }

    #[test]
    #[should_panic]
    fn panics_on_construction_if_checks_are_out_of_bound() {
        let checks = vec![vec![0, 1, 5], vec![2, 3, 4]];
        ParityCheckMatrix::with_block_size(5).with_checks(checks);
    }

    #[test]
    fn syndrome() {
        let parity_check =
            ParityCheckMatrix::with_block_size(3).with_checks(vec![vec![0, 1], vec![1, 2]]);
        let mut error_block = ErrorBlock::new();
        error_block.insert_at(0, BinaryError::Flipped);
        let syndrome = vec![CheckState::Unsatisfied, CheckState::Satisfied].into();
        assert_eq!(parity_check.syndrome_of(&error_block), syndrome);
    }

    #[test]
    fn check_if_error_is_equivalent_to_identity() {
        let parity_check =
            ParityCheckMatrix::with_block_size(3).with_checks(vec![vec![0, 1], vec![1, 2]]);
        let no_error = ErrorBlock::new();
        assert!(parity_check.is_equivalent_to_identity(&no_error));
        let mut error_block = ErrorBlock::new();
        error_block.insert_at(0, BinaryError::Flipped);
        assert!(parity_check.is_not_equivalent_to_identity(&error_block));
    }
}

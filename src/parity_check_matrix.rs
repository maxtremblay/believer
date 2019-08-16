use crate::GF2;

#[derive(Debug, PartialEq)]
pub struct ParityCheckMatrix {
    row_ranges: Vec<usize>,
    column_indices: Vec<usize>,
}

impl ParityCheckMatrix {
    /// Computes the dot product between `self` and a binary vector.
    /// 
    /// # Example 
    /// 
    /// ```
    /// # use::believer::*;
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     (0, 0),
    ///     (0, 1),    
    ///     (1, 1),
    ///     (1, 2),
    /// ]);
    /// let vector = vec![GF2::B0, GF2::B1, GF2::B1];
    /// 
    /// assert_eq!(parity_check.dot(&vector), vec![GF2::B1, GF2::B0]);
    /// ```
    pub fn dot(&self, vector: &[GF2]) -> Vec<GF2> {
        self.rows_iter()
            .map(|row| row.dot(vector))
            .collect()
    }

    /// Creates a new `ParityCheckMatrix` from a list of `(row, col)` where
    /// each tuple is the position of a non zero element.
    ///
    /// # Example
    ///
    /// ```
    /// # use::believer::ParityCheckMatrix;
    /// // The parity check matrix of a 3 bits repetition code.
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     (0, 0),
    ///     (0, 1),    
    ///     (1, 1),
    ///     (1, 2),
    /// ]);
    /// ```
    pub fn new(positions: Vec<(usize, usize)>) -> Self {
        let mut column_indices = Vec::with_capacity(positions.len());
        let mut row_ranges = Vec::new();
        row_ranges.push(0);

        let mut active_row = 0;
        let mut row_lenght = 0;

        for (row, col) in positions.into_iter() {
            if row == active_row {
                row_lenght += 1;
            } else {
                while active_row < row {
                    row_ranges.push(*row_ranges.last().unwrap_or(&0) + row_lenght);
                    active_row += 1;
                }
                row_lenght = 1;
            }
            column_indices.push(col);
        }

        row_ranges.push(*row_ranges.last().unwrap_or(&0) + row_lenght);

        Self {
            row_ranges,
            column_indices,
        }
    }

    /// Returns `Some` slice of the given `row` in `self`. Returns `None` if 
    /// `row` is out of bound.
    /// 
    /// # Example 
    /// ```
    /// # use::believer::*;
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     (0, 0),
    ///     (0, 1),    
    ///     (1, 1),
    ///     (1, 2),
    /// ]);
    /// let slice = parity_check.row_slice(0).unwrap();
    /// let vector = vec![GF2::B1, GF2::B1, GF2::B0];
    /// 
    /// assert_eq!(slice.dot(&vector), GF2::B0);
    /// ```
    pub fn row_slice(&self, row: usize) -> Option<RowSlice> {
        self.row_ranges.get(row).and_then(|&row_start| {
            self.row_ranges.get(row + 1).map(|&row_end| {
                RowSlice {
                    positions: &self.column_indices[row_start..row_end]
                }
            })
        })
    }

    /// Returns an iterator that yields a slice for each row of `self`.
    /// 
    /// # Example 
    /// ```
    /// # use::believer::*;
    /// let parity_check = ParityCheckMatrix::new(vec![
    ///     (0, 0),
    ///     (0, 1),    
    ///     (1, 1),
    ///     (1, 2),
    /// ]);
    /// let mut iter = parity_check.rows_iter();
    /// 
    /// assert_eq!(iter.next(), parity_check.row_slice(0));
    /// assert_eq!(iter.next(), parity_check.row_slice(1));
    /// assert_eq!(iter.next(), None);
    /// 
    /// ```
    pub fn rows_iter(&self) -> RowsIter {
        RowsIter {
            matrix: &self,
            active_row: 0,
        }
    }
}

pub struct RowsIter<'a> {
    matrix: &'a ParityCheckMatrix,
    active_row: usize,
}

impl<'a> Iterator for RowsIter<'a> {
    type Item = RowSlice<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let slice = self.matrix.row_slice(self.active_row);
        self.active_row += 1;
        slice
    }
}

#[derive(Debug, PartialEq)]
pub struct RowSlice<'a> {
    positions: &'a [usize],
}

impl<'a> RowSlice<'a> {
    pub fn dot(&self, other: &[GF2]) -> GF2 {
        let mut total = GF2::B0;
        self.positions.iter().for_each(|&pos| {
            other.get(pos).map(|&value| total = total + value);
        });
        total
    }

    pub fn positions(&self) -> &[usize] {
        self.positions
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn dot_product() {
        let parity_check = ParityCheckMatrix::new(vec![
            (0, 0),
            (0, 1),    
            (1, 1),
            (1, 2),
        ]);   
        let bits = vec![GF2::B0, GF2::B1, GF2::B1];
        
        assert_eq!(parity_check.row_slice(0).unwrap().dot(&bits), GF2::B1);
        assert_eq!(parity_check.row_slice(1).unwrap().dot(&bits), GF2::B0);
        assert_eq!(parity_check.dot(&bits), vec![GF2::B1, GF2::B0]);
    }
}
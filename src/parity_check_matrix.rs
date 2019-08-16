pub struct ParityCheckMatrix {
    row_ranges: Vec<usize>,
    column_indices: Vec<usize>,
}

impl ParityCheckMatrix {
    /// Creates a new `ParityCheckMatrix` from a list of `(row, col)` where
    /// each tuple is the position of a non zero element.
    ///
    /// # Example
    ///
    /// ```
    /// # use::believer::ParityCheckMatrix;
    /// // The parity check matrix of a 3 bits repetition code.
    /// let parity_check = ParityCheckMatrix(vec![
    ///     (0, 0),
    ///     (0, 1),    
    ///     (1, 1),
    ///     (1, 2),
    /// ])
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
}


use crate::ParityCheckMatrix;

pub(crate) struct SparseMatrix<'a> {
    values: Vec<f64>,
    row_ranges: &'a [usize],
    column_indices: &'a [usize],
}

impl<'a> SparseMatrix<'a> {
    pub(crate) fn from_parity_check(parity_check: &'a ParityCheckMatrix, values: Vec<f64>) -> Self {
        Self {
            values,
            row_ranges: parity_check.row_ranges(),
            column_indices: parity_check.column_indices(),
        }
    }

    pub(crate) fn rows_iter(&self) -> RowsIter {
        RowsIter {
            matrix: &self,
            active_row: 0,
        }
    }

    pub(crate) fn row_slice(&self, row: usize) -> Option<RowSlice> {
        self.row_ranges.get(row).and_then(|&row_start| {
            self.row_ranges.get(row + 1).map(|&row_end| RowSlice {
                values: &self.values[row_start..row_end],
                positions: &self.column_indices[row_start..row_end],
                active: 0
            })
        })
    }
}

pub(crate) struct RowsIter<'a> {
    matrix: &'a SparseMatrix<'a>,
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

pub(crate) struct RowSlice<'a> {
    values: &'a [f64],
    positions: &'a [usize],
    active: usize,
}

impl<'a> Iterator for RowSlice<'a> {
    type Item = (&'a f64, &'a usize);

    fn next(&mut self) -> Option<Self::Item> {
        let val_pos = self.values.get(self.active).and_then(|val| {
            self.positions.get(self.active).map(|pos| (val, pos))
        });
        self.active += 1;
        val_pos
    }
}

pub(crate) struct Transposer {
    indices: Vec<usize>,
    row_ranges: Vec<usize>,
    column_indices: Vec<usize>,
}

impl Transposer {
    pub(crate) fn new(parity_check: &ParityCheckMatrix) -> Self {
        unimplemented!()
    }

    pub(crate) fn transpose(&self, matrix: &SparseMatrix) -> SparseMatrix {
        unimplemented!()
    }
}

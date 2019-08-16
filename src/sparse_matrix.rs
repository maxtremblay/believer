use crate::ParityCheckMatrix;

pub(crate) struct SparseMatrix<'a> {
    values: Vec<f64>,
    row_ranges: &'a [usize],
    column_indices: &'a [usize],
}

impl<'a> SparseMatrix<'a> {
    pub(crate) fn from_parity_check(parity_check: &ParityCheckMatrix, values: Vec<f64>) -> Self {
        unimplemented!()
    }

    pub(crate) fn rows_iter(&self) -> RowsIter {
        unimplemented!()
    }

    pub(crate) fn row_slice(&self, row: usize) -> RowSlice {
        unimplemented!()
    }

    pub(crate) fn values_iter(&self) -> std::slice::Iter<f64> {
        unimplemented!()
    }
}

pub(crate) struct RowsIter<'a> {
    matrix: &'a SparseMatrix<'a>,
    active_row: usize,
}

impl<'a> Iterator for RowsIter<'a> {
    type Item = RowSlice<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        unimplemented!()
    }
}

pub(crate) struct RowSlice<'a> {
    values: &'a [f64],
    positions: &'a [usize],
}

impl<'a> Iterator for RowSlice<'a> {
    type Item = (f64, usize);

    fn next(&mut self) -> Option<Self::Item> {
        unimplemented!()
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

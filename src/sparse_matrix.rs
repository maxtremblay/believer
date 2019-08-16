use crate::ParityCheckMatrix;

struct SparseMatrix<'a> {
    values: Vec<f64>,
    row_ranges: &'a [usize],
    column_indices: &'a [usize],
}

impl<'a> SparseMatrix<'a> {
    fn from_parity_check(parity_check: &ParityCheckMatrix) -> Self {
        unimplemented!()
    }

    fn get_transposer(&self) -> Transposer {
        unimplemented!()
    }

    fn rows_iter(&self) -> RowsIter {
        unimplemented!()
    }

    fn row_slice(&self, row: usize) -> RowSlice {
        unimplemented!()
    }

    fn values_iter(&self) -> std::slice::Iter<f64> {
        unimplemented!()
    }
}

struct RowsIter<'a> {
    matrix: &'a SparseMatrix<'a>,
    active_row: usize,
}

impl<'a> Iterator for RowsIter<'a> {
    type Item = RowSlice<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        unimplemented!()
    }
}

struct RowSlice<'a> {
    values: &'a [f64],
    positions: &'a [usize],
}


impl <'a> Iterator for RowSlice<'a> {
    type Item = (f64, usize);

    fn next(&mut self) -> Option<Self::Item> {
        unimplemented!()
    }
}

struct Transposer {
    indices: Vec<usize>,
    row_ranges: Vec<usize>,
    column_indices: Vec<usize>,
}

impl Transposer {
    fn transpose(&self, matrix: &SparseMatrix) -> SparseMatrix {
        unimplemented!()
    }
}
use crate::ParityCheckMatrix;
use itertools::Itertools;
use std::cmp::Ordering;

#[derive(Debug)]
pub(crate) struct SparseMatrix<'a> {
    values: Vec<f64>,
    row_ranges: &'a [usize],
    column_indices: &'a [usize],
}

impl<'a> SparseMatrix<'a> {
    pub(crate) fn from_parity_check(parity_check: &'a ParityCheckMatrix, values: Vec<f64>) -> Self {
        if parity_check.len() != values.len() {
            panic!("wrong number of values");
        }
        Self {
            values,
            row_ranges: parity_check.check_ranges(),
            column_indices: parity_check.bit_indices(),
        }
    }

    pub(crate) fn rows_iter(&self) -> RowsIter {
        RowsIter {
            matrix: &self,
            active_row: 0,
        }
    }

    pub(crate) fn row_slice(&self, row: usize) -> Option<RowSlice> {
        println!("*****");
        self.row_ranges.get(row).and_then(|&row_start| {
            self.row_ranges.get(row + 1).map(|&row_end| {
                println!("Row start and end: {} & {}", row_start, row_end);
                println!("Val, col,: {:?}, {:?}", self.values, self.column_indices);
                RowSlice {
                    values: &self.values[row_start..row_end],
                    positions: &self.column_indices[row_start..row_end],
                    active: 0,
                }
            })
        })
    }

    pub(crate) fn values(&self) -> &[f64] {
        &self.values
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

#[derive(Debug)]
pub(crate) struct RowSlice<'a> {
    values: &'a [f64],
    positions: &'a [usize],
    active: usize,
}

impl<'a> Iterator for RowSlice<'a> {
    type Item = (&'a f64, &'a usize);

    fn next(&mut self) -> Option<Self::Item> {
        let val_pos = self
            .values
            .get(self.active)
            .and_then(|val| self.positions.get(self.active).map(|pos| (val, pos)));
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
        let mut indices = Vec::with_capacity(parity_check.n_bits());
        let mut column_indices = Vec::with_capacity(parity_check.n_bits());
        let mut row_ranges = Vec::new();
        row_ranges.push(0);

        let mut active_col = 0;
        let mut row_lenght = 0;

        parity_check
            .positions_iter()
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

        Transposer {
            indices,
            row_ranges,
            column_indices,
        }
    }

    pub(crate) fn transpose(&self, matrix: &SparseMatrix) -> SparseMatrix {
        SparseMatrix {
            values: self
                .indices
                .iter()
                .map(|idx| matrix.values()[*idx])
                .collect(),
            row_ranges: &self.row_ranges,
            column_indices: &self.column_indices,
        }
    }
}

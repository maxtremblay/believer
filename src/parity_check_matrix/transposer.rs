use super::{Check, CheckView};

// A tool to help transpose a parity check matrix.
pub(super) struct Transposer {
    checks: Vec<Check>,
    active_check: usize,
}

impl Transposer {
    pub(super) fn new() -> Self {
        Self { checks: Vec::new(), active_check: 0 }
    }

    pub(super) fn insert_bits_from(&mut self, check: CheckView) {
        check.iter().for_each(|bit| self.insert(*bit));
        self.active_check += 1;
    }

    fn insert(&mut self, bit: usize) {
        self.checks[bit].push(self.active_check);
    }

    pub(super) fn get_checks(self) -> Vec<Check> {
        self.checks
    }

}
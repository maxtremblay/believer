use crate::GF2;
use std::cmp::min;

/// Represent one check in a parity check matrix.
pub type Check = Vec<usize>;

/// A reference to a check.
pub type CheckSlice<'a> = &'a [usize];

pub fn get_dot_product(left_check: CheckSlice, right_check: CheckSlice) -> GF2 {
    Dotter::from_checks(left_check, right_check).get_dot_product()
}

struct Dotter<'a> {
    left_check: CheckSlice<'a>,
    right_check: CheckSlice<'a>,
    left_index: usize,
    right_index: usize,
    sum: GF2,
}

impl<'a> Dotter<'a> {
    fn from_checks(left_check: CheckSlice<'a>, right_check: CheckSlice<'a>) -> Self {
        Self {
            left_check,
            right_check,
            left_index: 0,
            right_index: 0,
            sum: GF2::B0,
        }
    }

    fn get_dot_product(mut self) -> GF2 {
        while self.has_not_reach_the_end_of_a_check() {
            self.update_sum_from_active_bits();
            self.go_to_next_bits();
        }
        self.sum
    }

    fn has_not_reach_the_end_of_a_check(&self) -> bool {
        self.left_index < self.left_check.len() && self.right_index < self.right_check.len()
    }

    fn update_sum_from_active_bits(&mut self) {
        if self.get_active_left_bit() == self.get_active_right_bit() {
            self.sum = self.sum + GF2::B1;
        }
    }

    fn get_active_left_bit(&self) -> usize {
        self.left_check[self.left_index]
    }

    fn get_active_right_bit(&self) -> usize {
        self.right_check[self.right_index]
    }

    fn go_to_next_bits(&mut self) {
        if self.get_active_left_bit() == self.get_active_right_bit() {
            self.go_to_next_left_bit();
            self.go_to_next_right_bit();
        } else if self.get_active_left_bit() <= self.get_active_right_bit() {
            self.go_to_next_left_bit();
        } else {
            self.go_to_next_right_bit();
        }
    }

    fn go_to_next_left_bit(&mut self) {
        self.left_index += 1;
    }

    fn go_to_next_right_bit(&mut self) {
        self.right_index += 1;
    }
}

pub fn get_bitwise_sum(left_check: CheckSlice, right_check: CheckSlice) -> Check {
    BitwiseSummer::from_checks(left_check, right_check).get_bitwise_sum()
}

// pub(super) fn get_bitwise_sum(left_check: CheckSlice, right_check: CheckSlice) -> Check {
//     let mut sum = Vec::new();
//     BinarySparseOperator::with_operation(|left_bit, right_bit| {
//         if left_bit != right_bit {
//             println!("Left {} // Right {}", left_bit, right_bit);
//             sum.push(min(left_bit, right_bit));
//             println!("Sum: {:?}", sum);
//         }
//     })
//     .apply_over_checks(left_check, right_check);
//     println!("Final sum: {:?}", sum);
//     sum
// }

struct BitwiseSummer<'a> {
    left_check: CheckSlice<'a>,
    right_check: CheckSlice<'a>,
    left_index: usize,
    right_index: usize,
    sum: Check,
}

impl<'a> BitwiseSummer<'a> {
    fn from_checks(left_check: CheckSlice<'a>, right_check: CheckSlice<'a>) -> Self {
        Self {
            left_check,
            right_check,
            left_index: 0,
            right_index: 0,
            sum: Check::new(),
        }
    }

    fn get_bitwise_sum(mut self) -> Check {
        while self.has_not_reach_the_end_of_a_check() {
            self.update_sum_from_active_bits();
            self.go_to_next_bits();
        }
        self.fill_sum_with_leftovers();
        self.sum
    }

    fn has_not_reach_the_end_of_a_check(&self) -> bool {
        self.left_index < self.left_check.len() && self.right_index < self.right_check.len()
    }

    fn update_sum_from_active_bits(&mut self) {
        let left_bit = self.get_active_left_bit();
        let right_bit = self.get_active_right_bit();
        if left_bit != right_bit {
            self.sum.push(min(left_bit, right_bit));
        }
    }

    fn get_active_left_bit(&self) -> usize {
        self.left_check[self.left_index]
    }

    fn get_active_right_bit(&self) -> usize {
        self.right_check[self.right_index]
    }

    fn go_to_next_bits(&mut self) {
        if self.get_active_left_bit() == self.get_active_right_bit() {
            self.go_to_next_left_bit();
            self.go_to_next_right_bit();
        } else if self.get_active_left_bit() <= self.get_active_right_bit() {
            self.go_to_next_left_bit();
        } else {
            self.go_to_next_right_bit();
        }
    }

    fn go_to_next_left_bit(&mut self) {
        self.left_index += 1;
    }

    fn go_to_next_right_bit(&mut self) {
        self.right_index += 1;
    }

    fn fill_sum_with_leftovers(&mut self) {
        if self.left_index < self.left_check.len() {
            self.sum
                .extend_from_slice(&self.left_check[self.left_index..])
        } else {
            self.sum
                .extend_from_slice(&self.right_check[self.right_index..])
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::GF2;

    #[test]
    fn dot_product() {
        let left_check = vec![0, 1, 3, 5, 8];
        let right_check = vec![1, 2, 5, 6, 8];
        assert_eq!(get_dot_product(&left_check, &right_check), GF2::B1);

        let left_check = vec![0, 1, 5, 6, 8];
        let right_check = vec![1, 2, 5, 6, 8];
        assert_eq!(get_dot_product(&left_check, &right_check), GF2::B0);

        let left_check = vec![0, 1, 3, 5, 8];
        let right_check = vec![];
        assert_eq!(get_dot_product(&left_check, &right_check), GF2::B0);
    }

    #[test]
    fn bitwise_sum() {
        let left_check = vec![0, 1, 3, 5, 8];
        let right_check = vec![1, 2, 5, 6, 8];
        assert_eq!(get_bitwise_sum(&left_check, &right_check), vec![0, 2, 3, 6]);

        let left_check = vec![0, 1, 5, 6, 8];
        let right_check = vec![1, 2, 5, 6, 8];
        assert_eq!(get_bitwise_sum(&left_check, &right_check), vec![0, 2]);

        let left_check = vec![0, 1, 3, 5, 8];
        let right_check = vec![];
        assert_eq!(
            get_bitwise_sum(&left_check, &right_check),
            vec![0, 1, 3, 5, 8]
        );
    }
}

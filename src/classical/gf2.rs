//! A simple implementation of a binary field with addition, multiplication and comparison.
//!
//! # Rules
//!
//! ```
//! # use believer::GF2;
//! assert_eq!(GF2::B0 + GF2::B0, GF2::B0);
//! assert_eq!(GF2::B0 + GF2::B1, GF2::B1);
//! assert_eq!(GF2::B1 + GF2::B0, GF2::B1);
//! assert_eq!(GF2::B1 + GF2::B1, GF2::B0);
//!
//! assert_eq!(GF2::B0 * GF2::B0, GF2::B0);
//! assert_eq!(GF2::B1 * GF2::B1, GF2::B1);
//! assert_eq!(GF2::B0 * GF2::B1, GF2::B0);
//! assert_eq!(GF2::B1 * GF2::B0, GF2::B0);
//! ```
use num::{One, Zero};
use std::ops::{Add, Mul};

#[repr(u8)]
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub enum GF2 {
    B0,
    B1,
}

impl GF2 {
    /// Converts a `0_u8` or `1_u8` to the corresponding `GF2` element.
    ///
    /// # Panics
    /// If the input is neither `0_u8` or a `1_u8`.
    pub fn from_u8(b: u8) -> Self {
        if b == 0 {
            GF2::B0
        } else if b == 1 {
            GF2::B1
        } else {
            panic!("Input must be 0 or 1.");
        }
    }
}

impl Add for GF2 {
    type Output = Self;
    fn add(self, other: Self) -> Self::Output {
        if self == other {
            GF2::B0
        } else {
            GF2::B1
        }
    }
}

impl Mul for GF2 {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        if self == GF2::B1 && other == GF2::B1 {
            GF2::B1
        } else {
            GF2::B0
        }
    }
}

impl One for GF2 {
    fn one() -> GF2 {
        GF2::B1
    }
}

impl Zero for GF2 {
    fn is_zero(&self) -> bool {
        self == &GF2::zero()
    }

    fn zero() -> GF2 {
        GF2::B0
    }
}

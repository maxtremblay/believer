use std::ops::Mul;

/// The Pauli operators.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy)]
pub enum Pauli {
    I,
    X,
    Y,
    Z,
}

impl Pauli {
    /// Returns the GF4 representation of the Pauli operator.
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let paulis = vec![Pauli::I, Pauli::X, Pauli::Y, Pauli::Z];
    /// let gf4_paulis = paulis.iter().map(|p| p.as_gf4()).collect();
    /// assert_eq!(gf4_paulis[0], (0, 0)); // I
    /// assert_eq!(gf4_paulis[1], (1, 0)); // X
    /// assert_eq!(gf4_paulis[2], (1, 1)); // Y
    /// assert_eq!(gf4_paulis[3], (0, 1)); // Z
    /// ```
    pub fn as_gf4(&self) -> (usize, usize) {
        match self {
            Self::I => (0, 0),
            Self::X => (1, 0),
            Self::Y => (1, 1),
            Self::Z => (0, 1),
        }
    }

    /// Returns `+1` if `self` and `other` are commuting and `-1` otherwise.
    /// 
    /// # Example 
    /// 
    /// ```
    /// # use believer::*;
    /// let paulis = vec![Pauli::I, Pauli::X, Pauli::Y, Pauli::Z];
    /// // Every pauli commute with identity
    /// paulis.iter().for_each(|p| assert_eq!(p.commutator_with(Pauli::I), 1));
    /// 
    /// // Every pauli commute with itself.
    /// paulis.iter().for_each(|p| assert_eq!(p.commutator_with(*p), 1));
    /// 
    /// // X, Y and Z anticommute
    /// [Pauli::Y, Pauli::Z].iter().for_each(|p| assert_eq!(p.commutator_with(Pauli::X), -1));
    /// assert_eq!(Pauli::Y.commutator_with(Pauli::Z), -1);
    /// ```
    pub fn commutator_with(self, other: Self) -> i32 {
        match (self, other) {
            (Self::I, _) => 1,
            (_, Self::I) => 1,
            (a, b) => if a == b { 1 } else { -1 }
        }
    }
}

/// Returns the group product of 2 pauli operators.
///
/// # Example
///
/// ```
/// # use believer::*;
/// [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z].into_iter().for_each(|&p| {
///     // All paulis square to identity.
///     assert_eq!(p * p, Pauli::I);
///
///     // Any pauli multiply by identity is itself.
///     assert_eq!(p * Pauli::I, p);
///     assert_eq!(Pauli::I * p, p);     
/// });
///
/// // We have the following commution relations.
/// assert_eq!(Pauli::X * Pauli::Y, Pauli::Z);
/// assert_eq!(Pauli::Y * Pauli::X, Pauli::Z);
///
/// assert_eq!(Pauli::Y * Pauli::Z, Pauli::X);
/// assert_eq!(Pauli::Z * Pauli::Y, Pauli::X);
///
/// assert_eq!(Pauli::Z * Pauli::X, Pauli::Y);
/// assert_eq!(Pauli::X * Pauli::Z, Pauli::Y);
/// ```
impl Mul for Pauli {
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        match (self, other) {
            (Self::I, a) => a,
            (Self::X, Self::Y) => Self::Z,
            (Self::Y, Self::Z) => Self::X,
            (Self::Z, Self::X) => Self::Y,
            (a, b) => {
                if a == b {
                    Pauli::I
                } else {
                    b * a
                }
            }
        }
    }
}

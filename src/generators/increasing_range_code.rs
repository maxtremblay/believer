use super::{CodeGenerator, RandomCheckGenerator};
use crate::ParityCheckMatrix;
use rand::{thread_rng, Rng};

/// A generator for increasing range code.
///
/// An increasing range code is a ldpc code where half the checks fit in a box of length
/// `initial_range`, a quarter in a box of length `2 * initial_range`, an eight in a box of length
/// `4 * initial_range` and so on.
///
/// The way to create a `IncreasingRangeCodeGenerator` is to use the `IRCodeGenBuilder`
/// interface.
///
/// # Example
///
/// ```
/// # use believer::*;
/// let n_bits = 4;
/// let n_checks = 3;
/// let generator = IRCodeGenBuilder::new(n_bits, n_checks)
///     .with_max_bit_degree(3)
///     .with_max_check_degree(4)
///     .build();
///
/// let code = generator.generate();
///
/// assert_eq!(code.n_bits(), 4);
/// assert_eq!(code.n_checks(), 3);
///
/// // Bit degrees are at most 3.
/// code.bit_degrees().iter().for_each(|&d| assert_eq!(d <= 3, true));
///
/// // Check degrees are at most 4.
/// code.check_degrees().iter().for_each(|&d| assert_eq!(d <= 4, true));
/// ```
pub struct IncreasingRangeCodeGenerator {
    initial_range: usize,
    n_bits: usize,
    n_checks: usize,
    max_bit_degree: usize,
    max_check_degree: usize,
    minimal_girth: usize,
}

impl IncreasingRangeCodeGenerator {
    // Get a random range in which to pick the next check.
    fn get_random_range(&self, range_length: usize) -> Vec<usize> {
        if range_length >= self.n_bits {
            (0..self.n_bits).collect()
        } else {
            let mut rng = thread_rng();
            let range_start = rng.gen_range(0, self.n_bits - range_length);
            (range_start..range_start + range_length).collect()
        }
    }
}

impl CodeGenerator for IncreasingRangeCodeGenerator {
    /// Returns a code generated by the generator.
    fn generate(&self) -> ParityCheckMatrix {
        let mut check_generator =
            RandomCheckGenerator::new(vec![self.max_bit_degree; self.n_bits], self.minimal_girth);
        let mut checks = Vec::with_capacity(self.n_checks);
        let mut rng = thread_rng();
        let mut range_length = self.initial_range;
        let mut last_update = 0;
        for c in 0..self.n_checks {
            check_generator.over_bits(self.get_random_range(range_length));
            if let Some(check) = check_generator.generate(self.max_check_degree) {
                checks.push(check)
            }
            if (self.n_checks - last_update) <= 2 * c {
                last_update = c;
                range_length *= 2;
            }
        }
        ParityCheckMatrix::new(checks, self.n_bits)
    }
}

/// An interface to build `IncreasingRangeCodeGenerator`.
///
/// All parameters are optional except the number of bits and checks. It is possible to set the
/// maximal degrees of both bits and checks, the initial range (the first half of the checks must
/// fit in a box of length `initial_range`, the next quarter in a box of length `2 * initial_range`
/// and so on) and the minimal girth.
///
/// # Example
///
/// ```
/// # use believer::*;
/// // Build a code generator by setting all parameters.
///
/// let n_bits = 16;
/// let n_checks = 12;
/// let code_generator = IRCodeGenBuilder::new(n_bits, n_checks)
///     .with_max_bit_degree(5)
///     .with_max_check_degree(5)
///     .with_initial_range(10)
///     .with_minimal_girth(8)
///     .build();
///
/// // Build a code generator by setting the size and the degrees. This is equivalent to an
/// // standard LDPC code with no locality contraints.
///
/// let bit_deg = 3;
/// let check_deg = 4;
/// let code_generator = IRCodeGenBuilder::new(n_bits, n_checks)
///     .with_max_degrees(bit_deg, check_deg)
///     .build();
///
/// // It is possible to reuse a builder.
///
/// let mut builder = IRCodeGenBuilder::new(n_bits, n_checks);
/// let generator_a = builder.build();
///
/// // Set the maximal degrees.
/// builder.with_max_degrees(bit_deg, check_deg);
/// let generator_b = builder.build();
///
/// // Set the other parameters. The maximal degrees are still the same.
/// builder.with_initial_range(10).with_minimal_girth(8);
/// let generator_c = builder.build();
///
/// // Change the maximal bit degree.
/// builder.with_max_bit_degree(5);
/// let generator_d = builder.build();
/// ```
pub struct IRCodeGenBuilder {
    n_bits: usize,
    n_checks: usize,
    initial_range: Option<usize>,
    max_bit_degree: Option<usize>,
    max_check_degree: Option<usize>,
    minimal_girth: Option<usize>,
}

impl IRCodeGenBuilder {
    // **************
    // Public methods
    // **************

    /// Builds an `IncreasingRangeCodeGenerator` from the given parameters.
    pub fn build(&self) -> IncreasingRangeCodeGenerator {
        IncreasingRangeCodeGenerator {
            initial_range: self.initial_range(),
            n_bits: self.n_bits(),
            n_checks: self.n_checks(),
            max_bit_degree: self.max_bit_degree(),
            max_check_degree: self.max_check_degree(),
            minimal_girth: self.minimal_girth(),
        }
    }

    /// Creates an `IRCodeGenBuilder` that will build `IncreasingRangeCodeGenerator` over `n_bits`
    /// and `n_checks`.
    pub fn new(n_bits: usize, n_checks: usize) -> Self {
        Self {
            n_bits,
            n_checks,
            initial_range: None,
            max_bit_degree: None,
            max_check_degree: None,
            minimal_girth: None,
        }
    }

    /// Fixes the `max_bit_degree` of the code generator.
    ///
    /// If not fixed, there will be no restriction on the degree of the bits.
    pub fn with_max_bit_degree(&mut self, max_bit_degree: usize) -> &mut Self {
        self.max_bit_degree = Some(max_bit_degree);
        self
    }

    /// Fixes the `max_check_degree` of the code generator.
    ///
    /// If not fixed, there will be no restriction on the degree of the checks.
    pub fn with_max_check_degree(&mut self, max_check_degree: usize) -> &mut Self {
        self.max_check_degree = Some(max_check_degree);
        self
    }

    /// Fixes the `max_bit_degree` and the `max_check_degree` of the code generator.
    ///
    /// If not fixed, there will be no restriction on the degree of the bits and the checks.
    pub fn with_max_degrees(
        &mut self,
        max_bit_degree: usize,
        max_check_degree: usize,
    ) -> &mut Self {
        self.max_bit_degree = Some(max_bit_degree);
        self.max_check_degree = Some(max_check_degree);
        self
    }

    /// Fixes the `initial_range` of the code generator.
    ///
    /// If not fixed, it will be set to the maximal check degree.
    pub fn with_initial_range(&mut self, initial_range: usize) -> &mut Self {
        self.initial_range = Some(initial_range);
        self
    }

    /// Fixes the `minimal_girth` of the code generator.
    ///
    /// If not fixed, there will be no restriction on the minimal girth.
    pub fn with_minimal_girth(&mut self, minimal_girth: usize) -> &mut Self {
        self.minimal_girth = Some(minimal_girth);
        self
    }

    // ***************
    // Private methods
    // ***************

    // Returns the initial range if fixed or default to the number of bits if not.
    fn initial_range(&self) -> usize {
        self.initial_range.unwrap_or_else(|| self.n_bits())
    }

    // Returns the number of bits.
    fn n_bits(&self) -> usize {
        self.n_bits
    }

    // Returns the number of checks.
    fn n_checks(&self) -> usize {
        self.n_checks
    }

    // Returns the maximal bit degree if fixed or default to the number of checks if not.
    fn max_bit_degree(&self) -> usize {
        self.max_bit_degree.unwrap_or(self.n_checks())
    }

    // Returns the maximal check degree if fixed or default to the number of bits if not.
    fn max_check_degree(&self) -> usize {
        self.max_check_degree.unwrap_or(self.n_bits())
    }

    // Returns the minimal girth if fixed or default to 0 if not.
    fn minimal_girth(&self) -> usize {
        self.minimal_girth.unwrap_or(0)
    }
}

// #[cfg(test)]
// mod test {
//     use super::*;

//     #[test]
//     fn general() {
//         let IncreasingRangeCodeGenerator::new(4, 16, 12, 4, 4, 0);

//     }
// }

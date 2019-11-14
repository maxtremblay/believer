use super::CodeGenerator;
use crate::ParityCheckMatrix;
use rand::Rng;

pub struct HierarchicalCodeGenerator {
    bit_degree: usize,
    check_degree: usize,
    n_layers: usize,
    block_length: usize,
}

impl HierarchicalCodeGenerator {
    // ************
    // Construction
    // ************

    /// Creates a `HierarchicalCodeGenerator` that will generate an empty code.
    /// 
    /// It is intended to be chain with other construction functions.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use believer::{HierarchicalCodeGenerator, ParityCheckMatrix};
    /// 
    /// ```
    pub fn new() -> Self {
        Self {
            bit_degree: 0,
            check_degree: 0,
            n_layers: 0,
            block_length: 0,
        }
    }

    pub fn with_bit_degree(mut self, degree: usize) -> Self {
        self.bit_degree = degree;
        self
    }

    pub fn with_check_degree(mut self, degree: usize) -> Self {
        self.check_degree = degree;
        self
    }

    pub fn with_n_layers(mut self, n_layers: usize) -> Self {
        self.n_layers = n_layers;
        self
    }

    pub fn with_block_length(mut self, length: usize) -> Self {
        self.block_length = length;
        self
    }

    // *******
    // Setters
    // *******

    pub fn set_bit_degree(&mut self, degree: usize) -> &mut Self {
        self.bit_degree = degree;
        self
    }

    pub fn set_check_degree(&mut self, degree: usize) -> &mut Self {
        self.check_degree = degree;
        self
    }

    pub fn set_n_layers(&mut self, n_layers: usize) -> &mut Self {
        self.n_layers = n_layers;
        self
    }

    pub fn set_block_length(&mut self, length: usize) -> &mut Self {
        self.block_length = length;
        self
    }

    // *********************
    // Code generation tools
    // *********************

    pub fn will_generate_an_empty_code(&self) -> bool {
        self.bit_degree == 0
            || self.check_degree == 0
            || self.n_layers == 0
            || self.block_length == 0
    }
}

impl CodeGenerator for HierarchicalCodeGenerator {
    fn generate<R: Rng + ?Sized>(&self, rng: &mut R) -> ParityCheckMatrix {
        if self.will_generate_an_empty_code() {
            ParityCheckMatrix::empty()
        } else {
            unimplemented!()
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn default_constructor_generate_empty_code() {
        let generator = HierarchicalCodeGenerator::new();
        assert!(generator.will_generate_an_empty_code());
        assert_eq!(generator.generate(&mut thread_rng()), ParityCheckMatrix::empty());
    }
}

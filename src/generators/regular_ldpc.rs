use super::CodeGenerator;
use super::RandomCheckGenerator;
use crate::ParityCheckMatrix;

pub struct RegularLDPCCodeGenerator {
    bit_degree: usize,
    check_degree: usize,
    scale: usize,
    minimal_girth: usize,
}

impl CodeGenerator for RegularLDPCCodeGenerator {
    fn generate(&self) -> ParityCheckMatrix {
        let mut check_generator =
            RandomCheckGenerator::new(vec![self.bit_degree; self.n_bits()], self.minimal_girth);
        let mut checks = Vec::with_capacity(self.n_checks());
        for _ in 0..self.n_checks() {
            if let Some(check) = check_generator.generate(self.check_degree) {
                checks.push(check);
            }
        }
        ParityCheckMatrix::new(checks, self.n_bits())
    }
}

impl RegularLDPCCodeGenerator {
    pub fn n_bits(&self) -> usize {
        self.scale * self.check_degree
    }

    pub fn n_checks(&self) -> usize {
        self.scale * self.bit_degree
    }

    pub fn new(bit_degree: usize, check_degree: usize, scale: usize, minimal_girth: usize) -> Self {
        Self {
            bit_degree,
            check_degree,
            scale,
            minimal_girth,
        }
    }
}
use super::CodeGenerator;
use super::random_check::Generator as RandomCheckGenerator;
use crate::ParityCheckMatrix;
use rand::Rng;

pub struct RegularLDPCCodeGenerator {
    bit_degree: usize,
    check_degree: usize,
    scale: usize,
    minimal_girth: usize,
}

impl CodeGenerator for RegularLDPCCodeGenerator {
    fn generate<R: Rng + ?Sized>(&self, rng: &mut R) -> ParityCheckMatrix {
        let mut check_generator =
            RandomCheckGenerator::new()
                .with_max_bit_degrees(self.bit_degree)
                .with_n_bits(self.n_bits())
                .with_minimal_girth(self.minimal_girth);   
        let mut checks = Vec::with_capacity(self.n_checks());
        for _ in 0..self.n_checks() {
            if let Some(check) = check_generator.generate(self.check_degree, rng) {
                checks.push(check);
            }
        }
        ParityCheckMatrix::new(checks, self.n_bits())
    }
}

impl RegularLDPCCodeGenerator {
    pub fn n_bits(&self) -> usize {
        (self.scale * self.check_degree) as usize
    }

    pub fn n_checks(&self) -> usize {
        (self.scale * self.bit_degree) as usize
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

use super::random_checks::Generator as RandomCheckGenerator;
use super::CodeGenerator;
use crate::ParityCheckMatrix;
use rand::Rng;

pub struct RegularLDPCCodeGenerator {
    bit_degree: usize,
    check_degree: usize,
    scale: usize,
    minimal_girth: usize,
}

impl CodeGenerator for RegularLDPCCodeGenerator {
    fn generate_with_rng<R: Rng>(&self, rng: &mut R) -> ParityCheckMatrix {
        let mut check_generator = RandomCheckGenerator::with_n_bits(self.n_bits())
            .with_random_number_generator(rng);
        check_generator
            .set_maximal_bit_degree(self.bit_degree)
            .set_minimal_girth(self.minimal_girth);
        let mut checks = Vec::with_capacity(self.n_checks());
        for _ in 0..self.n_checks() {
            let check = check_generator
                .set_target_check_degree(self.check_degree)
                .get_random_check();
            if let Some(check) = check { checks.push(check) }
        }
        ParityCheckMatrix::with_n_bits(self.n_bits()).with_checks(checks)
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

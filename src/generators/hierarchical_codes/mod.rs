use super::random_checks::Generator as CheckGenerator;
use super::CodeGenerator;
use crate::ParityCheckMatrix;

mod code_constructor;
use code_constructor::CodeConstructor;

use rand::Rng;

pub struct Generator {
    bit_degree: usize,
    check_degree: usize,
    n_layers: u32,
    initial_block_length: Option<usize>,
    n_checks_per_block: usize,
    n_blocks_per_layer: usize,
}

impl Generator {
    // ***** Construction *****

    pub fn new() -> Self {
        Self {
            bit_degree: 0,
            check_degree: 0,
            n_layers: 1,
            initial_block_length: None,
            n_checks_per_block: 1,
            n_blocks_per_layer: 1,
        }
    }

    pub fn with_bit_and_check_degree(bit_degree: usize, check_degree: usize) -> Self {
        Self::new()
            .with_bit_degree(bit_degree)
            .with_check_degree(check_degree)
    }

    fn with_bit_degree(mut self, degree: usize) -> Self {
        self.bit_degree = degree;
        self
    }

    fn with_check_degree(mut self, degree: usize) -> Self {
        self.check_degree = degree;
        self
    }

    // ***** Setters *****

    pub fn set_bit_degree(&mut self, bit_degree: usize) -> &mut Self {
        self.bit_degree = bit_degree;
        self
    }

    pub fn set_check_degree(&mut self, check_degree: usize) -> &mut Self {
        self.check_degree = check_degree;
        self
    }

    pub fn set_n_layers(&mut self, n_layers: u32) -> &mut Self {
        self.n_layers = n_layers;
        self
    }

    pub fn set_initial_block_length(&mut self, length: usize) -> &mut Self {
        self.initial_block_length = Some(length);
        self
    }

    pub fn set_n_checks_per_block(&mut self, n_checks: usize) -> &mut Self {
        self.n_checks_per_block = n_checks;
        self
    }

    pub fn set_n_blocks_per_layer(&mut self, n_blocks: usize) -> &mut Self {
        self.n_blocks_per_layer = n_blocks;
        self
    }

    // ***** Getters *****

    pub fn get_n_bits(&self) -> usize {
        self.get_initial_block_length().pow(self.n_layers as u32)
    }

    pub fn get_n_checks(&self) -> usize {
        self.n_checks_per_block.pow(self.n_layers as u32) as usize
    }

    fn get_initial_block_length(&self) -> usize {
        self.initial_block_length.unwrap_or(self.check_degree)
    }

    // ***** Methods for trait CodeGenerator *****

    fn will_generate_an_empty_code(&self) -> bool {
        self.bit_degree == 0
            || self.check_degree == 0
            || self.n_layers == 0
            || self.get_initial_block_length() == 0
            || self.n_checks_per_block == 0
            || self.n_blocks_per_layer == 0
    }

    fn generate_non_empty_code_with_rng<R: Rng>(&self, rng: R) -> ParityCheckMatrix {
        self.get_code_constructor_with_rng(rng).get_random_code()
    }

    fn get_code_constructor_with_rng<R: Rng>(&self, rng: R) -> CodeConstructor<R> {
        CodeConstructor {
            check_generator: self.initialize_check_generator_with_rng(rng),
            checks: Vec::with_capacity(self.get_n_checks()),
            active_layer: 0,
            n_layers: self.n_layers,
            n_bits: self.get_n_bits(),
            n_blocks_per_layer: self.n_blocks_per_layer,
            n_checks_per_block: self.n_checks_per_block,
        }
    }

    fn initialize_check_generator_with_rng<R: Rng>(&self, rng: R) -> CheckGenerator<R> {
        let mut generator =
            CheckGenerator::with_n_bits(self.get_n_bits()).with_random_number_generator(rng);
        generator
            .set_maximal_bit_degree(self.bit_degree)
            .set_target_check_degree(self.check_degree)
            .allow_checks_of_degree_at_least(2);
        generator
    }
}

impl CodeGenerator for Generator {
    fn generate_with_rng<R: Rng>(&self, rng: R) -> ParityCheckMatrix {
        if self.will_generate_an_empty_code() {
            ParityCheckMatrix::new()
        } else {
            self.generate_non_empty_code_with_rng(rng)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn default_constructor_generate_empty_code() {
        let generator = Generator::new();
        assert_eq!(generator.generate(), ParityCheckMatrix::new());
    }

    #[test]
    fn generate_single_check_code() {
        let generator = Generator::with_bit_and_check_degree(3, 4);

        let expected_code = ParityCheckMatrix::with_n_bits(4).with_checks(vec![vec![0, 1, 2, 3]]);
        assert_eq!(generator.generate(), expected_code);
    }

    #[test]
    fn generate_single_layer_code() {
        let mut generator = Generator::with_bit_and_check_degree(3, 3);
        generator
            .set_initial_block_length(8)
            .set_n_checks_per_block(6);
        let code = generator.generate_with_rng(ChaCha8Rng::seed_from_u64(10));
        println!("{}", code);
    }
}

use super::random_check::Generator as RandomCheckGenerator;
use super::CodeGenerator;
use crate::ParityCheckMatrix;
use rand::Rng;

pub struct HierarchicalCodeGeneratorBuilder {
    bit_degree: usize,
    check_degree: usize,
    n_layers: Option<usize>,
    initial_block_length: Option<usize>,
    n_checks_per_block: Option<usize>,
    n_blocks_per_layer: Option<usize>,
}

impl HierarchicalCodeGeneratorBuilder {
    // ************
    // Construction
    // ************

    pub fn new() -> Self {
        Self {
            bit_degree: 0,
            check_degree: 0,
            n_layers: None,
            initial_block_length: None,
            n_checks_per_block: None,
            n_blocks_per_layer: None,
        }
    }

    pub fn with_bit_and_check_degree(mut self, bit_degree: usize, check_degree: usize) -> Self {
        self.bit_degree = bit_degree;
        self.check_degree = check_degree;
        self
    }

    pub fn with_n_layers(mut self, n_layers: usize) -> Self {
        self.n_layers = Some(n_layers);
        self
    }

    pub fn with_initial_block_length(mut self, length: usize) -> Self {
        self.initial_block_length = Some(length);
        self
    }

    pub fn with_n_checks_per_block(mut self, n_checks: usize) -> Self {
        self.n_checks_per_block = Some(n_checks);
        self
    }

    pub fn with_n_blocks_per_layer(mut self, n_blocks: usize) -> Self {
        self.n_blocks_per_layer = Some(n_blocks);
        self
    }

    // ********
    // Building
    // ********

    pub fn build(&self) -> HierarchicalCodeGenerator {
        HierarchicalCodeGenerator {
            bit_degree: self.bit_degree,
            check_degree: self.check_degree,
            n_layers: self.get_n_layers(),
            initial_block_length: self.get_initial_block_length(),
            n_checks_per_block: self.get_n_checks_per_block(),
            n_blocks_per_layer: self.get_n_blocks_for_layer(),
        }
    }

    fn get_n_layers(&self) -> usize {
        self.n_layers.unwrap_or(1)
    }

    fn get_initial_block_length(&self) -> usize {
        self.initial_block_length.unwrap_or(self.check_degree)
    }

    fn get_n_checks_per_block(&self) -> usize {
        self.n_checks_per_block.unwrap_or(1)
    }

    fn get_n_blocks_for_layer(&self) -> usize {
        self.n_blocks_per_layer.unwrap_or(1)
    }
}

pub struct HierarchicalCodeGenerator {
    bit_degree: usize,
    check_degree: usize,
    n_layers: usize,
    initial_block_length: usize,
    n_checks_per_block: usize,
    n_blocks_per_layer: usize,
}

impl HierarchicalCodeGenerator {
    // ************
    // Construction
    // ************

    /// Creates a `HierarchicalCodeGenerator` that will generate an empty code.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::hierarchical_codes::HierarchicalCodeGenerator;
    /// let generator = HierarchicalCodeGenerator::new()
    /// ```
    pub fn new() -> Self {
        Self {
            bit_degree: 0,
            check_degree: 0,
            n_layers: 1,
            initial_block_length: 0,
            n_checks_per_block: 1,
            n_blocks_per_layer: 1,
        }
    }

    // *********************
    // Code generation tools
    // *********************

    fn will_generate_an_empty_code(&self) -> bool {
        self.bit_degree == 0
            || self.check_degree == 0
            || self.n_layers == 0
            || self.initial_block_length == 0
            || self.n_checks_per_block == 0
            || self.n_blocks_per_layer == 0
    }

    fn generate_non_empty_code<R: Rng + ?Sized>(&self, rng: &mut R) -> ParityCheckMatrix {
        let mut check_generator = self.initialize_check_generator();
        let mut checks = Vec::with_capacity(self.get_n_checks());
        for layer in 0..self.n_layers {
            checks.extend(self.generate_layer_checks(layer, &mut check_generator, rng));
        }
        ParityCheckMatrix::with_n_bits(self.get_n_bits()).with_checks(checks)
    }

    fn initialize_check_generator(&self) -> RandomCheckGenerator {
        RandomCheckGenerator::new()
            .with_n_bits(self.get_n_bits())
            .with_max_bit_degrees(self.bit_degree)
    }

    fn generate_layer_checks<R: Rng + ?Sized>(
        &self,
        layer: usize,
        check_generator: &mut RandomCheckGenerator,
        rng: &mut R,
    ) -> Vec<Vec<usize>> {
        let blocks = self.get_blocks_for_layer(layer);
        blocks
            .into_iter()
            .flat_map(|block| self.generate_checks_on_block(block, check_generator, rng))
            .collect()
    }

    fn get_blocks_for_layer(&self, layer: usize) -> Vec<Vec<usize>> {
        let n_blocks = self.get_n_blocks_for_layer(layer);
        let block_length = self.get_block_length_for_layer(layer);
        (0..n_blocks)
            .map(|block| {
                (0..block_length)
                    .map(|b| b + block * block_length)
                    .collect()
            })
            .collect()
    }

    fn get_n_blocks_for_layer(&self, layer: usize) -> usize {
        self.n_blocks_per_layer.pow((self.n_layers - layer) as u32)
    }

    fn get_block_length_for_layer(&self, layer: usize) -> usize {
        self.get_n_bits() / self.get_n_blocks_for_layer(layer) as usize
    }

    fn generate_checks_on_block<R: Rng + ?Sized>(
        &self,
        block: Vec<usize>,
        check_generator: &mut RandomCheckGenerator,
        rng: &mut R,
    ) -> Vec<Vec<usize>> {
        check_generator.over_bits(block);
        (0..self.n_checks_per_block)
            .filter_map(|_| check_generator.generate(self.check_degree as usize, rng))
            .collect()
    }

    // *********
    // Code size
    // *********

    pub fn get_n_bits(&self) -> usize {
        self.initial_block_length.pow(self.n_layers as u32)
    }

    pub fn get_n_checks(&self) -> usize {
        self.n_checks_per_block.pow(self.n_layers as u32) as usize
    }
}

impl CodeGenerator for HierarchicalCodeGenerator {
    fn generate<R: Rng + ?Sized>(&self, rng: &mut R) -> ParityCheckMatrix {
        if self.will_generate_an_empty_code() {
            ParityCheckMatrix::empty()
        } else {
            self.generate_non_empty_code(rng)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{thread_rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn default_constructor_generate_empty_code() {
        let generator = HierarchicalCodeGenerator::new();
        assert_eq!(
            generator.generate(&mut thread_rng()),
            ParityCheckMatrix::empty()
        );
    }

    #[test]
    fn generate_single_check_code() {
        let generator = HierarchicalCodeGeneratorBuilder::new()
            .with_bit_and_check_degree(3, 4)
            .build();

        let expected_code = ParityCheckMatrix::with_n_bits(4).with_checks(vec![vec![0, 1, 2, 3]]);
        assert_eq!(generator.generate(&mut thread_rng()), expected_code);
    }

    #[test]
    fn generate_single_layer_code() {
        let generator = HierarchicalCodeGeneratorBuilder::new()
            .with_bit_and_check_degree(3, 3)
            .with_initial_block_length(8)
            .with_n_checks_per_block(6)
            .build();

        let mut rng = ChaCha8Rng::seed_from_u64(10);
        let code = generator.generate(&mut rng);
        println!("{}", code);
    }
}

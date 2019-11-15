use super::{CodeGenerator, RandomCheckGenerator, BitAndCheckDegrees};
use crate::ParityCheckMatrix;
use rand::Rng;

pub struct HierarchicalCodeGeneratorBuilder {
    degrees: BitAndCheckDegrees,
    n_layers: Option<u32>,
    initial_block_length: Option<u32>,
    n_checks_per_block: Option<u32>,
    n_blocks_per_layer: Option<u32>,
}

impl HierarchicalCodeGeneratorBuilder {
    // ************
    // Construction
    // ************

    pub fn new(degrees: BitAndCheckDegrees) -> Self {
        Self {
            degrees,
            n_layers: None,
            initial_block_length: None,
            n_checks_per_block: None,
            n_blocks_per_layer: None,
        }        
    }

    pub fn with_n_layers(mut self, n_layers: u32) -> Self {
        self.n_layers = Some(n_layers);
        self
    }

    pub fn with_initial_block_length(mut self, length: u32) -> Self {
        self.initial_block_length = Some(length);
        self
    }

    pub fn with_n_checks_per_block(mut self, n_checks: u32) -> Self {
        self.n_checks_per_block = Some(n_checks);
        self
    }

    pub fn with_n_blocks_per_layer(mut self, n_blocks: u32) -> Self {
        self.n_blocks_per_layer = Some(n_blocks);
        self
    }

    // ********
    // Building
    // ********

    pub fn build(&self) -> HierarchicalCodeGenerator {
        HierarchicalCodeGenerator {
            degrees: self.degrees,
            n_layers: self.get_n_layers(),
            initial_block_length: self.get_initial_block_length(),
            n_checks_per_block: self.get_n_checks_per_block(),
            n_blocks_per_layer: self.get_n_blocks_for_layer(),
        }
    }

    fn get_n_layers(&self) -> u32 {
        self.n_layers.unwrap_or(1)
    }

    fn get_initial_block_length(&self) -> u32 {
        self.initial_block_length.unwrap_or(self.degrees.get_check_degrees())
    }

    fn get_n_checks_per_block(&self) -> u32 {
        self.n_checks_per_block.unwrap_or(1)
    }

    fn get_n_blocks_for_layer(&self) -> u32 {
        self.n_blocks_per_layer.unwrap_or(1)
    }
}


pub struct HierarchicalCodeGenerator {
    degrees: BitAndCheckDegrees,
    n_layers: u32,
    initial_block_length: u32,
    n_checks_per_block: u32,
    n_blocks_per_layer: u32,
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
            degrees: BitAndCheckDegrees::new(0, 0),
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
        self.degrees.will_generate_an_empty_code()
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
        let bit_degrees = vec![self.degrees.get_bit_degrees() as usize; self.get_n_bits()];
        RandomCheckGenerator::new(bit_degrees, 0)
    }

    fn generate_layer_checks<R: Rng + ?Sized>(
        &self, 
        layer: u32, 
        check_generator: &mut RandomCheckGenerator,
        rng: &mut R
    ) -> Vec<Vec<usize>> {
        let blocks = self.get_blocks_for_layer(layer);
        blocks
            .into_iter()
            .flat_map(|block| self.generate_checks_on_block(block, check_generator, rng))
            .collect()
    }

    fn get_blocks_for_layer(&self, layer: u32) -> Vec<Vec<usize>> {
        let n_blocks = self.get_n_blocks_for_layer(layer);
        let block_length = self.get_block_length_for_layer(layer);
        (0..n_blocks).map(|block| {
            (0..block_length).map(|b| b + block * block_length).collect()
        })
        .collect()
    }

    fn get_n_blocks_for_layer(&self, layer: u32) -> usize
     {
        self.n_blocks_per_layer.pow(self.n_layers - layer) as usize
    }

    fn get_block_length_for_layer(&self, layer: u32) -> usize {
        self.get_n_bits() / self.get_n_blocks_for_layer(layer) as usize
    }

    fn generate_checks_on_block<R: Rng + ?Sized>(
        &self,
        block: Vec<usize>, 
        check_generator: &mut RandomCheckGenerator,
        rng: &mut R
    ) -> Vec<Vec<usize>> {
        check_generator.over_bits(block);
        (0..self.n_checks_per_block)
            .filter_map(|_| {
                check_generator.generate(self.degrees.get_check_degrees() as usize, rng)
            })
            .collect()
    }

    // *********
    // Code size
    // *********

    pub fn get_n_bits(&self) -> usize {
        self.initial_block_length.pow(self.n_layers) as usize
    }

    pub fn get_n_checks(&self) -> usize {
        self.n_checks_per_block.pow(self.n_layers) as usize
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
    use crate::BitAndCheckDegrees;
    use rand::{SeedableRng, thread_rng};
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn default_constructor_generate_empty_code() {
        let generator = HierarchicalCodeGenerator::new();
        assert_eq!(generator.generate(&mut thread_rng()), ParityCheckMatrix::empty());
    }

    #[test]
    fn generate_single_check_code() {
        let degrees = BitAndCheckDegrees::new(3, 4);
        let generator = HierarchicalCodeGeneratorBuilder::new(degrees).build();

        let expected_code = ParityCheckMatrix::with_n_bits(4).with_checks(vec![vec![0, 1, 2, 3]]);
        assert_eq!(generator.generate(&mut thread_rng()), expected_code);
    }

    #[test]
    fn generate_single_layer_code() {
        let degrees = BitAndCheckDegrees::new(3, 3);
        let generator = HierarchicalCodeGeneratorBuilder::new(degrees)
            .with_initial_block_length(8)
            .with_n_checks_per_block(6)
            .build();

        let mut rng = ChaCha8Rng::seed_from_u64(10);    
        let code = generator.generate(&mut rng);
        println!("{}", code);

    }
}

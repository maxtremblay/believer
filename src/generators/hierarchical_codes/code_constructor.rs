use crate::random_checks::Generator as CheckGenerator;
use crate::ParityCheckMatrix;
use rand::Rng;

pub(super) struct CodeConstructor<R: Rng> {
    pub(super) check_generator: CheckGenerator<R>,
    pub(super) checks: Vec<Vec<usize>>,
    pub(super) active_layer: u32,
    pub(super) n_layers: u32,
    pub(super) n_bits: usize,
    pub(super) n_blocks_per_layer: usize,
    pub(super) n_checks_per_block: usize,
}

impl<R: Rng> CodeConstructor<R> {
    pub(super) fn get_random_code(mut self) -> ParityCheckMatrix {
        (0..self.n_layers).for_each(|_| {
            self.generate_active_layer();
            self.active_layer += 1;
        });
        self.get_parity_check_matrix()
    }

    fn generate_active_layer(&mut self) {
        let blocks = self.get_blocks_for_active_layer();
        blocks
            .into_iter()
            .for_each(|block| self.generate_checks_over(block));
    }

    fn get_blocks_for_active_layer(&self) -> Blocks {
        let n_blocks = self.get_n_blocks_for_active_layer();
        let block_length = self.get_block_length_for_active_layer();
        Blocks {
            n_blocks,
            block_length,
            active_block: 0,
        }
    }

    fn get_n_blocks_for_active_layer(&self) -> usize {
        self.n_blocks_per_layer
            .pow(self.n_layers - self.active_layer)
    }

    fn get_block_length_for_active_layer(&self) -> usize {
        self.n_bits / self.get_n_blocks_for_active_layer()
    }

    fn generate_checks_over(&mut self, block: Vec<usize>) {
        self.check_generator.set_over_bits(block);
        (0..self.n_checks_per_block).for_each(|_| self.add_check_to_code());
    }

    fn add_check_to_code(&mut self) {
        if let Some(check) = self.check_generator.get_random_check() {
            self.checks.push(check);
        }
    }

    fn get_parity_check_matrix(self) -> ParityCheckMatrix {
        ParityCheckMatrix::with_n_bits(self.n_bits).with_checks(self.checks)
    }
}

struct Blocks {
    block_length: usize,
    n_blocks: usize,
    active_block: usize,
}

impl Blocks {
    fn get_bits_for_active_block(&self) -> Vec<usize> {
        (0..self.block_length)
            .map(|bit| bit + self.active_block * self.block_length)
            .collect()
    }
}

impl Iterator for Blocks {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.active_block < self.n_blocks {
            let bits = self.get_bits_for_active_block();
            self.active_block += 1;
            Some(bits)
        } else {
            None
        }
    }
}

// fn generate_layer_checks<R: Rng>(
//     &self,
//     layer: usize,
//     check_generator: &mut CheckGenerator<R>
// ) -> Vec<Vec<usize>> {
//     let blocks = self.get_blocks_for_layer(layer);
//     blocks
//         .into_iter()
//         .flat_map(|block| self.generate_checks_on_block(block, check_generator, rng))
//         .collect()
// }

// fn get_blocks_for_layer(&self, layer: usize) -> Vec<Vec<usize>> {
//     let n_blocks = self.get_n_blocks_for_layer(layer);
//     let block_length = self.get_block_length_for_layer(layer);
//     (0..n_blocks)
//         .map(|block| {
//             (0..block_length)
//                 .map(|b| b + block * block_length)
//                 .collect()
//         })
//         .collect()
// }

// fn get_n_blocks_for_layer(&self, layer: usize) -> usize {
//     self.n_blocks_per_layer.pow((self.n_layers - layer) as u32)
// }

// fn get_block_length_for_layer(&self, layer: usize) -> usize {
//     self.get_n_bits() / self.get_n_blocks_for_layer(layer)
// }

// fn generate_checks_on_block<R: Rng>(
//     &self,
//     block: Vec<usize>,
//     check_generator: &mut RandomCheckGenerator,
//     rng: &mut R,
// ) -> Vec<Vec<usize>> {
//     check_generator.set_over_bits(block);
//     (0..self.n_checks_per_block)
//         .filter_map(|_| check_generator.get_random_check())
//         .collect()
// }

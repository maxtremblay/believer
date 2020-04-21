use super::{Decoder, DecodingResult, SimulationResult};
use rand::distributions::Standard;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

pub(super) struct NIterationsSimulator<'a, D: Decoder> {
    decoder: &'a mut D,
    n_iterations: usize,
    n_successes: usize,
    random_seeds: Vec<u64>,
}

impl<'a, D: Decoder> NIterationsSimulator<'a, D> {
    pub(super) fn from(decoder: &'a mut D) -> Self {
        Self {
            decoder,
            n_iterations: 0,
            n_successes: 0,
            random_seeds: Vec::new(),
        }
    }

    pub(super) fn simulate_n_iterations_with_rng<R: Rng>(
        mut self,
        n_iterations: usize,
        rng: &mut R,
    ) -> Self {
        self.initialize_simulation_with_n_iterations_and_rng(n_iterations, rng);
        self.run_the_simulation();
        self
    }

    fn initialize_simulation_with_n_iterations_and_rng<R: Rng>(
        &mut self,
        n_iterations: usize,
        rng: &mut R,
    ) {
        self.n_iterations = n_iterations;
        self.initialize_random_seeds_with_rng(rng);
    }

    fn initialize_random_seeds_with_rng<R: Rng>(&mut self, rng: &mut R) {
        self.random_seeds = rng.sample_iter(Standard).take(self.n_iterations).collect()
    }

    fn run_the_simulation(&mut self) {
        self.n_successes = (0..self.n_iterations)
            // .into_par_iter()
            .filter(|thread_index| {
                let mut rng = self.get_thread_rng(*thread_index);
                self.decoder
                    .decode_random_error_with_rng(&mut rng)
                    .is_success()
            })
            .count();
    }

    // Yep, I'm imposing ChaCha8Rng for each thread.
    // I don't have a better solution for now that preserve reproductability.
    fn get_thread_rng(&self, thread_index: usize) -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(self.random_seeds[thread_index])
    }

    pub(super) fn get_result(&self) -> SimulationResult {
        let n_failures = self.n_iterations - self.n_successes;
        SimulationResult::with_n_successes_and_failures(self.n_successes as u64, n_failures as u64)
    }
}

// #[cfg(test)]
// mod test {
//     use super::super::ErasureDecoder;
//     use super::*;
//     use crate::ParityCheckMatrix;
//     use rand::SeedableRng;
//     use rand_chacha::ChaCha8Rng;

//     #[test]
//     fn there_is_n_iterations() {
//         let code = ParityCheckMatrix::with_n_bits(3).with_checks(vec![vec![0, 1], vec![1, 2]]);

//         let decoder = ErasureDecoder::with_prob(0.5).for_code(code);
//         let rng = ChaCha8Rng::seed_from_u64(123);
//         let n_iterations = NIterationsSimulator::from(&decoder)
//             .simulate_n_iterations_with_rng(1000, &mut rng.clone())
//             .get_result()
//             .get_n_iterations();
//         assert_eq!(n_iterations, 1000);
//     }

//     #[test]
//     fn reproductibility_for_repetition_code() {
//         let code = ParityCheckMatrix::with_n_bits(3).with_checks(vec![vec![0, 1], vec![1, 2]]);

//         let decoder = ErasureDecoder::with_prob(0.5).for_code(code);
//         let rng = ChaCha8Rng::seed_from_u64(123);
//         let result_0 = NIterationsSimulator::from(&decoder)
//             .simulate_n_iterations_with_rng(1000, &mut rng.clone())
//             .get_result()
//             .get_success_rate();

//         let result_1 = NIterationsSimulator::from(&decoder)
//             .simulate_n_iterations_with_rng(1000, &mut rng.clone())
//             .get_result()
//             .get_success_rate();

//         assert!((result_0 - result_1).abs() < 1e-6);
//     }

//     #[test]
//     fn reproductibility_for_hamming_code() {
//         let code = ParityCheckMatrix::with_n_bits(7).with_checks(vec![
//             vec![0, 1, 2, 4],
//             vec![0, 1, 3, 5],
//             vec![0, 2, 3, 6],
//         ]);

//         let decoder = ErasureDecoder::with_prob(0.5).for_code(code);
//         let rng = ChaCha8Rng::seed_from_u64(123);
//         let result_0 = NIterationsSimulator::from(&decoder)
//             .simulate_n_iterations_with_rng(1000, &mut rng.clone())
//             .get_result()
//             .get_success_rate();

//         let result_1 = NIterationsSimulator::from(&decoder)
//             .simulate_n_iterations_with_rng(1000, &mut rng.clone())
//             .get_result()
//             .get_success_rate();

//         assert!((result_0 - result_1).abs() < 1e-6);
//     }
// }

use super::{Decoder, SimulationResult};
use rand::distributions::Standard;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

pub(super) struct NEventsSimulator<'a, D> {
    decoder: &'a mut D,
    n_events: usize,
    result: SimulationResult,
    random_seeds: Vec<u64>,
}

impl<'a, D: Decoder> NEventsSimulator<'a, D> {
    pub(super) fn from(decoder: &'a mut D) -> Self {
        Self {
            decoder,
            n_events: 0,
            result: SimulationResult::new(),
            random_seeds: Vec::new(),
        }
    }

    pub(super) fn simulate_until_n_events_are_found_with_rng<R: Rng>(
        mut self,
        n_events: usize,
        rng: &mut R,
    ) -> Self {
        self.initialize_simulation_with_n_events_and_rng(n_events, rng);
        self.run_the_simulation();
        self
    }

    fn initialize_simulation_with_n_events_and_rng<R: Rng>(
        &mut self,
        n_events: usize,
        rng: &mut R,
    ) {
        self.n_events = n_events;
        self.initialize_random_seeds_with_rng(rng);
    }

    fn run_the_simulation(&mut self) {
        let results = (0..self.n_events)
            // .into_par_iter()
            .map(|thread_index| self.simulate_thread_until_one_event_is_found(thread_index));

        let mut n_successes: u64 = 0;
        let mut n_failures: u64 = 0;

        for simres in results {
            n_successes += simres.get_n_successes();
            n_failures += simres.get_n_failures();
        }

        self.result = SimulationResult::with_n_successes_and_failures(n_successes, n_failures);
    }

    fn simulate_thread_until_one_event_is_found(&mut self, thread_index: usize) -> SimulationResult {
        let mut rng = self.get_thread_rng(thread_index);
        let mut result = SimulationResult::new();
        while result.has_not_at_least_one_success_and_one_failure() {
            let decoding_result = self.decoder.decode_random_error_with_rng(&mut rng);
            result.add_decoding_result(decoding_result);
        }
        result
    }

    fn initialize_random_seeds_with_rng<R: Rng>(&mut self, rng: &mut R) {
        self.random_seeds = rng.sample_iter(Standard).take(self.n_events).collect()
    }

    // Yep, I'm imposing ChaCha8Rng with different seeds for each thread.
    // I don't have a better solution for now that preserve reproductability.
    fn get_thread_rng(&self, thread_index: usize) -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(self.random_seeds[thread_index])
    }

    pub(super) fn get_result(self) -> SimulationResult {
        self.result
    }
}

// #[cfg(test)]
// mod test {
//     use super::*;
//     use super::super::ErasureDecoder;
//     use crate::ParityCheckMatrix;
//     use rand_chacha::ChaCha8Rng;
//     use rand::SeedableRng;

//     #[test]
//     fn 
//     there_are_at_least_n_failures_and_n_erasures() {
//         let code = ParityCheckMatrix::with_n_bits(3)
//         .with_checks(vec![vec![0, 1], vec![1, 2]]);

//         let decoder = ErasureDecoder::with_prob(0.5).for_code(code);
//         let rng = ChaCha8Rng::seed_from_u64(123);
        
//         let result = NEventsSimulator::from(&decoder)
//             .simulate_until_n_events_are_found_with_rng(10, &mut rng.clone())
//             .get_result();

//         assert!(result.get_n_failures() >= 10);
//         assert!(result.get_n_successes() >= 10);

//     }

//     #[test]
//     fn reproductibility_for_repetition_code() {
//         let code = ParityCheckMatrix::with_n_bits(3)
//             .with_checks(vec![vec![0, 1], vec![1, 2]]);

//         let decoder = ErasureDecoder::with_prob(0.5).for_code(code);
//         let rng = ChaCha8Rng::seed_from_u64(123);
        
//         let result_0 = NEventsSimulator::from(&decoder)
//             .simulate_until_n_events_are_found_with_rng(10, &mut rng.clone())
//             .get_result()
//             .get_success_rate();

//         let result_1 = NEventsSimulator::from(&decoder)
//             .simulate_until_n_events_are_found_with_rng(10, &mut rng.clone())
//             .get_result()
//             .get_success_rate();

//         assert!((result_0 - result_1).abs() < 1e-6); 
//     }

//     #[test]
//     fn reproductibility_for_hamming_code() {
//         let code = ParityCheckMatrix::with_n_bits(7)
//             .with_checks(vec![vec![0, 1, 2, 4], vec![0, 1, 3, 5], vec![0, 2, 3, 6]]);

//         let decoder = ErasureDecoder::with_prob(0.5).for_code(code);
//         let rng = ChaCha8Rng::seed_from_u64(123);
        
//         let result_0 = NEventsSimulator::from(&decoder)
//             .simulate_until_n_events_are_found_with_rng(10, &mut rng.clone())
//             .get_result()
//             .get_success_rate();

//         let result_1 = NEventsSimulator::from(&decoder)
//             .simulate_until_n_events_are_found_with_rng(10, &mut rng.clone())
//             .get_result()
//             .get_success_rate();

//         assert!((result_0 - result_1).abs() < 1e-6); 
//     }
// }
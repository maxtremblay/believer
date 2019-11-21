use super::{Decoder, DecodingResult, SimulationResult};
use rayon::prelude::*;

pub(super) struct NIterationsSimulator<'a, D: Decoder> {
    decoder: &'a D,
    n_iterations: u64,
    n_successes: u64,
}

impl<'a, D: Decoder> NIterationsSimulator<'a, D> {
    pub(super) fn from(decoder: &'a D) -> Self {
        Self { decoder, n_iterations: 0, n_successes: 0 }
    }
    
    pub(super) fn simulate_n_iterations(mut self, n_iterations: u64) -> Self {
        self.n_iterations = n_iterations;
        self.n_successes = (0..n_iterations)
            .into_par_iter()
            .filter(|_| self.decoder.decode_random_error().is_success())
            .count() as u64;
        self
    }

    pub(super) fn get_result(&self) -> SimulationResult {
        let n_failures = self.n_iterations - self.n_successes;
        SimulationResult::with_n_successes_and_failures(self.n_successes, n_failures)
    }
}

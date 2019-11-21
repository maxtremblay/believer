use super::{Decoder, DecodingResult, SimulationResult};
use rayon::prelude::*;

pub(super) struct NEventsSimulator<'a, D> {
    decoder: &'a D,
    n_events: u64,
    result: SimulationResult,
}

impl<'a, D: Decoder> NEventsSimulator<'a, D> {
    pub(super) fn from(decoder: &'a D) -> Self {
        Self { decoder, n_events: 0, result: SimulationResult::new() }
    }
    
    pub(super) fn simulate_until_n_events_are_found(mut self, n_events: u64) -> Self {
        self.n_events = n_events;
        self.result = (0..n_events)
            .into_par_iter()
            .map(|_| self.simulate_until_one_event_is_found())
            .reduce(
                || SimulationResult::new(), 
                |result_accumulator, result| result_accumulator.combine_with(result)
            );
        self
    }

    fn simulate_until_one_event_is_found(&self) -> SimulationResult {
        let mut n_successes = 0;
        let mut n_failures = 0;
        while n_successes == 0 || n_failures == 0 {
            if self.decoder.decode_random_error().is_success() {
                n_successes += 1;
            } else {
                n_failures += 1;
            }
        }
        SimulationResult::with_n_successes_and_failures(n_successes, n_failures)
    }

    pub(super) fn get_result(self) -> SimulationResult {
        self.result
    }
}

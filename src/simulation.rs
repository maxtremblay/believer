//! Toolbox to simulate error correction codes performance.

use crate::{Decoder, DecodingResult};
use rayon::prelude::*;

pub struct Simulator<'a, D>
where
    D: Decoder,
{
    decoder: &'a D,
}

impl<'a, D> Simulator<'a, D>
where
    D: Decoder,
{
    pub fn new(decoder: &'a D) -> Self {
        Self { decoder }
    }

    pub fn simulate_n_iterations(&self, n_iterations: usize) -> SimulationResult {
        let successes: usize = (0..n_iterations)
            .into_par_iter()
            .map::<_, D::Result>(|_| {
                let error = self.decoder.random_error();
                self.decoder.decode(&error)
            })
            .filter(|decoding| decoding.succeed())
            .count();

        SimulationResult {
            successes,
            failures: n_iterations - successes,
        }
    }

    pub fn simulate_until_failures_are_found(
        &self,
        n_threads: usize,
        n_failures_per_thread: usize,
    ) -> SimulationResult {
        let successes = (0..n_threads)
            .into_par_iter()
            .map::<_, usize>(|_| {
                (0..n_failures_per_thread)
                    .into_par_iter()
                    .map(|_| {
                        let mut successes = 0;
                        let mut has_failed = false;
                        while !has_failed {
                            let error = self.decoder.random_error();
                            if self.decoder.decode(&error).succeed() {
                                successes += 1;
                            } else {
                                has_failed = true;
                            }
                        }
                        successes
                    })
                    .sum()
            })
            .sum();

        SimulationResult {
            successes,
            failures: n_failures_per_thread * n_threads,
        }
    }
}

pub struct SimulationResult {
    successes: usize,
    failures: usize,
}

impl SimulationResult {
    pub fn failure_rate(&self) -> f64 {
        self.failures() as f64 / self.total() as f64
    }

    pub fn failures(&self) -> usize {
        self.failures
    }

    pub fn success_rate(&self) -> f64 {
        self.successes() as f64 / self.total() as f64
    }

    pub fn successes(&self) -> usize {
        self.successes
    }

    pub fn total(&self) -> usize {
        self.failures() + self.successes()
    }
}

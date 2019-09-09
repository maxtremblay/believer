//! Toolbox to simulate error correction codes performance.

use crate::{BPDecoder, BPResult, BinaryChannel, ErasureDecoder, ErasureResult, ParityCheckMatrix};
use rayon::prelude::*;

pub trait Simulator {
    fn simulate(&self) -> SimulationResult;
}

// pub struct BinaryChannelSimulator<'a, C>
// where
//     C: BinaryChannel,
// {
//     decoder: BPDecoder<'a, C>,

// }

pub struct ErasureSimulator<'a> {
    decoder: ErasureDecoder<'a>,
    n_threads: usize,
    n_failures_per_thread: usize,
    erasure_prob: f64,
}

impl<'a> ErasureSimulator<'a> {
    pub fn new(
        checks: &'a ParityCheckMatrix,
        n_threads: usize,
        n_failures_per_thread: usize,
        erasure_prob: f64,
    ) -> Self {
        if erasure_prob < 0.0 || erasure_prob > 1.0 {
            panic!("invalid probability");
        }
        Self {
            decoder: ErasureDecoder::new(checks),
            n_threads,
            n_failures_per_thread,
            erasure_prob,
        }
    }
}

impl<'a> Simulator for ErasureSimulator<'a> {
    fn simulate(&self) -> SimulationResult {
        let successes = (0..self.n_threads)
            .into_par_iter()
            .map::<_, usize>(|_| {
                (0..self.n_failures_per_thread)
                    .into_par_iter()
                    .map(|_| {
                        let mut successes = 0;
                        let mut has_failed = false;
                        while !has_failed {
                            let erasure_locations = self.decoder.random_erasures(self.erasure_prob);
                            if self.decoder.decode(&erasure_locations) == ErasureResult::Success {
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
            failures: self.n_failures_per_thread * self.n_threads,
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

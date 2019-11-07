//! Toolbox to simulate error correction codes performance.

use crate::{Decoder, DecodingResult};
use rayon::prelude::*;

pub struct Simulator<D: Decoder> {
    decoder: D
}

impl<D: Decoder> Simulator<D> {
    fn new(decoder: D) -> Self {
        Self { decoder }
    }

    fn simulate_n_iterations(&self, n_iterations: usize) -> SimulationResult {
        let decoder = self.decoder();
        let successes: usize = (0..n_iterations)
            .map(|_| {
                let error = decoder.random_error();
                decoder.decode(&error)
            })
            .filter(|decoding| decoding.succeed())
            .count();

        SimulationResult {
            successes,
            failures: n_iterations - successes,
        }
    }

    fn simulate_until_events_are_found(&self, n_events: usize) -> SimulationResult {
        let (succ, fail) = (0..n_events)
            .into_par_iter()
            .map(|_| {
                let mut successes = 0;
                let mut failures = 0;
                while successes == 0 || failures == 0 {
                    let error = self.decoder.random_error();
                    if self.decoder.decode(&error).succeed() {
                        successes += 1;
                    } else {
                        failures += 1;
                    }
                }
                (successes, failures)
            })
            .reduce(|| (0, 0), |(a_s, a_f), (s, f)| (a_s + s, a_f + f));

        SimulationResult {
            successes: succ,
            failures: fail,
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

    pub fn worse_result() -> Self {
        Self {
            successes: 0,
            failures: 1,
        }
    }
}

// pub struct CSSSimulationResult {
//     successes: usize,
//     x_failures: usize,
//     z_failures: usize,
//     x_and_z_failures: usize,
// }

// impl CSSSimulationResult {
//     pub fn x_failures(&self) -> usize {
//         self.x_failures
//     }

//     pub fn x_failure_rate(&self) -> f64 {
//         self.x_failures() as f64 / self.total() as f64
//     }

//     pub fn z_failures(&self) -> usize {
//         self.z_failures
//     }

//     pub fn z_failure_rate(&self) -> f64 {
//         self.z_failures() as f64 / self.total() as f64
//     }

//     pub fn x_and_z_failures(&self) -> usize {
//         self.x_and_z_failures
//     }

//     pub fn x_and_z_failure_rate(&self) -> f64 {
//         self.x_and_z_failures() as f64 / self.total() as f64
//     }

//     pub fn total_failures(&self) -> usize {
//         self.x_failures + self.z_failures + self.x_and_z_failures
//     }

//     pub fn total_failure_rate(&self) -> f64 {
//         self.total_failures() as f64 / self.total() as f64
//     }

//     pub fn successes(&self) -> usize {
//         self.successes
//     }

//     pub fn total(&self) -> usize {
//         self.total_failures() + self.successes()
//     }
// }

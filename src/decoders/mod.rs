//! Toolbox for decoding.

use rand::{Rng, thread_rng};
use rayon::prelude::*;

pub mod simulation_results;
pub use simulation_results::SimulationResult;

mod n_iterations_simulator;
use n_iterations_simulator::NIterationsSimulator;

/// An interface to deal with decoders
///
/// This is the global decoder trait. For more details, see each decoder implementation.
pub trait Decoder: Send + Sync + Sized {
    /// The type of code the decoder is using.
    type Code;

    /// The type of error that the decoder can decode.
    type Error;

    /// The type of result the decoder is returning.
    type Result: DecodingResult;

    /// Set self to use `code` without changing the other parameter. This consume `code`
    fn set_code(&mut self, code: Self::Code);

    /// Takes the `code` out of the decoder leaving an empty set of code instead.
    fn take_code(&mut self) -> Self::Code;

    /// Tries to decode a given error.
    fn decode(&self, error: &Self::Error) -> Self::Result;

    /// Generates a random error.
    fn get_random_error(&self) -> Self::Error;

    /// Generates and decodes a random error.
    fn decode_random_error(&self) -> Self::Result {
        self.decode(&self.get_random_error())
    }

    /// Simulates the decoder for `n_iterations`. 
    fn simulate_n_iterations(&self, n_iterations: u64) -> SimulationResult {
        NIterationsSimulator::from(self)
            .simulate_n_iterations(n_iterations)
            .get_result()
    }

    // /// 
    // fn simulate_until_n_events_are_found(&self, n_events: usize) -> SimulationResult {
    //     NEventsSimulator::from(self)
    //         .until_n_events_are_found(n_events)
    //         .get_simulation_result()

    //     // let (succ, fail) = (0..n_events)
    //     //     .into_par_iter()
    //     //     .map(|_| {
    //     //         let mut successes = 0;
    //     //         let mut failures = 0;
    //     //         while successes == 0 || failures == 0 {
    //     //             let error = self.get_random_error();
    //     //             if self.decode(&error).is_success() {
    //     //                 successes += 1;
    //     //             } else {
    //     //                 failures += 1;
    //     //             }
    //     //         }
    //     //         (successes, failures)
    //     //     })
    //     //     .reduce(|| (0, 0), |(a_s, a_f), (s, f)| (a_s + s, a_f + f));

    //     // SimulationResult::with_n_successes_and_failures(succ, fail)
    // }
}

/// An interface for decoder outcome.
///
/// Decoding can either succeed or fail. However, it is possible that there are many kind of
/// success and failures.
pub trait DecodingResult: Send + Sync {
    /// Returns `true` if the decoding procedure succeed, `false` otherwise.
    fn is_success(&self) -> bool;

    /// Returns `false` if the decoding procedure succeed, `true` otherwise.
    fn is_failure(&self) -> bool {
        !self.is_success()
    }
}

// pub mod belief_propagation;
// pub use belief_propagation::*;

// pub mod erasure;
// pub use erasure::*;

// pub mod quantum_erasure;
// pub use quantum_erasure::{QuantumErasureDecoder, QuantumErasureDecoderBuilder};

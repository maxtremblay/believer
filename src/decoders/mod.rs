//! Toolbox for decoding.

use rand::{Rng, thread_rng};

pub mod simulation_results;
pub use simulation_results::SimulationResult;

mod n_iterations_simulator;
use n_iterations_simulator::NIterationsSimulator;

mod n_events_simulator;
use n_events_simulator::NEventsSimulator;

// pub mod belief_propagation;
// pub use belief_propagation::*;

pub mod erasure;
pub use erasure::*;

// pub mod quantum_erasure;
// pub use quantum_erasure::{QuantumErasureDecoder, QuantumErasureDecoderBuilder};


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

    /// Simulates decoding random error using `self` for `n_iterations`. 
    fn simulate_n_iterations_with_rng<R: Rng>(
        &self, 
        n_iterations: u64,
        rng: &mut R
    ) -> SimulationResult {
        NIterationsSimulator::from(self)
            .simulate_n_iterations(n_iterations)
            .get_result()
    }

    fn simulate_n_iterations(&self, n_iterations: u64) -> SimulationResult {
        self.simulate_n_iterations_with_rng(n_iterations, &mut thread_rng())
    }

    /// Simulates the decoder until `n_events` are found. 
    /// 
    /// That is, simulate until `n_events` successes and `n_events` are found.
    fn simulate_until_n_events_are_found(&self, n_events: u64) -> SimulationResult {
        NEventsSimulator::from(self)
            .simulate_until_n_events_are_found(n_events)
            .get_result()
    }
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

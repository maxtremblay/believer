//! Toolbox for decoding.

use rand::{Rng, thread_rng};

/// An interface to deal with decoders
///
/// This is the global decoder trait. For more details, see each decoder implementation.
pub trait Decoder: Send + Sync {
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

    fn get_random_error_with_rng<R: Rng>(&self, rng: &mut R) -> Self::Error;

    /// Generates a random error.
    fn get_random_error(&self) -> Self::Error {
        self.get_random_error_with_rng(&mut thread_rng())
    }  

}

/// An interface for decoder outcome.
///
/// Decoding can either succeed or fail. However, it is possible that there are many kind of
/// success and failures.
pub trait DecodingResult: Send + Sync {
    /// Returns `true` if the decoding procedure succeed, `false` otherwise.
    fn succeed(&self) -> bool;

    /// Returns `false` if the decoding procedure succeed, `true` otherwise.
    fn failed(&self) -> bool {
        !self.succeed()
    }
}

// pub mod belief_propagation;
// pub use belief_propagation::*;

pub mod erasure;
pub use erasure::*;

// pub mod quantum_erasure;
// pub use quantum_erasure::{QuantumErasureDecoder, QuantumErasureDecoderBuilder};

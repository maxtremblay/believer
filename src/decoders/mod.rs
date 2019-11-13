// use crate::ParityCheckMatrix;

/// An interface to deal with decoders
///
/// This is the global decoder trait. For more details, see each decoder implementation.
pub trait Decoder: Send + Sync {
    /// The type of error that the decoder can decode.
    type Error;

    /// The type of result the decoder is returning.
    type Result: DecodingResult;

    /// The type of checks the decoder is using.
    type Checks;

    /// Tries to decode a given error.
    fn decode(&self, error: &Self::Error) -> Self::Result;

    /// Generates a random error.
    fn random_error(&self) -> Self::Error;

    /// Takes the `checks` out of the decoder leaving an empty set of checks instead.
    fn take_checks(&mut self) -> Self::Checks;
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

/// A tool to build decoder.
///
/// It is often useful to construct decoders with the same parameters by varying only the checks.
/// A `DecoderBuilder` can be set from a set of parameters and generate a decoder when some checks
/// are specified. See different implementations for more details.
pub trait DecoderBuilder {
    /// The type of checks the decoder is using.
    type Checks;

    /// The type of decoder to build.
    type Decoder;

    /// Returns a `Decoder` build from the given `checks`.
    fn build_from(&self, checks: Self::Checks) -> Self::Decoder;
}

pub mod belief_propagation;
pub use belief_propagation::*;

pub mod erasure;
pub use erasure::*;

// pub mod quantum_erasure;
// pub use quantum_erasure::{QuantumErasureDecoder, QuantumErasureDecoderBuilder};

// use crate::ParityCheckMatrix;

/// TO COMMENT
pub trait Decoder: Send + Sync {
    type Error;
    type Result: DecodingResult;
    type Checks;

    /// TO COMMENT
    fn decode(&self, error: &Self::Error) -> Self::Result;

    /// TO COMMENT
    fn random_error(&self) -> Self::Error;

    /// TO COMMENT
    fn take_checks(&mut self) -> Self::Checks;
}

/// TO COMMENT
pub trait DecodingResult: Send + Sync {
    /// TO COMMENT
    fn succeed(&self) -> bool;

    /// TO COMMENT
    fn failed(&self) -> bool {
        !self.succeed()
    }
}

/// TO COMMENT
pub trait DecoderBuilder {
    /// TO COMMENT
    type Code;

    /// TO COMMENT
    type Decoder;

    /// TO COMMENT
    fn from_code(&self, code: Self::Code) -> Self::Decoder;
}

pub mod belief_propagation;
pub use belief_propagation::{BPDecoder, BPResult};

pub mod erasure;
pub use erasure::{ErasureDecoder, ErasureDecoderBuilder, ErasureResult};

pub mod css_erasure;
pub use css_erasure::{CSSErasureDecoder, CSSErasureDecoderBuilder};

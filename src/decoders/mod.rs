pub trait Decoder: Send + Sync {
    type Error;
    type Result: DecodingResult;
    fn decode(&self, error: &Self::Error) -> Self::Result;
    fn random_error(&self) -> Self::Error;
}

pub trait DecodingResult: Send + Sync {
    fn succeed(&self) -> bool;

    fn failed(&self) -> bool {
        !self.succeed()
    }
}

pub mod belief_propagation;
pub use belief_propagation::{BPDecoder, BPResult};

pub mod erasure;
pub use erasure::{ErasureDecoder, ErasureResult};

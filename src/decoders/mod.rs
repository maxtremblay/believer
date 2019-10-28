
use crate::ParityCheckMatrix;

pub trait Checks {
}

impl Checks for ParityCheckMatrix {}

impl Checks for (ParityCheckMatrix,ParityCheckMatrix) {}

pub trait Code<Cks: Checks> {
    fn get_checks(&self) -> Cks;
}


pub trait Decoder: Send + Sync {
    type Error;
    type Result: DecodingResult;
    fn decode(&self, error: &Self::Error) -> Self::Result;
    fn random_error(&self) -> Self::Error;
}
//    fn new(checks: &I, erasure_prob: f64);

pub trait DecodingResult: Send + Sync {
    fn succeed(&self) -> bool;

    fn failed(&self) -> bool {
        !self.succeed()
    }
}

pub trait DecoderBuilder<C, D> where
C: Checks,
D: Decoder{

    fn from_code(&self,code: C) -> D;

}

pub mod belief_propagation;
pub use belief_propagation::{BPDecoder, BPResult};

pub mod erasure;
pub use erasure::{ErasureDecoder, ErasureResult};

pub mod css_erasure;
pub use css_erasure::{CSSErasureDecoder, CSSErasureResult, CSSEDBuilder};
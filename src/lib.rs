//! A belief propapagation decoder for classical and quantum sparse error correcting codes.

pub mod channel;
pub use channel::BinaryChannel;

pub mod decoders;
pub use decoders::{BPDecoder, BPResult, ErasureDecoder, ErasureResult};

pub mod gf2;
pub use gf2::GF2;

pub mod parity_check_matrix;
pub use parity_check_matrix::ParityCheckMatrix;

pub mod simulation;

mod sparse_matrix;

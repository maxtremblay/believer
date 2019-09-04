//! A belief propapagation decoder for classical and quantum sparse error correcting codes.

pub mod channel;
pub use channel::BinaryChannel;

pub mod decoder;
pub use decoder::{Decoder, DecodingResult};

pub mod gf2;
pub use gf2::GF2;

pub mod parity_check_matrix;
pub use parity_check_matrix::ParityCheckMatrix;

mod sparse_matrix;

pub mod random_codes;
pub use random_codes::CodeGenerator;

//! A belief propapagation decoder for classical and quantum sparse error correcting codes.

pub mod channel;
pub use channel::*;

pub mod decoders;
pub use decoders::*;

pub mod gf2;
pub use gf2::*;

// pub mod gf4_stabilizers;
// pub use gf4_stabilizers::*;

// pub mod generators;
// pub use generators::*;

pub mod parity_check_matrix;
pub use parity_check_matrix::*;

// pub mod paulis;
// pub use paulis::*;

mod sparse_matrix;

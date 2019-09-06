//! A classical erasure decoder.

use crate::ParityCheckMatrix;

pub struct ErasureDecoder<'a> {
    parity_check: &'a ParityCheckMatrix,
}

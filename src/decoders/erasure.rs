//! A classical erasure decoder.

use super::{Decoder, DecodingResult};
use crate::ParityCheckMatrix;
use rand::Rng;

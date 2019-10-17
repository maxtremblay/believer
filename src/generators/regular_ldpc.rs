use crate::ParityCheckMatrix;
use super::CodeGenerator;
use rayon::prelude::*;

pub struct RegularLdpcGenerator {
    bit_degree: usize,
    check_degree: usize,
    scale: usize,
    min_girth: usize,
}

impl RegularLdpcGenerator {
    // *************
    // Public method
    // *************

    /// Creates a new generator that will generate (`bit_degree`, `check_degree`)-regular LDPC
    /// codes with `scale * check_degree` bits and `scale * bit_degree` checks. The codes are
    /// garanteed to have girth at least `min_girth`.
    pub fn new(bit_degree: usize, check_degree: usize, scale: usize, min_girth: usize) -> Self {
        Self { bit_degree, check_degree, scale, min_girth }
    }
}

impl CodeGenerator for RegularLdpcGenerator {
    fn generate(&self) -> ParityCheckMatrix {
        unimplemented!()
    }
}
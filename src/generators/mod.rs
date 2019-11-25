use crate::ParityCheckMatrix;
use rand::{thread_rng, Rng};

pub mod best_code_finder;
pub use best_code_finder::BestCodeFinderFromErasure;

pub mod random_checks;

pub mod regular_ldpc;
pub use regular_ldpc::*;

pub mod hierarchical_codes;

pub mod increasing_range_code;
pub use increasing_range_code::*;

pub trait CodeGenerator: Sync + Send {
    /// Returns a code generated using the given random number generator `rng`.
    fn generate_with_rng<R: Rng>(&self, rng: &mut R) -> ParityCheckMatrix;

    /// Returns a code generated using the thread rng.
    fn generate(&self) -> ParityCheckMatrix {
        self.generate_with_rng(&mut thread_rng())
    }
}
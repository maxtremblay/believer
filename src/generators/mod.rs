use crate::ErasureDecoder;
use crate::ParityCheckMatrix;
use crate::SimulationResult;
use crate::Simulator;
use rand::{thread_rng, Rng};

pub trait CodeGenerator {
    /// Returns a code generated using the given random number generator `rng`.
    fn generate_with_rng<R: Rng>(&self, rng: R) -> ParityCheckMatrix;

    fn generate(&self) -> ParityCheckMatrix {
        self.generate_with_rng(thread_rng())
    }

    /// Returns the code with the best performance among the `n_codes` generated codes. That is,
    /// the code that took the most iterations to obtain `n_failures_per_code` failures using an
    /// erausure decoder with probability `erasure_prob`.
    fn find_best_code_from_erasure(
        &self,
        erasure_prob: f64,
        n_codes: usize,
        n_failures_per_code: usize,
    ) -> Option<(ParityCheckMatrix, SimulationResult)> {
        let mut best_code = None;
        let mut best_performance = SimulationResult::worse_result();
        for _ in 0..n_codes {
            let code = self.generate();
            let decoder = ErasureDecoder::new(code.clone(), erasure_prob);
            let simulator = Simulator::new(decoder);
            let result = simulator.simulate_until_events_are_found(n_failures_per_code);
            if result.failure_rate() < best_performance.failure_rate() {
                best_performance = result;
                best_code = Some(code);
            }
        }
        best_code.map(|code| (code, best_performance))
    }
}

pub mod random_checks;

// pub mod regular_ldpc;
// pub use regular_ldpc::*;

pub mod hierarchical_codes;
// pub use hierarchical_code::*;

pub mod increasing_range_code;
pub use increasing_range_code::*;

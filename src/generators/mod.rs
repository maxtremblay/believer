use crate::ParityCheckMatrix;
use crate::ErasureDecoder;
use crate::Simulator;
use crate::SimulationResult;

pub trait CodeGenerator {
    fn generate(&self) -> ParityCheckMatrix;

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
            let decoder = ErasureDecoder::new(&code, erasure_prob);
            let simulator = Simulator::new(&decoder);
            let result = simulator.simulate_until_failures_are_found(n_failures_per_code, 1);
            if result > best_performance {
                best_performance = result;
                best_code = Some(code);
            }
        }
        best_code.map(|code| (code, best_performance))
    }
}
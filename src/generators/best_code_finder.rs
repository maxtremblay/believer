use super::CodeGenerator;
use crate::{Decoder, ParityCheckMatrix, SimulationResult, ErasureDecoder};
use rand::{Rng, SeedableRng};
use rand::distributions::Standard;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

type CodeAndResult = (Option<ParityCheckMatrix>, SimulationResult);

pub struct BestCodeFinderFromErasure<'a, G: CodeGenerator> {
    generator: &'a G,
    erasure_prob: f64,
    n_codes_to_try: usize,
}

impl<'a, G: CodeGenerator> BestCodeFinderFromErasure<'a, G> {
    pub fn from_generator(generator: &'a G) -> Self {
        Self {
            generator,
            erasure_prob: 0.5,
            n_codes_to_try: 0,
        }
    }

    pub fn among_n_codes(mut self, n_codes: usize) -> Self {
        self.n_codes_to_try = n_codes;
        self
    }

    pub fn with_erasure_prob(mut self, prob: f64) -> Self {
        if prob < 0.0 || prob > 1.0 { panic!("prob is not between 0 and 1") }
        self.erasure_prob = prob;
        self
    }

    pub fn find_best_code_simulating_n_iterations_with_rng<R: Rng>(
        &self,
        n_iterations: usize,
        rng: &mut R,
    ) -> (Option<ParityCheckMatrix>, SimulationResult) {
        NIterationsBestCodeFinderFromErasure::from(self)
            .with_n_iterations(n_iterations)
            .find_with_rng(rng)  
    }
}

struct NIterationsBestCodeFinderFromErasure<'a, G: CodeGenerator> {
    code_finder: &'a BestCodeFinderFromErasure<'a, G>,
    n_iterations: usize,
    random_seeds: Vec<u64>,
}

impl<'a, G: CodeGenerator> NIterationsBestCodeFinderFromErasure<'a, G> {
    fn from(code_finder: &'a BestCodeFinderFromErasure<'a, G>) -> Self {
        Self { code_finder, n_iterations: 0, random_seeds: Vec::new() }
    }

    fn with_n_iterations(mut self, n_iterations: usize) -> Self {
        self.n_iterations = n_iterations;
        self
    }

    fn find_with_rng<R: Rng>(mut self, rng: &mut R) -> CodeAndResult {
        self.initialize_random_seeds_with_rng(rng);
        (0..self.code_finder.n_codes_to_try)
            .into_par_iter()
            .map(|code_index| {
                let mut rng = self.get_rng_for(code_index);
                self.simulate_one_code_with_rng(&mut rng)
            })
            .reduce(
                || (None, SimulationResult::worse_result()),
                |accumulator, code_and_result| Self::get_best_between(accumulator, code_and_result),
        )
    }

    fn initialize_random_seeds_with_rng<R: Rng>(&mut self, rng: &mut R) {
        self.random_seeds = rng
            .sample_iter(Standard)
            .take(self.code_finder.n_codes_to_try)
            .collect()
    }

    fn get_rng_for(&self, index: usize) -> ChaCha8Rng {
        ChaCha8Rng::seed_from_u64(self.random_seeds[index])
    }

    fn simulate_one_code_with_rng<R: Rng>(&self, rng: &mut R) -> CodeAndResult {
        let code = self.code_finder.generator.generate_with_rng(rng);
        let mut decoder = ErasureDecoder::with_prob(self.code_finder.erasure_prob).for_code(code);
        let result = decoder.simulate_n_iterations_with_rng(self.n_iterations, rng);
        (Some(decoder.take_code()), result)
    }

    fn get_best_between(first: CodeAndResult, second: CodeAndResult) -> CodeAndResult {
        if first.1.is_better_than(&second.1) { first } else { second }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;
    use super::super::RegularLDPCCodeGenerator;

    #[test]
    fn reproductibility_for_finding_best_ldpc_code_using_n_iterations() {
        let rng = ChaCha8Rng::seed_from_u64(123);
        let generator = RegularLDPCCodeGenerator::new(3, 4, 2, 4);

        let code_finder = BestCodeFinderFromErasure::from_generator(&generator)
            .with_erasure_prob(0.5)
            .among_n_codes(10);

        let code_and_result_0 = code_finder
            .find_best_code_simulating_n_iterations_with_rng(50, &mut rng.clone());
        
        let code_and_result_1 = code_finder
            .find_best_code_simulating_n_iterations_with_rng(50, &mut rng.clone());

        assert_eq!(code_and_result_0, code_and_result_1);
    }
}

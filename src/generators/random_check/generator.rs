use super::adjacency::Adjacency;
use rand::distributions::WeightedIndex;
use rand::Rng;

/// A `Generator` helps generating checks for a code while respecting some global
/// constraints.
///
/// Constraints are bit and check degrees and minimal girth. A `Generator` is consumed
/// while generating checks and need to be reset before being use to generate checks for another
/// code.
pub struct Generator {
    max_bit_degrees: usize,
    bit_degrees: Vec<usize>,
    adjacency: Adjacency,
    active_bits: Vec<usize>,
    distribution: Vec<f64>,
}

impl Generator {
    // ************ 
    // Construction
    // ************

    pub fn new() -> Self {
        Self {
            max_bit_degrees: 0,
            bit_degrees: Vec::new(),
            adjacency: Adjacency::new(),
            active_bits: Vec::new(),
            distribution: Vec::new(),
        }
    }

    pub fn with_minimal_girth(mut self, girth: usize) -> Self {
        self.adjacency.set_recursion_depth_from_girth(girth);
        self
    }

    pub fn with_max_bit_degrees(mut self, degrees: usize) -> Self {
        self.max_bit_degrees = degrees;
        self
    }

    pub fn with_n_bits(mut self, n_bits: usize) -> Self {
        self.set_bit_degrees(n_bits);
        self.adjacency.initialize_adjacencies(n_bits);
        self.set_active_bits(n_bits);
        self.set_distribution(n_bits);
        self
    }

    fn set_bit_degrees(&mut self, n_bits: usize) {
        self.bit_degrees = vec![0; n_bits];
    }

    fn set_active_bits(&mut self, n_bits: usize) {
        self.active_bits = (0..n_bits).collect();
    }

    fn set_distribution(&mut self, n_bits: usize) {
        self.distribution = vec![1.0 / n_bits as f64; n_bits];
    }

    // **************
    // Public methods
    // **************

    /// Returns the list of bits adjacent to `bit` given the minimal girth of the generator.
    ///
    /// Two bits are adjacent if connecting them to the same check will create a cycle smaller than
    /// the minimal girth.
    pub fn adjacent_to(&self, bit: usize) -> Vec<usize> {
        self.adjacency.get_bits_adjacent_to(bit)
    }

    /// Generates a random check of degree `check_degree` using the random number generator `rng`.
    pub fn generate<R: Rng + ?Sized>(
        &mut self,
        check_degree: usize,
        rng: &mut R,
    ) -> Option<Vec<usize>> {
        let mut check = Vec::with_capacity(check_degree);
        for _ in 0..check_degree {
            self.add_random_bit_to_check(&mut check, rng);
        }
        if check.len() >= 2 {
            self.update_adjacency_from_check(&check);
            check.iter().for_each(|bit| self.bit_degrees[*bit] += 1);
            check.sort();
            Some(check)
        } else {
            None
        }
    }

    pub fn get_n_bits(&self) -> usize {
        self.bit_degrees.len()
    }

    pub fn over_all_bits(&mut self) -> &mut Self {
        self.active_bits = (0..self.get_n_bits()).collect();
        self
    }

    pub fn over_bits(&mut self, mut bits: Vec<usize>) -> &mut Self {
        bits.sort();
        bits.dedup();
        self.active_bits = bits;
        self
    }

    pub fn reset(&mut self) {
        self.set_bit_degrees(self.get_n_bits());
        self.adjacency.initialize_adjacencies(self.get_n_bits());
    }


    pub fn with_distribution(&mut self, probs: Vec<f64>) -> &mut Self {
        if probs.len() != self.get_n_bits() {
            panic!("wrong number of probabilities");
        }
        self.distribution = probs;
        self
    }

    pub fn with_uniform_distribution(&mut self) -> &mut Self {
        self.distribution = vec![1.0 / self.get_n_bits() as f64; self.get_n_bits()];
        self
    }

    pub fn without(&mut self, bits: &[usize]) -> &mut Self {
        bits.iter().for_each(|bit| {
            if let Some(index) = self.active_bits.iter().position(|b| b == bit) {
                self.active_bits.swap_remove(index);
            }
        });
        self
    }

    // ***************
    // Private methods
    // ***************

    fn add_random_bit_to_check<R: Rng + ?Sized>(&self, check: &mut Vec<usize>, rng: &mut R) {
        let availables: Vec<usize> = self
            .active_bits
            .iter()
            .filter(|bit| self.bit_degrees[**bit] < self.max_bit_degrees)
            .filter(|bit| check.iter().all(|b| self.are_not_adjacent(b, bit)))
            .cloned()
            .collect();
        let probs: Vec<f64> = availables
            .iter()
            .map(|bit| self.distribution[*bit])
            .collect();
        if probs.len() > 0 && probs.iter().sum::<f64>() > 0.0 {
            let distribution = WeightedIndex::new(&probs).unwrap();
            check.push(availables[rng.sample(distribution)]);
        }
    }

    fn are_not_adjacent(&self, bit_0: &usize, bit_1: &usize) -> bool {
        !self.adjacent_to(*bit_0).contains(bit_1)
    }

    fn update_adjacency_from_check(&mut self, check: &[usize]) {
        self.adjacency.update_from_check(check)
    }
}
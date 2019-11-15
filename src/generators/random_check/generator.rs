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
        self.adjacency.set_recursion_depth(girth / 2);
        self
    }

    pub fn with_max_bit_degrees(mut self, degrees: usize) -> Self {
        self.max_bit_degrees = degrees;
        self
    }

    pub fn with_n_bits(mut self, n_bits: usize) -> Self {
        self.set_bit_degrees(n_bits);
        self.set_adjacency(n_bits);
        self.set_active_bits(n_bits);
        self.set_distribution(n_bits);
        self
    }

    fn set_bit_degrees(&mut self, n_bits: usize) {
        self.bit_degrees = vec![0; n_bits];
    }

    fn set_adjacency(&mut self, n_bits: usize) {
        self.adjacency = Adjacency::with_n_bits(n_bits);
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
        let n_bits = self.get_n_bits();
        self.set_bit_degrees(n_bits);
        self.set_adjacency(n_bits);
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

#[cfg(test)]
mod test {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn doesnt_include_same_bit_twice() {
        let mut rng = ChaCha8Rng::seed_from_u64(10);

        let mut generator = Generator::new()
            .with_n_bits(3)
            .with_max_bit_degrees(2);

        let first_check = generator.over_bits(vec![0, 1]).generate(2, &mut rng);
        assert_eq!(first_check, Some(vec![0, 1]));

        // Can't generate a degree 4 check over only 2 bits.
        let second_check = generator.generate(4, &mut rng);
        assert_eq!(second_check, None);
    }

    #[test]
    fn doesnt_exceed_bit_maximal_degree() {
        let mut rng = ChaCha8Rng::seed_from_u64(10);

        let mut generator = Generator::new()
            .with_n_bits(3)
            .with_max_bit_degrees(2);

        let first_check = generator.without(&[2]).generate(2, &mut rng);
        assert_eq!(first_check, Some(vec![0, 1]));

        let second_check = generator
            .over_all_bits()
            .without(&[0])
            .generate(2, &mut rng);
        assert_eq!(second_check, Some(vec![1, 2]));

        // We already have checks [0,1] and [1, 2]. Degree of bit 1 is 2 and it can't
        // be included in another check.

        let third_check = generator.over_all_bits().generate(3, &mut rng);
        assert_eq!(third_check, None);

        let fourth_check = generator.generate(2, &mut rng);
        assert_eq!(fourth_check, Some(vec![0, 2]));

        // Every bit has max degree. Can't generate anymore check
        assert_eq!(generator.generate(1, &mut rng), None);
    }

    #[test]
    fn doesnt_create_cycle_smaller_than_minimal_girth() {
        let mut rng = ChaCha8Rng::seed_from_u64(10);

        // Minimal girth 6
        let mut generator = Generator::new()
            .with_n_bits(3)
            .with_max_bit_degrees(2)
            .with_minimal_girth(6);

        let first_check = generator.generate(3, &mut rng);
        assert_eq!(first_check, Some(vec![0, 1, 2]));

        // Any check of degree 2 will create a 4-cycle.
        let second_check = generator.generate(2, &mut rng);
        assert_eq!(second_check, None);

        // A degree 1 check will not create a 4-cycle.
        let third_check = generator.generate(1, &mut rng);
        assert_eq!(third_check.is_some(), true);

        // Minimal girth 8
        let mut generator = Generator::new()
            .with_n_bits(5)
            .with_max_bit_degrees(2)
            .with_minimal_girth(8);

        let first_check = generator.over_bits(vec![0, 1, 2]).generate(3, &mut rng);
        assert_eq!(first_check, Some(vec![0, 1, 2]));

        let second_check = generator.over_bits(vec![2, 3]).generate(2, &mut rng);
        assert_eq!(second_check, Some(vec![2, 3]));

        // A check over [0, 3] will create a 6-cycle.
        let third_check = generator.over_bits(vec![0, 3]).generate(2, &mut rng);
        assert_eq!(third_check, None);

        // Possible checks are [0, 4] or [3, 4]
        let fourth_check = generator.over_bits(vec![0, 3, 4]).generate(2, &mut rng);
        assert_eq!(fourth_check.clone().unwrap().contains(&4), true);
        assert_eq!(fourth_check.unwrap().len(), 2);
    }

    #[test]
    fn generate_bit_according_to_distribution() {
        let mut rng = ChaCha8Rng::seed_from_u64(10);

        let mut generator = Generator::new()
            .with_n_bits(5)
            .with_max_bit_degrees(2);
            
        generator
            .with_distribution(vec![0.25, 0.25, 0.0, 0.25, 0.25])
            .over_bits(vec![0, 1, 2]);

        // Can't generate 3 bits from this distribution over the first 3.
        assert_eq!(generator.generate(3, &mut rng), None);
        assert_eq!(generator.generate(2, &mut rng), Some(vec![0, 1]));
        assert_eq!(generator.generate(2, &mut rng), Some(vec![0, 1]));

        // Degree of the first 2 bits is 2.
        assert_eq!(generator.generate(2, &mut rng), None);

        generator.over_all_bits();
        assert_eq!(generator.generate(2, &mut rng), Some(vec![3, 4]));

        // Can't pick a degree 3 check because probability of bit 2 is 0.
        assert_eq!(generator.generate(3, &mut rng), None);

        // Reset distribution.
        generator.with_uniform_distribution();
        assert_eq!(generator.generate(3, &mut rng), Some(vec![2, 3, 4]));
    }
}
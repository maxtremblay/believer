use rand::distributions::WeightedIndex;
use rand::Rng;
use std::collections::{BTreeMap, BTreeSet};

/// A `RandomCheckGenerator` helps generating checks for a code while respecting some global
/// constraints.
///
/// Constraints are bit and check degrees and minimal girth. A `RandomCheckGenerator` is consumed
/// while generating checks and need to be reset before being use to generate checks for a other
/// code.
pub struct RandomCheckGenerator {
    max_bit_degrees: Vec<usize>,
    bit_degrees: Vec<usize>,
    adjacencies: Vec<BTreeSet<usize>>,
    adjacency_depth: usize,
    active_bits: Vec<usize>,
    distribution: Vec<f64>,
}

impl RandomCheckGenerator {
    // **************
    // Public methods
    // **************

    /// Returns the list of bits adjacent to `bit` given the minimal girth of the generator.
    ///
    /// Two bits are adjacent if connecting them to the same check will create a cycle smaller than
    /// the minimal girth.
    pub fn adjacent_to(&self, bit: usize) -> Vec<usize> {
        if bit >= self.n_bits() {
            Vec::new()
        } else if self.adjacency_depth == 0 {
            vec![bit]
        } else {
            let mut adjacents = BTreeMap::new();
            self.adjacent_to_recursion(bit, self.adjacency_depth, &mut adjacents);
            adjacents.keys().cloned().collect()
        }
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
        if check.len() == check_degree {
            self.update_adjacency(&check);
            check.iter().for_each(|bit| self.bit_degrees[*bit] += 1);
            check.sort();
            Some(check)
        } else {
            None
        }
    }

    pub fn n_bits(&self) -> usize {
        self.bit_degrees.len()
    }

    pub fn new(max_bit_degrees: Vec<usize>, minimal_girth: usize) -> Self {
        let n_bits = max_bit_degrees.len();
        Self {
            max_bit_degrees,
            bit_degrees: vec![0; n_bits],
            adjacency_depth: minimal_girth / 2,
            adjacencies: (0..n_bits).map(|b| [b].iter().cloned().collect()).collect(),
            active_bits: (0..n_bits).collect(),
            distribution: vec![1.0 / n_bits as f64; n_bits],
        }
    }

    pub fn over_all_bits(&mut self) -> &mut Self {
        self.active_bits = (0..self.n_bits()).collect();
        self
    }

    pub fn over_bits(&mut self, mut bits: Vec<usize>) -> &mut Self {
        bits.sort();
        bits.dedup();
        self.active_bits = bits;
        self
    }

    pub fn reset(&mut self) {
        self.bit_degrees = vec![0; self.n_bits()];
        self.adjacencies = (0..self.n_bits())
            .map(|b| [b].iter().cloned().collect())
            .collect();
    }

    pub fn with_distribution(&mut self, probs: Vec<f64>) -> &mut Self {
        if probs.len() != self.n_bits() {
            panic!("wrong number of probabilities");
        }
        self.distribution = probs;
        self
    }

    pub fn with_uniform_distribution(&mut self) -> &mut Self {
        self.distribution = vec![1.0 / self.n_bits() as f64; self.n_bits()];
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
            .filter(|bit| self.bit_degrees[**bit] < self.max_bit_degrees[**bit])
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

    // Helper function for `self.adjacent_to(...)`.
    fn adjacent_to_recursion(
        &self,
        bit: usize,
        depth: usize,
        adjacents: &mut BTreeMap<usize, usize>,
    ) {
        adjacents.insert(bit, depth);
        if depth > 0 {
            self.adjacencies[bit].iter().for_each(|b| {
                if adjacents.get(&b).map(|d| *d < depth).unwrap_or(true) {
                    self.adjacent_to_recursion(*b, depth - 1, adjacents);
                }
            })
        }
    }

    fn are_not_adjacent(&self, bit_0: &usize, bit_1: &usize) -> bool {
        !self.adjacent_to(*bit_0).contains(bit_1)
    }

    fn update_adjacency(&mut self, check: &[usize]) {
        check.iter().for_each(|bit_0| {
            check.iter().for_each(|bit_1| {
                self.adjacencies[*bit_1].insert(*bit_0);
            });
        });
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_adjacencies() {
        let mut rng = ChaCha8Rng::seed_from_u64(10);

        // A 15 bit check generator with minimal girth 6;
        let mut generator = RandomCheckGenerator::new(vec![3; 15], 6);

        // Each bit is adjacent to itself.
        (0..15).for_each(|b| assert_eq!(generator.adjacent_to(b), vec![b]));

        // Generate some random weight 3 checks.
        assert_eq!(generator.generate(3, &mut rng), Some(vec![0, 8, 9]));
        assert_eq!(generator.generate(3, &mut rng), Some(vec![2, 10, 12]));
        assert_eq!(generator.generate(3, &mut rng), Some(vec![1, 4, 12]));
        assert_eq!(generator.generate(3, &mut rng), Some(vec![0, 1, 11]));
        assert_eq!(generator.generate(3, &mut rng), Some(vec![0, 5, 6]));

        // Check adjacencies
        assert_eq!(
            generator.adjacent_to(0),
            vec![0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12]
        );
        assert_eq!(
            generator.adjacent_to(1),
            vec![0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12]
        );
        assert_eq!(generator.adjacent_to(2), vec![0, 1, 2, 4, 10, 11, 12]);
        assert_eq!(generator.adjacent_to(3), vec![3]);
        assert_eq!(
            generator.adjacent_to(4),
            vec![0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12]
        );
        assert_eq!(generator.adjacent_to(5), vec![0, 1, 4, 5, 6, 8, 9, 11, 12]);
        assert_eq!(generator.adjacent_to(6), vec![0, 1, 4, 5, 6, 8, 9, 11, 12]);
        assert_eq!(generator.adjacent_to(7), vec![7]);
        assert_eq!(generator.adjacent_to(8), vec![0, 1, 4, 5, 6, 8, 9, 11, 12]);
        assert_eq!(generator.adjacent_to(9), vec![0, 1, 4, 5, 6, 8, 9, 11, 12]);
        assert_eq!(generator.adjacent_to(10), vec![0, 1, 2, 4, 10, 11, 12]);
        assert_eq!(
            generator.adjacent_to(11),
            vec![0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12]
        );
        assert_eq!(
            generator.adjacent_to(12),
            vec![0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12]
        );
        assert_eq!(generator.adjacent_to(13), vec![13]);
        assert_eq!(generator.adjacent_to(14), vec![14]);
    }

    #[test]
    fn doesnt_include_same_bit_twice() {
        let mut rng = ChaCha8Rng::seed_from_u64(10);

        let mut generator = RandomCheckGenerator::new(vec![2, 2, 2], 0);
        let first_check = generator.over_bits(vec![0, 1]).generate(2, &mut rng);
        assert_eq!(first_check, Some(vec![0, 1]));

        // Can't generate a degree 4 check over only 2 bits.
        let second_check = generator.generate(4, &mut rng);
        assert_eq!(second_check, None);
    }

    #[test]
    fn doesnt_exceed_bit_maximal_degree() {
        let mut rng = ChaCha8Rng::seed_from_u64(10);

        let mut generator = RandomCheckGenerator::new(vec![2, 2, 2], 0);

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
        let mut generator = RandomCheckGenerator::new(vec![2, 2, 2], 6);

        let first_check = generator.generate(3, &mut rng);
        assert_eq!(first_check, Some(vec![0, 1, 2]));

        // Any check of degree 2 will create a 4-cycle.
        let second_check = generator.generate(2, &mut rng);
        assert_eq!(second_check, None);

        // A degree 1 check will not create a 4-cycle.
        let third_check = generator.generate(1, &mut rng);
        assert_eq!(third_check.is_some(), true);

        // Minimal girth 8
        let mut generator = RandomCheckGenerator::new(vec![2; 5], 8);

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

        let mut generator = RandomCheckGenerator::new(vec![2; 5], 0);
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

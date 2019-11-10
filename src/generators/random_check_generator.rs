use rand::distributions::WeightedIndex;
use rand::Rng;

/// A `RandomCheckGenerator` helps generating checks for a code while respecting some global
/// constraints.
///
/// Constraints are bit and check degrees and minimal girth.
pub struct RandomCheckGenerator {
    max_bit_degrees: Vec<usize>,
    bit_degrees: Vec<usize>,
    adjacencies: Vec<Vec<usize>>,
    ajdacency_depth: usize,
    active_bits: Vec<usize>,
    distribution: Vec<f64>,
}

impl RandomCheckGenerator {
    pub fn generate<R: Rng + ?Sized>(&mut self, check_degree: usize, rng: &mut R) -> Option<Vec<usize>> {
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

    fn add_random_bit_to_check<R: Rng + ?Sized>(&self, check: &mut Vec<usize>, rng: &mut R) {
        let availables: Vec<usize> = self
            .active_bits
            .iter()
            .filter(|bit| self.bit_degrees[**bit] < self.max_bit_degrees[**bit])
            .filter(|bit| !check.contains(bit))
            .filter(|bit| check.iter().all(|b| self.are_not_adjacent(b, bit)))
            .map(|bit| *bit)
            .collect();
        let probs: Vec<f64> = availables
            .iter()
            .map(|bit| self.distribution[*bit])
            .collect();
        if probs.len() > 0 && probs.iter().sum::<f64>() > 0.0 {
            let distribution = WeightedIndex::new(Self::normalize(&probs)).unwrap();
            check.push(availables[rng.sample(distribution)]);
        }
    }

    fn are_not_adjacent(&self, bit_0: &usize, bit_1: &usize) -> bool {
        !self.adjacencies[*bit_0].contains(bit_1) && !self.adjacencies[*bit_1].contains(bit_0)
    }

    fn update_adjacency(&mut self, check: &[usize]) {
        check.iter().for_each(|bit_0| {
            check.iter().for_each(|bit_1| {
                let adjacent_bits = self.adjacent_bits(bit_0, self.ajdacency_depth);
                self.adjacencies[*bit_1].extend_from_slice(&adjacent_bits);
            });
        });
    }

    fn adjacent_bits(&self, bit: &usize, depth: usize) -> Vec<usize> {
        if depth == 0 {
            Vec::new()
        } else if depth == 1 {
            vec![*bit]
        } else {
            self.adjacencies[*bit]
                .iter()
                .flat_map(|b| self.adjacent_bits(b, depth - 1))
                .collect()
        }
    }

    pub fn new(max_bit_degrees: Vec<usize>, minimal_girth: usize) -> Self {
        let n_bits = max_bit_degrees.len();
        Self {
            max_bit_degrees,
            bit_degrees: vec![0; n_bits],
            ajdacency_depth: if minimal_girth > 4 {
                (minimal_girth - 4) / 2
            } else {
                0
            },
            adjacencies: (0..n_bits).map(|b| vec![b]).collect(),
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

    pub fn with_distribution(&mut self, probs: &[f64]) -> &mut Self {
        if probs.len() != self.n_bits() {
            panic!("wrong number of probabilities");
        }
        self.distribution = Self::normalize(probs);
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

    pub fn n_bits(&self) -> usize {
        self.bit_degrees.len()
    }

    fn normalize(probs: &[f64]) -> Vec<f64> {
        let sum: f64 = probs.iter().sum();
        if sum == 0.0 {
            probs.to_vec()
        } else {
            probs.iter().map(|p| p / sum).collect()
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand_chacha::ChaCha8Rng;
    use rand::SeedableRng;

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

        let second_check = generator.over_all_bits().without(&[0]).generate(2, &mut rng);
        assert_eq!(second_check, Some(vec![1, 2]));

        // We already have checks [0,1] and [1, 2]. Degree of bit 1 is 2 and it can't be include in
        // another check.

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
            .with_distribution(&[0.25, 0.25, 0.0, 0.25, 0.25])
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

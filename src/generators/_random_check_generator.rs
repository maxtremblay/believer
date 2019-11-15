use rand::distributions::WeightedIndex;
use rand::Rng;
use std::collections::{BTreeMap, BTreeSet};

/// A `RandomCheckGenerator` helps generating checks for a code while respecting some global
/// constraints.
///
/// Constraints are bit and check degrees and minimal girth. A `RandomCheckGenerator` is consumed
/// while generating checks and need to be reset before being use to generate checks for another
/// code.
pub struct RandomCheckGenerator {
    max_bit_degrees: usize,
    bit_degrees: Vec<usize>,
    adjacency: Adjacency,
    active_bits: Vec<usize>,
    distribution: Vec<f64>,
}

impl RandomCheckGenerator {
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

struct Adjacency {
    adjacencies: Vec<BTreeSet<usize>>,
    depth: usize,
}

impl Adjacency {

    // ***** Construction *****

    fn new() -> Self {
        Self {
            depth: 0,
            adjacencies: Vec::new(),
        }
    }

    fn initialize_adjacencies(&mut self, n_bits: usize) {
        self.adjacencies = (0..n_bits)
            .map(|b| {
                let mut set = BTreeSet::new();
                set.insert(b);
                set
            })
            .collect();
    }

    fn set_recursion_depth_from_girth(&mut self, girth: usize) {
        self.depth = girth / 2;
    }

    // ***** Get adjacent bits *****

    fn get_bits_adjacent_to(&self, bit: usize) -> Vec<usize> {
        if self.is_out_of_bound(bit) {
            Vec::new()
        } else {
            self.get_bits_adjacent_to_inbound_bit(bit)
        }
    }

    fn is_out_of_bound(&self, bit:usize) -> bool {
        bit >= self.get_n_bits()
    }

    fn get_n_bits(&self) -> usize {
        self.adjacencies.len()
    }

    fn get_bits_adjacent_to_inbound_bit(&self, bit: usize) -> Vec<usize> {
        if self.depth == 0 {
            vec![bit]
        } else {
            let mut getter = self.initialize_adjacent_bits_getter(bit);
            getter.get_adjacent_bits_at_depth(self.depth)
        }
    }

    fn initialize_adjacent_bits_getter(&self, source_bit: usize) -> AdjacentBitsGetter {
        AdjacentBitsGetter {
            source_bit,
            adjacent_bits: BTreeMap::new(),
            adjacencies: &self.adjacencies,
        }
    }

    // ***** Update adjacent bits *****

    fn update_from_check(&mut self, check: &[usize]) {
        check.iter().for_each(|bit| self.set_bits_in_check_adjacent_to(*bit, check));
    }

    fn set_bits_in_check_adjacent_to(&mut self, source_bit: usize, check: &[usize]) {
        check.iter().for_each(|adjacent_bit| {
            self.adjacencies[source_bit].insert(*adjacent_bit);
        });
    }
}

struct AdjacentBitsGetter<'a> {
    source_bit: usize,
    adjacent_bits: BTreeMap<usize, usize>,
    adjacencies: &'a [BTreeSet<usize>]
}

impl<'a> AdjacentBitsGetter<'a> {
    fn get_adjacent_bits_at_depth(&mut self, depth: usize) -> Vec<usize> {
        if self.adjacent_bits.is_empty() {
            self.recursively_set_adjacent_bits(self.source_bit, depth);
        }
        self.get_adjacent_bits_as_vec()
    }

    fn recursively_set_adjacent_bits(&mut self, bit: usize, depth: usize) {
        self.adjacent_bits.insert(bit, depth);
        if depth > 0 {
            self.adjacencies[bit].iter().for_each(|b| {
                if self.has_not_been_set_at_depth(b, depth) {
                    self.recursively_set_adjacent_bits(*b, depth - 1);
                }
            })
        }
    }

    fn has_not_been_set_at_depth(&self, bit: &usize, depth: usize) -> bool {
        self.adjacent_bits[bit] < depth
    }

    fn get_adjacent_bits_as_vec(&self) -> Vec<usize> {
        self.adjacent_bits.keys().cloned().collect()
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
        let mut generator = RandomCheckGenerator::new()
            .with_max_bit_degrees(3)
            .with_n_bits(6)
            .with_minimal_girth(6);

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

        let mut generator = RandomCheckGenerator::new()
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

        let mut generator = RandomCheckGenerator::new()
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
        let mut generator = RandomCheckGenerator::new()
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
        let mut generator = RandomCheckGenerator::new()
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

        let mut generator = RandomCheckGenerator::new()
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

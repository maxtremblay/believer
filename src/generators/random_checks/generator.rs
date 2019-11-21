use super::adjacency::Adjacency;
use rand::distributions::WeightedIndex;
use rand::rngs::ThreadRng;
use rand::{thread_rng, Rng};

/// A `Generator` helps generating checks for a code while respecting some global constraints.
///
/// Constraints are degrees and minimal girth. A `Generator` is consumed while generating checks.
/// Take a look at all the setters methods (the one that have `&mut self` as argument and output)
/// to see all the parameters that can be set.
///
/// # Example
///
/// ```
/// use believer::random_checks::Generator;
/// use rand::SeedableRng;
/// use rand_chacha::ChaCha8Rng;
///
/// // Create a check generator for 10 bits with the given random number generator.
/// let mut generator = Generator::with_n_bits(10)
///     .with_random_number_generator(ChaCha8Rng::seed_from_u64(123));
///
/// // Set some parameters for the generator.
/// generator
///     .set_maximal_bit_degree(3)
///     .set_minimal_girth(6)
///     .set_target_check_degree(4);
///
/// // Generate a check over the first 4 bits (this is in fact deterministic).
/// let first_check = generator.set_over_bits(vec![0, 1, 2, 3]).get_random_check();
/// assert_eq!(first_check, Some(vec![0, 1, 2, 3]));
///
/// // Use relative weight for the next check.
/// let distribution = vec![0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2];
/// let second_check = generator
///     .set_over_all_bits()
///     .set_distribution(distribution)
///     .get_random_check();
/// assert_eq!(second_check, Some(vec![6, 7, 8, 9]));
///
/// // By default, all generated checks of degree below target (4 in this example) are rejected.
/// // We can't generate a degree 4 check over only 3 bits.
/// let third_check = generator
///     .set_uniform_distribution()
///     .set_over_bits(vec![4, 5, 6])
///     .get_random_check();
/// assert!(third_check.is_none());
///
/// // But, we can relax this constraint. In this case, the generator will try to create a degree 4
/// // (the target check degree) check. Since it not possible, a degree 3 check will be returns
/// // instead.
/// generator.allow_checks_of_degree_at_least(3);
/// let fourth_check = generator.get_random_check();
/// assert_eq!(fourth_check, Some(vec![4, 5, 6]));
///
/// // It is also possible that a check is rejected if it would create a cycle smaller than the
/// // minimal girth or if all the bits have reached the maximum degree. In the following example
/// // having bit 4 and 6 together in an other check would create a length 4 cycle.
/// let fifth_check = generator.set_over_bits(vec![0, 4, 6]).get_random_check();
/// assert!(fifth_check.is_none());
///
/// // However, it is possible to generate a degree 2 check.
/// let sixth_check = generator.allow_checks_of_degree_at_least(2).get_random_check();
/// assert_eq!(sixth_check, Some(vec![0, 4]));
/// ```
pub struct Generator<R: Rng = ThreadRng> {
    maximal_bit_degree: usize,
    bit_degrees: Vec<usize>,
    adjacency: Adjacency,
    active_bits: Vec<usize>,
    distribution: Vec<f64>,
    target_check_degree: usize,
    check_degree_condition: CheckDegreeCondition,
    random_number_generator: R,
}

impl Generator<ThreadRng> {
    /// Creates a generator for empty code.
    pub fn new() -> Self {
        Self {
            maximal_bit_degree: 1,
            bit_degrees: Vec::new(),
            adjacency: Adjacency::new(),
            active_bits: Vec::new(),
            distribution: Vec::new(),
            target_check_degree: 2,
            check_degree_condition: CheckDegreeCondition::MustBeFullDegree,
            random_number_generator: thread_rng(),
        }
    }

    /// Creates a generator for `n_bits` without any restriction on the checks that are going
    /// to be generated and using the uniform distribution over all bits.
    pub fn with_n_bits(n_bits: usize) -> Self {
        let mut generator = Self::new()
            .initialize_bit_degrees(n_bits)
            .initialize_adjacency(n_bits);
        generator.set_over_all_bits().set_uniform_distribution();
        generator
    }

    fn initialize_bit_degrees(mut self, n_bits: usize) -> Self {
        self.bit_degrees = vec![0; n_bits];
        self
    }

    fn initialize_adjacency(mut self, n_bits: usize) -> Self {
        self.adjacency = Adjacency::with_n_bits(n_bits);
        self
    }

    // fn initialize_active_bits(mut self, n_bits: usize) -> Self {
    //     self.active_bits = (0..n_bits).collect();
    //     self
    // }

    // fn initialize_distribution(mut self, n_bits: usize) -> Self {
    //     self.distribution = vec![1.0 / n_bits as f64; n_bits];
    //     self
    // }
}

impl<R: Rng> Generator<R> {
    // ***** Construction *****

    /// Create a copy of `self` with the given `rng` consuming self .
    ///
    /// If not set, random numbers will be generate using `rand::thread_rng`. This method consume
    /// `self` and create a new `Generator` from it taking ownership of most of it's parameters. It
    /// is designed to be use during construction.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::random_checks::Generator;
    ///
    /// // Some random number generator.
    /// use rand::SeedableRng;
    /// use rand_chacha::ChaCha8Rng;
    ///
    /// let mut rng = ChaCha8Rng::seed_from_u64(123);
    /// let generator = Generator::new().with_random_number_generator(rng);
    /// ```
    pub fn with_random_number_generator<S: Rng>(self, rng: S) -> Generator<S> {
        Generator {
            maximal_bit_degree: 0,
            bit_degrees: self.bit_degrees,
            adjacency: self.adjacency,
            active_bits: self.active_bits,
            distribution: self.distribution,
            target_check_degree: self.target_check_degree,
            check_degree_condition: self.check_degree_condition,
            random_number_generator: rng,
        }
    }

    // ***** Setters *****

    /// Set the minimal girth of `self`.
    ///
    /// The checks generated after this won't create any cycle smaller than the minimal girth.
    /// If this is set before generating the first checks, the generated checks will induce
    /// a code with girth at least `minimal_girth`.
    pub fn set_minimal_girth(&mut self, minimal_girth: usize) -> &mut Self {
        let depth = if minimal_girth >= 2 {
            (minimal_girth - 2) / 2
        } else {
            0
        };
        self.adjacency.set_recursion_depth(depth);
        self
    }

    /// Set the maximal bit degree of `self`.
    ///
    /// A given bit will not be part of any new check if it was generated in
    /// `maximal_bit_degree` previous checks. The default value is 1.
    pub fn set_maximal_bit_degree(&mut self, degree: usize) -> &mut Self {
        self.maximal_bit_degree = degree;
        self
    }

    /// Set `self` to generate the following checks without any restriction on the bits.
    ///
    /// This is the default behavior.
    pub fn set_over_all_bits(&mut self) -> &mut Self {
        self.active_bits = (0..self.get_n_bits()).collect();
        self
    }

    /// Set `self` to generate the following checks using only the given `bits`.
    pub fn set_over_bits(&mut self, mut bits: Vec<usize>) -> &mut Self {
        bits.sort();
        bits.dedup();
        self.active_bits = bits;
        self
    }

    /// Set `self` to generate the following checks without the given `bits`.
    ///
    /// This won't reset the target bits.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::random_checks::Generator;
    ///
    /// let mut generator = Generator::with_n_bits(10);
    ///
    /// // This limit the following checks over bits 2, 3, 4, ..., 9.
    /// generator.set_without_bits(vec![0, 1]);
    ///
    /// // This limit the following checks over bits 7, 8 and 9 since 0 and 1 are already removed.
    /// generator.set_without_bits(vec![2, 3, 4, 5, 6]);
    ///
    /// // This limit the following checks over bits 0 and 1.
    /// generator.set_over_bits(vec![0, 1, 2]).set_without_bits(vec![2]);
    /// ```
    pub fn set_without_bits(&mut self, bits: Vec<usize>) -> &mut Self {
        bits.into_iter().for_each(|bit| {
            if let Some(index) = self.active_bits.iter().position(|b| *b == bit) {
                self.active_bits.swap_remove(index);
            }
        });
        self
    }

    /// Set `self` to generate bits in checks according to the given distribution.
    ///
    /// The distribution does not need to be normalize. It will be normalize to sum to 1 except if
    /// the sum is 0. If not set, the uniform distribution is used.
    ///
    /// # Panic
    ///
    /// Panics if the distribution lenght is not the same as the number of bits in `self`. Also
    /// panics if there are some negative probabilities.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::random_checks::Generator;
    /// let mut generator = Generator::with_n_bits(5);
    ///
    /// // It is more like that the first 3 bits are selected in the following checks.
    /// generator.set_distribution(vec![0.3, 0.3, 0.3, 0.05, 0.05]);
    ///
    /// // It can be set with relative weights.
    /// generator.set_distribution(vec![30.0, 30.0, 30.0, 5.0, 5.0]);
    ///
    /// // If we limit the bits to a subset, the distribution is updated in consequence. In the
    /// // following, the distribution is updated to [0.75, 0.0, 0.0, 12.5, 12.5].
    /// generator.set_without_bits(vec![1, 2]);
    /// ```
    pub fn set_distribution(&mut self, distribution: Vec<f64>) -> &mut Self {
        if distribution.len() != self.get_n_bits() {
            panic!("wrong number of probabilities");
        }
        if distribution.iter().any(|prob| *prob < 0.0) {
            panic!("there are some negative probabilities");
        }
        self.distribution = distribution;
        self
    }

    /// Set `self` to generate bits in checks according to the uniform distribution.
    ///
    /// This is the default behavior.
    pub fn set_uniform_distribution(&mut self) -> &mut Self {
        self.distribution = vec![1.0 / self.get_n_bits() as f64; self.get_n_bits()];
        self
    }

    /// Set the target check degree of `self`.
    ///
    /// By default, the target degree is 2 and a check will be rejected if, due to other
    /// constraints, it can't be generated with the target degree. This can be change using
    /// `self.allow_checks_of_degree_at_least(degree)` or `self.allow_only_full_degree_check()`.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::random_checks::Generator;
    ///
    /// let mut generator = Generator::with_n_bits(5);
    /// generator.set_target_check_degree(4);
    ///
    /// // Can't generator a degree 4 check over only bit 0 and 1.
    /// assert_eq!(generator.set_over_bits(vec![0, 1]).get_random_check(), None);
    /// ```
    pub fn set_target_check_degree(&mut self, degree: usize) -> &mut Self {
        self.target_check_degree = degree;
        self
    }

    /// Set `self` to allow all checks of degree at least `degree`.
    ///
    /// By default, `self` only accept check of target degree (which is 2 by default).
    ///
    /// # Example
    ///
    /// ```
    /// use believer::random_checks::Generator;
    ///
    /// let mut generator = Generator::with_n_bits(5);
    ///
    /// // Can't generator a degree 2 (default) check over only bit 0.
    /// generator.set_over_bits(vec![0]);
    /// assert_eq!(generator.get_random_check(), None);
    ///
    /// // It is possible to allow check of degree 1 checks.
    /// generator.allow_checks_of_degree_at_least(1);
    /// assert_eq!(generator.get_random_check(), Some(vec![0]));
    /// ```
    pub fn allow_checks_of_degree_at_least(&mut self, degree: usize) -> &mut Self {
        self.check_degree_condition = CheckDegreeCondition::MustBeAtLeastOfDegree(degree);
        self
    }

    /// Set `self` to allow only check of target degree (2 by default).
    ///
    /// This is the default behavior.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::random_checks::Generator;
    ///
    /// let mut generator = Generator::with_n_bits(5);
    /// generator.set_target_check_degree(4);
    ///
    /// // Not necessary, this is the default.
    /// generator.allow_only_check_of_target_degree();
    ///
    /// // A degree 3 check is not allowed.
    /// generator.set_over_bits(vec![0, 1, 2]);
    /// assert_eq!(generator.get_random_check(), None);
    ///
    /// // It is possible to generate a check of degree 4 over 5 bits.
    /// generator.set_over_all_bits();
    /// assert!(generator.get_random_check().is_some());
    /// ```
    pub fn allow_only_check_of_target_degree(&mut self) -> &mut Self {
        self.check_degree_condition = CheckDegreeCondition::MustBeFullDegree;
        self
    }

    // ***** Getters *****

    /// Returns the list of bits adjacent to `bit` given the minimal girth of the generator.
    ///
    /// Two bits are adjacent if connecting them to the same check will create a cycle smaller than
    /// the minimal girth.
    pub fn get_bits_adjacent_to(&self, bit: usize) -> Vec<usize> {
        self.adjacency.get_bits_adjacent_to(bit)
    }

    /// Returns the number of bits in the code that `self` is generating checks for.
    pub fn get_n_bits(&self) -> usize {
        self.bit_degrees.len()
    }

    /// Generates a random check.
    ///
    /// The behavior of this function can be change using the setter methods.
    ///
    /// # Example
    ///
    /// ```
    /// use believer::random_checks::Generator;
    /// let mut generator = Generator::with_n_bits(10);
    /// let check = generator.get_random_check();
    /// ```
    pub fn get_random_check(&mut self) -> Option<Vec<usize>> {
        let mut candidate_check = self.get_candidate_check();
        if self.is_valid_check(&candidate_check) {
            candidate_check.sort();
            self.update_from_check(&candidate_check);
            Some(candidate_check)
        } else {
            None
        }
    }

    fn get_candidate_check(&mut self) -> Vec<usize> {
        let mut check = Vec::with_capacity(self.target_check_degree);
        for _ in 0..self.target_check_degree {
            self.add_random_bit_to_check(&mut check);
        }
        check
    }

    fn add_random_bit_to_check(&mut self, check: &mut Vec<usize>) {
        self.get_random_bit_generator_for_check(check)
            .add_random_bit_to_check(check);
    }

    fn get_random_bit_generator_for_check(&mut self, check: &[usize]) -> RandomBitGenerator<R> {
        let availables = self.get_available_bits_for_check(check);
        let distribution = self.get_distribution_over(&availables);
        RandomBitGenerator {
            availables,
            distribution,
            random_number_generator: &mut self.random_number_generator,
        }
    }

    fn get_available_bits_for_check(&self, check: &[usize]) -> Vec<usize> {
        self.active_bits
            .iter()
            .filter(|bit| self.is_available(**bit))
            .filter(|bit| self.is_not_adjacent_to_check(bit, check))
            .cloned()
            .collect()
    }

    fn is_available(&self, bit: usize) -> bool {
        self.bit_degrees[bit] < self.maximal_bit_degree
    }

    fn is_not_adjacent_to_check(&self, bit: &usize, check: &[usize]) -> bool {
        check.iter().all(|b| self.are_not_adjacent(*b, bit))
    }

    fn are_not_adjacent(&self, bit_0: usize, bit_1: &usize) -> bool {
        !self.get_bits_adjacent_to(bit_0).contains(bit_1)
    }

    fn get_distribution_over(&self, bits: &[usize]) -> Vec<f64> {
        bits.iter().map(|bit| self.distribution[*bit]).collect()
    }

    fn is_valid_check(&self, check: &[usize]) -> bool {
        match self.check_degree_condition {
            CheckDegreeCondition::MustBeAtLeastOfDegree(min_degree) => check.len() >= min_degree,
            CheckDegreeCondition::MustBeFullDegree => check.len() == self.target_check_degree,
        }
    }

    fn update_from_check(&mut self, check: &[usize]) {
        self.update_degrees_from_check(check);
        self.update_adjacency_from_check(check);
    }

    fn update_degrees_from_check(&mut self, check: &[usize]) {
        check.iter().for_each(|bit| self.bit_degrees[*bit] += 1);
    }

    fn update_adjacency_from_check(&mut self, check: &[usize]) {
        self.adjacency.update_from_check(check)
    }
}

enum CheckDegreeCondition {
    MustBeFullDegree,
    MustBeAtLeastOfDegree(usize),
}

struct RandomBitGenerator<'a, R: Rng> {
    availables: Vec<usize>,
    distribution: Vec<f64>,
    random_number_generator: &'a mut R,
}

impl<'a, R: Rng> RandomBitGenerator<'a, R> {
    fn add_random_bit_to_check(self, check: &mut Vec<usize>) -> Self {
        if let Ok(distribution) = WeightedIndex::new(&self.distribution) {
            let sample = self.random_number_generator.sample(distribution);
            check.push(self.availables[sample]);
        }
        self
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn an_empty_generator_does_not_generate_checks() {
        let mut generator = Generator::new();
        assert!(generator
            .set_target_check_degree(2)
            .get_random_check()
            .is_none());
    }

    #[test]
    fn a_generator_does_not_generate_checks_of_degree_less_than_minimal_check_degree() {
        let mut generator = Generator::with_n_bits(5);
        generator.allow_checks_of_degree_at_least(2);
        assert!(generator
            .set_target_check_degree(0)
            .get_random_check()
            .is_none());
        assert!(generator
            .set_target_check_degree(1)
            .get_random_check()
            .is_none());
        assert!(generator
            .set_target_check_degree(2)
            .get_random_check()
            .is_some());
    }

    #[test]
    fn a_generator_does_not_include_the_same_bit_twice_in_a_check() {
        let mut generator = Generator::with_n_bits(3);
        generator.set_maximal_bit_degree(2).set_minimal_girth(4);

        let first_check = generator
            .set_over_bits(vec![0, 1])
            .set_target_check_degree(2)
            .get_random_check();
        assert_eq!(first_check, Some(vec![0, 1]));

        // Can't include both bit 0 and 1 in the next check.
        let second_check = generator.set_target_check_degree(3).get_random_check();
        assert_eq!(second_check, None);
    }

    #[test]
    fn non_full_check_can_be_set_to_be_allowed() {
        let mut generator = Generator::with_n_bits(3);
        generator
            .allow_checks_of_degree_at_least(2)
            .set_maximal_bit_degree(2)
            .set_minimal_girth(4)
            .set_target_check_degree(2);

        // Generate a check over the first 2 bits and one over the 2 last.
        let first_check = generator.set_over_bits(vec![0, 1]).get_random_check();
        let second_check = generator.set_over_bits(vec![1, 2]).get_random_check();
        assert_eq!(first_check, Some(vec![0, 1]));
        assert_eq!(second_check, Some(vec![1, 2]));

        // Bit 2 has degree 2. Thus a check with target degree 3 will be generated over bit 0 and 2
        // only.
        let third_check = generator
            .set_over_all_bits()
            .set_target_check_degree(3)
            .get_random_check();
        assert_eq!(third_check, Some(vec![0, 2]));
    }

    #[test]
    fn a_generator_does_not_exceed_bit_maximal_degree() {
        let mut generator = Generator::with_n_bits(3);
        generator
            .set_maximal_bit_degree(2)
            .set_target_check_degree(2);

        let first_check = generator.set_without_bits(vec![2]).get_random_check();
        assert_eq!(first_check, Some(vec![0, 1]));

        let second_check = generator
            .set_over_all_bits()
            .set_without_bits(vec![0])
            .get_random_check();
        assert_eq!(second_check, Some(vec![1, 2]));

        // We already have checks [0,1] and [1, 2]. Degree of bit 1 is 2 and it can't
        // be included in another check.

        let third_check = generator
            .set_over_all_bits()
            .set_target_check_degree(3)
            .get_random_check();
        assert_eq!(third_check, None);

        let fourth_check = generator.set_target_check_degree(2).get_random_check();
        assert_eq!(fourth_check, Some(vec![0, 2]));

        // Every bit has max degree. Can't generate anymore check
        assert_eq!(generator.get_random_check(), None);
    }

    #[test]
    fn a_generator_does_not_create_cycle_of_length_two_and_four_if_minimal_girth_is_six() {
        let mut generator =
            Generator::with_n_bits(3).with_random_number_generator(ChaCha8Rng::seed_from_u64(10));

        generator.set_maximal_bit_degree(2).set_minimal_girth(6);

        let first_check = generator.set_target_check_degree(3).get_random_check();
        assert_eq!(first_check, Some(vec![0, 1, 2]));

        // Any check of degree 2 will create a 4-cycle.
        let second_check = generator.set_target_check_degree(2).get_random_check();
        assert_eq!(second_check, None);

        // A degree 1 check will not create a 4-cycle.
        let third_check = generator.set_target_check_degree(1).get_random_check();
        assert_eq!(third_check.is_some(), true);
    }

    #[test]
    fn a_generator_does_not_create_cycle_of_length_two_four_and_six_if_minimal_girth_is_eight() {
        let mut generator =
            Generator::with_n_bits(5).with_random_number_generator(ChaCha8Rng::seed_from_u64(10));

        generator.set_maximal_bit_degree(2).set_minimal_girth(8);

        let first_check = generator
            .set_over_bits(vec![0, 1, 2])
            .set_target_check_degree(3)
            .get_random_check();
        assert_eq!(first_check, Some(vec![0, 1, 2]));

        generator.set_target_check_degree(2);
        let second_check = generator.set_over_bits(vec![2, 3]).get_random_check();
        assert_eq!(second_check, Some(vec![2, 3]));

        // A check over [0, 3] will create a 6-cycle.
        let third_check = generator.set_over_bits(vec![0, 3]).get_random_check();
        assert_eq!(third_check, None);

        // Possible checks are [0, 4] or [3, 4]
        let fourth_check = generator.set_over_bits(vec![0, 3, 4]).get_random_check();
        assert_eq!(fourth_check.clone().unwrap().contains(&4), true);
        assert_eq!(fourth_check.unwrap().len(), 2);
    }

    #[test]
    fn a_generator_can_generate_check_according_to_a_bit_distribution() {
        let mut generator =
            Generator::with_n_bits(5).with_random_number_generator(ChaCha8Rng::seed_from_u64(10));

        generator
            .set_maximal_bit_degree(2)
            .set_distribution(vec![0.25, 0.25, 0.0, 0.25, 0.25])
            .set_over_bits(vec![0, 1, 2])
            .set_target_check_degree(3);

        // Can't generate 3 bits from this distribution over the first 3.
        assert_eq!(generator.get_random_check(), None);

        generator.set_target_check_degree(2);
        assert_eq!(generator.get_random_check(), Some(vec![0, 1]));
        assert_eq!(generator.get_random_check(), Some(vec![0, 1]));

        // Degree of the first 2 bits is 2.
        assert_eq!(generator.get_random_check(), None);

        generator.set_over_all_bits();
        assert_eq!(generator.get_random_check(), Some(vec![3, 4]));

        generator.set_target_check_degree(3);

        // Can't pick a degree 3 check because probability of bit 2 is 0.
        assert_eq!(generator.get_random_check(), None);

        // Reset distribution.
        generator.set_uniform_distribution();
        assert_eq!(generator.get_random_check(), Some(vec![2, 3, 4]));
    }
}

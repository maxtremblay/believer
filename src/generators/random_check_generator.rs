/// A `RandomCheckGenerator` helps generating checks for a code while respecting some global 
/// constraints.
/// 
/// Constraints are bit and check degrees and minimal girth.
pub struct RandomCheckGenerator {

}

impl RandomCheckGenerator {
    pub fn generate(&self, check_degree: usize) -> Option<Vec<usize>> {
        unimplemented!()
    }

    pub fn new(bit_degrees: &[usize], minimal_girth: usize) -> Self {
        unimplemented!()
    }

    pub fn over_all_bits(&mut self) -> Self {
        unimplemented!()
    }

    pub fn over_bits(&mut self, bits: &[usize]) -> Self{
        unimplemented!()
    }

    pub fn with_distribution(&mut self, probs: &[f64]) -> Self {
        unimplemented!()
    }

    pub fn with_uniform_distribution(&mut self) -> Self {

    }

    pub fn without(&mut self, bits: &[usize]) -> Self {
        unimplemented!()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn doesnt_include_same_bit_twice() {
        let mut generator = RandomCheckGenerator::new(&[2, 2, 2], 0);
        let first_check = generator.over_bits(&[0, 1]).generate(2);
        assert_eq!(first_check, Some(vec![0, 1]));

        // Can't generate a degree 4 check over only 2 bits.
        let second_check = generator.generate(4);
        assert_eq!(second_check, None);
    }

    #[test]
    fn doesnt_exceed_bit_maximal_degree() {
        let mut generator = RandomCheckGenerator::new(&[2, 2, 2], 0);
        
        let first_check = generator.without(&[2]).generate(2);
        assert_eq!(first_check, Some(vec![0, 1]));

        let second_check = generator.over_all_bits().without(&[0]).generate(2);
        assert_eq!(second_check, Some(vec![1, 2]));

        // We already have checks [0,1] and [1, 2]. Degree of bit 1 is 2 and it can't be include in
        // another check.
        
        let third_check = generator.over_all_bits().generate(3);
        assert_eq!(third_check, None);

        let fourth_check = generator.generate(2);
        assert_eq!(fourth_check, Some(vec![0, 2]));

        // Every bit has max degree. Can't generate anymore check
        assert_eq!(generator.generate(1), None);
    }

    #[test]
    fn doesnt_create_cycle_smaller_than_minimal_girth() {
        // Minimal girth 6
        let generator = RandomCheckGenerator::new(&[2, 2, 2], 6);
        
        let first_check = generator.generate(3);
        assert_eq!(first_check, Some(vec![0, 1, 2]));

        // Any check of degree 2 will create a 4-cycle. 
        let second_check = generator.generate(2);
        assert_eq!(second_check, None);

        // A degree 1 check will not create a 4-cycle.
        let third_check = generator.generate(1);
        assert_eq!(third_check.is_some(), true);

        // Minimal girth 8
        let mut generator = RandomCheckGenerator::new(&[2; 5], 8);

        let first_check = generator.over_bits(&[0, 1, 2]).generate(3);
        assert_eq!(first_check, Some(vec![0, 1, 2]));

        let second_check = generator.over_bits(&[2, 3]).generate(2);
        assert_eq!(second_check, Some(vec![2, 3]));

        // A check over [0, 3] will create a 6-cycle.
        let third_check = generator.over_bits(&[0, 3]).generate(2);
        assert_eq!(third_check, None);

        // Possible checks are [0, 4] or [3, 4]
        let fourth_check = generator.over_bits(&[0, 3, 4]).generate(2);
        assert_eq!(fourth_check.clone().unwrap().contains(&4), true);
        assert_eq!(fourth_check.unwrap().len(), 2);
    }

    #[test]
    fn generate_bit_according_to_distribution() {
        let generator = RandomCheckGenerator::new(&[2; 5], 0)
            .with_distribution(&[0.25, 0.25, 0.0, 0.25, 0.25]);
        
        // Generate all degree 2 checks from the given distribution
        for _ in 0..6 {
            assert_eq!(generator.generate(2).is_some(), true);
        }
        assert_eq!(generator.generate(2), None);

        let mut generator = RandomCheckGenerator::new(&[2; 5], 0)
            .with_distribution(&[0.25, 0.25, 0.0, 0.25, 0.25])
            .over_bits(&[0, 1, 2]);

        // Can't generate 3 bits from this distribution over the first 3.
        assert_eq!(generator.generate(3), None);
        assert_eq!(generator.generate(2), Some(vec![0, 1]));
        assert_eq!(generator.generate(2), Some(vec![0, 1]));
        
        // Degree of the first 2 bits is 2.
        assert_eq!(generator.generate(2), None);

        generator.over_all_bits();
        assert_eq!(generator.generate(2), Some(vec![3, 4]));

        // Can't pick a degree 3 check because probability of bit 2 is 0.
        assert_eq!(generator.generate(3), None);

        // Reset distribution. 
        generator.with_uniform_distribution();
        assert_eq!(generator.generate(3), Some(vec![2, 3, 4]));
    }
}
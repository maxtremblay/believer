//! The binary symmetric channel.

use qecsim::noise::binary::BinaryError;
use qecsim::noise::IndependentNoiseModel;
use qecsim::random::RandomNumberGenerator;
use rand::distributions::{Bernoulli, Distribution};

/// A binary symmetric channel caracterized by an error probability p.
/// Any input is either affected by the identity error with probability (1 - p) or
/// is flipped with probability p.
///
/// To use the full features of this struct, the `IndependentNoiseModel` trait of the crate qecsim
/// must be imported.
///
/// # Example
///
/// ```
/// # use believer::classical::channel::BinarySymmetricChannel;
/// use qecsim::random::rng_with_seed;
/// use qecsim::noise::IndependentNoiseModel;
/// // Create a bsc with an error probability of 0.2.
/// let bsc = BinarySymmetricChannel::with_probability(0.2);
/// let error_block = bsc.sample_error_block(1000, &mut rng_with_seed(72));
/// println!("{}", error_block.weight()); // Should be around 200.
/// ```
#[derive(Clone, Debug)]
pub struct BinarySymmetricChannel {
    distribution: Bernoulli,
    log_likelyhood: f64,
}

impl BinarySymmetricChannel {
    /// Creates a new binary symmetric channel from an given `probability` of error.
    ///
    /// # Panic
    ///
    /// Panic if `probability` is not between 0 and 1.
    pub fn with_probability(probability: f64) -> Self {
        let distribution = Bernoulli::new(probability).expect("probability is not between 0 and 1");
        Self {
            distribution,
            log_likelyhood: (probability / (1.0 - probability)).log2(),
        }
    }

    /// Computes the instrinsic likelyhood of the given `error`.
    ///
    /// If the probability of flipping a bit is p, then the intrinsic likelyhood of identity is
    /// log( p / (1 - p) ) and the intrinsic likelyhood of flip is log( (1 - p) / p ).
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::classical::channel::BinarySymmetricChannel;
    /// use qecsim::noise::binary::BinaryError;
    /// let bsc = BinarySymmetricChannel::with_probability(0.2);
    /// assert_eq!(bsc.intrinsic_likelyhood_of(BinaryError::Identity), -2.0);
    /// assert_eq!(bsc.intrinsic_likelyhood_of(BinaryError::Flipped), 2.0);
    /// ```
    pub fn intrinsic_likelyhood_of(&self, error: BinaryError) -> f64 {
        match error {
            BinaryError::Identity => self.log_likelyhood,
            BinaryError::Flipped => -1.0 * self.log_likelyhood,
        }
    }
}

impl IndependentNoiseModel for BinarySymmetricChannel {
    type Error = BinaryError;

    fn sample(&self, rng: &mut RandomNumberGenerator) -> Self::Error {
        if self.distribution.sample(rng) {
            BinaryError::Flipped
        } else {
            BinaryError::Identity
        }
    }
}

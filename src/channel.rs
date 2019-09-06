//! Transmission channels.

use crate::GF2;
use rand::distributions::Uniform;
use rand::{thread_rng, Rng};

/// A trait that represent a binary transmission channel. It takes an `GF2` element
/// and (randomly) map it to an other `GF2` element. It has an intrinsic likelyhood
/// for each output.
///
/// # Example
///
/// ```
/// # use believer::*;
/// // Create a bsc with error prob of 0.2.
/// let bsc = channel::BinarySymmetricChannel::new(0.2);
/// // Sample the channel by always sending 0.
/// let received = bsc.sample_uniform(GF2::B0, 1000);
/// let number_of_one = received.iter()
///     .filter(|&x| x == &GF2::B1)
///     .collect::<Vec<_>>()
///     .len();
/// println!("{}", number_of_one); // Should be around 200.
/// ```
pub trait BinaryChannel: Sync {
    /// For a given `output`, compute log(p(output|input = 1) / p(output|input = 0)).
    ///
    /// # Example
    ///
    /// ```
    /// # use believer::*;
    /// let bsc = channel::BinarySymmetricChannel::new(0.2);
    /// assert_eq!(bsc.intrinsic_likelyhood(GF2::B0), -2.0);
    /// assert_eq!(bsc.intrinsic_likelyhood(GF2::B1), 2.0);
    /// ```
    fn intrinsic_likelyhood(&self, output: GF2) -> f64;

    /// For a given `input`, returns an `output` according to some
    /// probability distribution depending on the channel.
    fn send(&self, input: GF2) -> GF2;

    // Returns the likelyhood of a given message.
    fn message_likelyhood(&self, output: &[GF2]) -> Vec<f64> {
        output
            .iter()
            .map(|x| self.intrinsic_likelyhood(*x))
            .collect()
    }

    /// Returns a `Vec` of outputs as long as the `inputs` where each
    /// output is sample using the `send` method.
    fn sample(&self, inputs: &[GF2]) -> Vec<GF2> {
        inputs.iter().map(|input| self.send(*input)).collect()
    }

    /// Returns a `Vec` of outputs where the `input` is send `n_inputs` times.
    fn sample_uniform(&self, input: GF2, n_inputs: usize) -> Vec<GF2> {
        (0..n_inputs).map(|_| self.send(input)).collect()
    }
}

/// A binary symmetric channel caracterize by its error probability `prob`.
/// That is, every time an input is send throught the channel, it is
/// flipped with probability `prob`.
pub struct BinarySymmetricChannel {
    prob: f64,
    log_likelyhood: f64,
}

impl BinarySymmetricChannel {
    /// Creates a new binary symmetric channel from an given `prob` of error.
    ///
    /// # Panic
    ///
    /// Panic if `prob` is not between 0 and 1.
    pub fn new(prob: f64) -> Self {
        if 0.0 <= prob && prob <= 1.0 {
            Self {
                prob,
                log_likelyhood: (prob / (1.0 - prob)).log2(),
            }
        } else {
            panic!("prob is not between 0 and 1")
        }
    }
}

impl BinaryChannel for BinarySymmetricChannel {
    fn intrinsic_likelyhood(&self, output: GF2) -> f64 {
        if output == GF2::B0 {
            self.log_likelyhood
        } else {
            -1.0 * self.log_likelyhood
        }
    }

    fn send(&self, input: GF2) -> GF2 {
        let rand = thread_rng().sample(Uniform::new(0.0, 1.0));
        if rand < self.prob {
            input + GF2::B1
        } else {
            input
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn binary_symmetric_channel() {
        let channel = BinarySymmetricChannel::new(0.2);

        assert_eq!(channel.intrinsic_likelyhood(GF2::B0), -2.0);
        assert_eq!(channel.intrinsic_likelyhood(GF2::B1), 2.0);
        assert_eq!(
            channel.message_likelyhood(&[GF2::B1, GF2::B0, GF2::B1]),
            vec![2.0, -2.0, 2.0]
        );
    }
}

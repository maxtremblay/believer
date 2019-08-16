use crate::GF2;
use rand::{Rng, thread_rng};
use rand::distributions::Uniform;

pub trait Channel {
    type Input: Copy;
    type Output;

    fn intrinsic_likelyhood(&self, output: Self::Output) -> f64;
    fn send(&self, input: Self::Input) -> Self::Output;
    
    fn sample(&self, inputs: &[Self::Input]) -> Vec<Self::Output> {
        inputs.iter().map(|input| self.send(*input)).collect()
    }

    fn sample_uniform(&self, input: Self::Input, n_inputs: u32) -> Vec<Self::Output> {
        (0..n_inputs).map(|_| self.send(input)).collect()
    }
}

pub struct BinarySymmetricChannel {
    prob: f64,
    log_likelyhood: f64,
}

impl BinarySymmetricChannel {
    pub fn new(prob: f64) -> Result<Self, &'static str> {
        if 0.0 <= prob && prob <= 1.0 {
            Ok(Self {
                prob,
                log_likelyhood: ((1.0 - prob) / prob).log2()
            })
        } else {
            Err("prob is not between 0 and 1")
        }
    }
}

impl Channel for BinarySymmetricChannel {
    type Input = GF2;
    type Output = GF2;

    fn send(&self, input: Self::Input) -> Self::Output {
        let rand = thread_rng().sample(Uniform::new(0.0, 1.0));
        if rand < self.prob {
            input + GF2::B1
        } else {
            input
        }
    }

    fn intrinsic_likelyhood(&self, output: Self::Output) -> f64 {
        if output == GF2::B0 {
            self.log_likelyhood
        } else {
            -1.0 * self.log_likelyhood
        }
    }
}
// Here we define the random number generators to use in the crate. We define 2 rngs, the
// one when no seed is provided and the seeded one when a u64 seed is provided. 
//
// IMPORTANT: Changing those will affecte the reproductability of the results.

use rand::{SeedableRng, thread_rng};
use rand::rngs::ThreadRng;
use rand_chacha::ChaCha8Rng;

// The random number generator when no seed are provided.
pub(crate) fn get_random_rng_without_seed() -> ThreadRng {
    thread_rng()
}

// The random number generator when a u64 seed is provided.
pub(crate) fn get_rng_with_seed(seed: u64) -> ChaCha8Rng {
    ChaCha8Rng::seed_from_u64(seed)
}

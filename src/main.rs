use believer::{BinaryChannel, Decoder, GF2, ParityCheckMatrix};
use believer::channel::BinarySymmetricChannel;
use rayon::prelude::*;
use std::fs;
// use std::io::prelude::*;

fn generate_parity_check(m: usize) -> ParityCheckMatrix {
    let mut checks = Vec::with_capacity(4 * m);
    for i in 0..(3 * m) {
        checks.push(vec![2 * i, 2 * i + 1, (2 * i + 2) % (6 * m)]);
    }

    for i in 0..m {
        checks.push(vec![2 * i + 1, 2 * i + 1 + 2 * m, 2 * i + 1 + 4 * m])
    }
    ParityCheckMatrix::new(checks)
}

fn write_vec_to_file<T: ToString>(vec: Vec<T>, file: &str) {
    let data = vec.iter()
        .fold(String::new(), |mut acc, x| {
            acc.push_str(&x.to_string());
            acc.push('\n');
            acc
        });
    fs::write(file, data).expect("Unable to write file");
}

fn main() {
    // *
    // Parameters
    // * 

    let min_m = 1;
    let max_m = 10;
    let n_errors_per_thread = 2;
    let n_thread = 2;
    let error_prob = 0.25;
    let max_iters = 50;

    // *
    // Setup
    // *

    let channel = BinarySymmetricChannel::new(error_prob);
    
    // *
    // Main loop
    // *

    let error_rates: Vec<f64> = (min_m..=max_m).map(|m| {
        let parity_check = generate_parity_check(m);
        let decoder = Decoder::new(&channel, &parity_check);

        (0..n_thread).into_par_iter()
            .map(|_| {
                let mut errors = 0;
                let mut total = 0;
                while errors < n_errors_per_thread {
                    total += 1; 
                    let received_message = channel.sample_uniform(GF2::B0, 6 * m);
                    let decoded_message = decoder.decode(&received_message, max_iters);
                    if decoded_message != Some(vec![GF2::B0; 6 * m]) {
                        errors += 1;
                    }
                }
                total
            })
            .sum()  
    })
    .map(|total: u32| (n_thread * n_errors_per_thread) as f64 / total as f64)
    .collect();

    write_vec_to_file(error_rates, "error_rates.txt");
    write_vec_to_file((min_m..=max_m).collect(), "repetitions.txt");
}





use believer::{BinaryChannel, Decoder, GF2, ParityCheckMatrix};
use believer::channel::BinarySymmetricChannel;
use rayon::prelude::*;
use std::fs;

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

// fn write_vec_to_file<T: ToString>(vec: &[T], file: &str) {
//     let data = vec.iter()
//         .fold(String::new(), |mut acc, x| {
//             acc.push_str(&x.to_string());
//             acc.push('\n');
//             acc
//         });
//     fs::write(file, data).expect("Unable to write file");
// }

fn write_data(total: usize, wrong_codeword: usize, stuck: usize, file: &str) {
    let to_write = format!("{} : {} / {}", wrong_codeword, stuck, total);
    fs::write(file, to_write).expect("Unable to write file");
}

fn main() {
    // *
    // Parameters
    // * 

    let min_m: usize = 17;
    let max_m: usize = 17;
    let n_errors_per_thread = 5;
    let n_threads = 4;
    let error_prob = 0.25;
    let max_iters = 200;

    // *
    // Setup
    // *

    let channel = BinarySymmetricChannel::new(error_prob);

    let directory = format!(
        "total_iters_prob_{}",
        error_prob
    );
    
    fs::create_dir(&directory).expect("Unable to create directory");

    // *
    // Main loop
    // *
   
    (min_m..=max_m).into_par_iter().for_each(|m: usize| {
        let parity_check = generate_parity_check(m);
        let decoder = Decoder::new(&channel, &parity_check);
        
        let local_directory = format!("{}/m_{}", directory, m);
        fs::create_dir(&local_directory).expect("Unable to create directory");

        (0..n_threads)
            .into_par_iter()
            .for_each(|idx| {
                let mut wrong_codeword = 0;
                let mut stuck = 0;
                let mut total = 0;
                while wrong_codeword + stuck < n_errors_per_thread {
                    total += 1; 
                    let received_message = channel.sample_uniform(GF2::B0, 6 * m);
                    let decoded_message = decoder.decode(&received_message, max_iters);
                    if let Some(b) = decoded_message {
                        if b != vec![GF2::B0; 6 * m] {
                            wrong_codeword += 1;
                        }
                    } else {
                        stuck += 1;
                    }
                }
                let file = format!("{}/{}.txt", local_directory, idx);
                write_data(total, wrong_codeword, stuck, &file);
            });
    });

}

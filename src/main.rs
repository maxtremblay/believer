use believer::channel::BinarySymmetricChannel;
use believer::{BinaryChannel, Decoder, ParityCheckMatrix, GF2, DecodingResult};
use rayon::prelude::*;
use std::fs;

// fn generate_parity_check(m: usize) -> ParityCheckMatrix {
//     let mut checks = Vec::with_capacity(4 * m);
//     for i in 0..(3 * m) {
//         checks.push(vec![2 * i, 2 * i + 1, (2 * i + 2) % (6 * m)]);
//     }

//     for i in 0..m {
//         checks.push(vec![2 * i + 1, 2 * i + 1 + 2 * m, 2 * i + 1 + 4 * m])
//     }

//     for i in 0..(3 * m / 2) {
//         checks.push(vec![4 * i + 1, 4 * i + 3, (4 * i + 5) % (6 * m)]);
//     }                   
    
//     ParityCheckMatrix::new(checks)
// }

fn generate_parity_check(m: usize) -> ParityCheckMatrix {
    let mut checks = Vec::with_capacity(11 * m);
    for i in 0..(4 * m) {
        checks.push(vec![3 * i, 3 * i + 1, 3 * i + 2, (3 * i + 3) % (12 * m)]);
        checks.push(vec![3 * i, (3 * i + 4) % (12 * m), (3 * i + 8) % (12 * m), 12 * m + i]);
    }
    for i in 0..m {
        checks.push(vec![1 + 3 * i, 1 + 3 * i + 3 * m, 1 + 3 * i + 6 * m, 1 + 3 * i + 6 * m]);
        checks.push(vec![2 + 3 * i, 2 + 3 * i + 3 * m, 2 + 3 * i + 6 * m, 2 + 3 * i + 6 * m]);
        checks.push(vec![12 * m + i , 13 * m + i, 14 * m + i, 15 * m + i]);
    }
    // println!("CHECKS: {:#?}", checks);
    ParityCheckMatrix::new(checks)
}

// fn rep_code(m: usize) -> ParityCheckMatrix {
//     let mut checks = Vec::with_capacity(m);
//     for i in 0..m {
//         checks.push(vec![i, (i + 1) % m]);
//     }
//     ParityCheckMatrix::new(checks)
// }

// fn write_vec_to_file<T: ToString>(vec: &[T], file: &str) {
//     let data = vec.iter()
//         .fold(String::new(), |mut acc, x| {
//             acc.push_str(&x.to_string());
//             acc.push('\n');
//             acc
//         });
//     fs::write(file, data).expect("Unable to write file");
// }

fn write_data(
    total: usize,
    error_type: &[char],
    weights: &[usize],
    inputs: &[Vec<GF2>],
    outputs: &[Vec<GF2>],
    file: &str,
) {
    let mut to_write = format!("Total: {}\n", total);
    to_write.push_str(&format!("Number of errors: {}\n\n", error_type.len()));
    error_type
        .iter()
        .zip(weights.iter())
        .zip(inputs.iter())
        .zip(outputs.iter())
        .for_each(|(((t, w), _), _)| {
            to_write.push_str(&format!("{}, {}\n", t, w));
        });
    fs::write(file, to_write).expect("Unable to write file");
}

fn main() {
    // *
    // Parameters
    // *

    let min_m: usize = 1;
    let max_m: usize = 10;
    let n_errors_per_thread = 5;
    let n_threads = 2;
    let error_prob = 0.01;
    let max_iters = 50;

    // *
    // Setup
    // *

    let channel = BinarySymmetricChannel::new(error_prob);

    let directory = format!("total_iters_prob_{}", error_prob);

    fs::create_dir(&directory).expect("Unable to create directory");

    // *
    // Main loop
    // *

    (min_m..=max_m).into_par_iter().map(|m| 10 * m).for_each(|m: usize| {
        let parity_check = generate_parity_check(m);
        // println!("N BITS: {}", parity_check.n_bits());
        // println!("M: {}", m);

        let decoder = Decoder::new(&channel, &parity_check);

        let local_directory = format!("{}/m_{}", directory, m);
        fs::create_dir(&local_directory).expect("Unable to create directory");

        (0..n_threads).into_par_iter().for_each(|idx| {
            let mut error_type = Vec::with_capacity(n_errors_per_thread);
            let mut total = 0;
            let mut weights: Vec<usize> = Vec::with_capacity(n_errors_per_thread);
            let mut inputs = Vec::with_capacity(n_errors_per_thread);
            let mut outputs = Vec::with_capacity(n_errors_per_thread);
            while error_type.len() < n_errors_per_thread {
                total += 1;
                let received_message = channel.sample_uniform(GF2::B0, 16 * m);
                let decoded_message = decoder.decode(&received_message, max_iters);
                if let DecodingResult::Codeword(b) = decoded_message {
                    if b != vec![GF2::B0; 16 * m] {
                        error_type.push('W');
                        weights.push(
                            received_message
                                .iter()
                                .map(|x| if x == &GF2::B0 { 0 } else { 1 })
                                .sum(),
                        );
                        inputs.push(received_message);
                        outputs.push(b);
                    }
                } else {
                    error_type.push('S');
                    weights.push(
                        received_message
                            .iter()
                            .map(|x| if x == &GF2::B0 { 0 } else { 1 })
                            .sum(),
                    );
                    inputs.push(received_message);
                    outputs.push(Vec::new());
                }
            }
            let file = format!("{}/{}.txt", local_directory, idx);
            write_data(total, &error_type, &weights, &inputs, &outputs, &file);
            // write_vec_to_file(&error_type, &format!("{}/{}_errors.txt", local_directory, idx));
            // write_vec_to_file(&weights, &format!("{}/{}_weights.txt", local_directory, idx));
        });
    });
}

// fn main() {
//     let channel = BinarySymmetricChannel::new(0.1);
//     println!("{:#?}", channel.sample_uniform(GF2::B0, 10));
// }

// use believer::*;

// fn main() {
//     let channel = channel::BinarySymmetricChannel::new(0.2);
//     let parity_check = ParityCheckMatrix::new(vec![
//         vec![0, 1],
//         vec![1, 2],
//         vec![2, 3],
//     ]);
//     let decoder = Decoder::new(&channel, &parity_check);

//     // Should be able to decode this message to the 1 codeword.
//     let easy_message = vec![GF2::B0, GF2::B1, GF2::B1, GF2::B1];
//     assert_eq!(decoder.decode(&easy_message, 10), Some(vec![GF2::B1; 4]));

//     // Should get stuck while decoding this message.
//     let impossible_message = vec![GF2::B0, GF2::B1, GF2::B0, GF2::B1];
//     assert_eq!(decoder.decode(&impossible_message, 10), None);
// }

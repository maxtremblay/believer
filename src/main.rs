use believer::channel::BinarySymmetricChannel;
use believer::{BinaryChannel, Decoder, ParityCheckMatrix, GF2};
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
        .for_each(|(((t, w), i), o)| {
            to_write.push_str(&format!("{}, {}:\n{:?}\n{:?}\n\n", t, w, i, o));
        });
    fs::write(file, to_write).expect("Unable to write file");
}

fn main() {
    // *
    // Parameters
    // *

    let min_m: usize = 5;
    let max_m: usize = 5;
    let n_errors_per_thread = 1;
    let n_threads = 1;
    let error_prob = 0.05;
    let max_iters = 1000;

    // *
    // Setup
    // *

    let channel = BinarySymmetricChannel::new(error_prob);

    let directory = format!("total_iters_prob_{}", error_prob);

    fs::create_dir(&directory).expect("Unable to create directory");

    // *
    // Main loop
    // *

    (min_m..=max_m).into_par_iter().for_each(|m: usize| {
        let parity_check = generate_parity_check(m);
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
                // let received_message = channel.sample_uniform(GF2::B0, 6 * m);
                let mut received_message = vec![GF2::B0; 23];
                received_message.append(&mut vec![GF2::B1; 3]);
                received_message.append(&mut vec![GF2::B0; 4]);
                let decoded_message = decoder.decode(&received_message, max_iters);
                if let Some(b) = decoded_message {
                    if b != vec![GF2::B0; 6 * m] {
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

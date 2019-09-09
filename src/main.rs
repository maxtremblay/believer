use believer::*;

fn main() {
    let parity_check = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2]]);

    let decoder = ErasureDecoder::new(&parity_check);

    println!("{:?}", decoder.decode(&[1]));
}

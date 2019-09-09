use believer::*;

fn main() {
    let checks = ParityCheckMatrix::new(vec![vec![0, 1], vec![1, 2], vec![2, 3], vec![3, 4]]);
    let decoder = ErasureDecoder::new(&checks, 0.25);
    let simulator = Simulator::new(&decoder);
    let result = simulator.simulate_until_failures_are_found(1, 1);
    println!("Failure rate: {}", result.total());
}

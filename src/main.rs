use believer::*;

fn main() {
    let parity_check = ParityCheckMatrix::new(vec![
        vec![0, 1],
        vec![1, 2],
        vec![2, 3],
        vec![3, 4],
    ]);
    let simulator = ErasureSimulator::new(&parity_check, 1, 1, 0.5);
    let result = simulator.simulate();
    println!("Failure rate: {}", result.failure_rate());
}


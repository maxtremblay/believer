use super::ParityCheckMatrix;

pub struct HyperGraphProduct<'a> {
    first_code: &'a ParityCheckMatrix,
    second_code: &'a ParityCheckMatrix,
}

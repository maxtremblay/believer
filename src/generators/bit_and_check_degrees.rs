#[derive(Clone, Copy)]
pub struct BitAndCheckDegrees {
    bit_degree: u32,
    check_degree: u32,
}

impl BitAndCheckDegrees {
    pub fn new(bit_degree: u32, check_degree: u32) -> Self {
        Self { bit_degree, check_degree }
    }

    pub fn get_bit_degrees(&self) -> u32 {
        self.bit_degree
    }

    pub fn get_check_degrees(&self) -> u32 {
        self.check_degree
    }

    pub fn will_generate_an_empty_code(&self) -> bool {
        self.bit_degree == 0 || self.check_degree <= 1
    }
}
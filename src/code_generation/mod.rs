use crate::ParityCheckMatrix;

pub trait CodeGenerator {
    fn generate(&self) -> ParityCheckMatrix;
}

pub mod hirarchical_code;

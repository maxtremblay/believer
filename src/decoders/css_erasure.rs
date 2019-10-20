use super::{Checks, Decoder, DecodingResult, DecoderBuilder};
use crate::ParityCheckMatrix;
use rand::{thread_rng, Rng};

pub struct CSSErasureDecoder{
    checks: (ParityCheckMatrix,ParityCheckMatrix),
    erasure_prob: f64,
}

impl CSSErasureDecoder {
    pub fn new(checks:  (ParityCheckMatrix,ParityCheckMatrix), erasure_prob: f64) -> Self {
        if erasure_prob < 0.0 || erasure_prob > 1.0 {
            panic!("invalid probability");
        }
        Self {
            checks,
            erasure_prob,
        }
    }
}

impl Decoder for CSSErasureDecoder{
    type Error = Vec<usize>;
    type Result = CSSErasureResult;

    

    fn decode(&self, error: &Self::Error) -> Self::Result {
        let erased_parity_x_check = self.checks.0.keep(error);
        let erased_parity_z_check = self.checks.1.keep(error);

        let lost_z_info = (error.len() - erased_parity_x_check.rank()) != 0;
        let lost_x_info = (error.len() - erased_parity_z_check.rank()) != 0;

        if lost_x_info && lost_z_info {
            CSSErasureResult::Failed(PauliErasure::XandZ)
        } else if lost_x_info {
            CSSErasureResult::Failed(PauliErasure::X)
        } else if lost_z_info {
            CSSErasureResult::Failed(PauliErasure::Z)
        } else {
            CSSErasureResult::Succeed
        }
    }

    fn random_error(&self) -> Vec<usize> {
        let mut rng = thread_rng();
        (0..self.checks.0.n_bits())
            .filter(|_| rng.gen::<f64>() < self.erasure_prob)
            .collect()
    }
}

pub struct CSSEDBuilder{
    erasure_prob: f64
}

impl CSSEDBuilder{
    pub fn new(erasure_prob: f64) -> Self {
        Self {
            erasure_prob
            }
    }
}

impl DecoderBuilder< (ParityCheckMatrix,ParityCheckMatrix), CSSErasureDecoder> for CSSEDBuilder where
{

    fn from_code(&self,code: (ParityCheckMatrix,ParityCheckMatrix)) -> CSSErasureDecoder { //ici ils demandent explicitement une lifetime
       CSSErasureDecoder::new(code, self.erasure_prob)
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum CSSErasureResult {
    Failed(PauliErasure),
    Succeed,
}

impl DecodingResult for CSSErasureResult {
    fn succeed(&self) -> bool {
        match self {
            Self::Succeed => true,
            _ => false 
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum PauliErasure {
    X,
    Z,
    XandZ
}
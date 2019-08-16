use crate::channel::BinaryChannel;
use crate::sparse_matrix::{SparseMatrix, Transposer};
use crate::ParityCheckMatrix;
use crate::GF2;

pub struct Decoder<'a> {
    channel: &'a BinaryChannel,
    parity_check: &'a ParityCheckMatrix,
    transposer: Transposer,
}

impl<'a> Decoder<'a> {
    pub fn decode(&self, message: Vec<GF2>, max_iters: u32) -> Option<Vec<GF2>> {
        let mut extrinsec_likelyhoods = 
            SparseMatrix::from_parity_check(self.parity_check, vec![0.0; self.parity_check.len()]);
        let intrinsec_likelyhoods = self.channel.message_likelyhood(&message);
        let mut total_likelyhoods = intrinsec_likelyhoods.clone();
        let mut iter = 0;
        let mut decoded_message = message.clone();

        while !self.is_codeword(&decoded_message) && iter < max_iters {
            iter += 1;

            extrinsec_likelyhoods =
                self.check_node_update(&extrinsec_likelyhoods, &total_likelyhoods);

            total_likelyhoods =
                self.bit_node_update(&intrinsec_likelyhoods, &extrinsec_likelyhoods);

            decoded_message = total_likelyhoods
                .iter()
                .map(|l| if l > &0.0 { GF2::B1 } else { GF2::B0 })
                .collect();
        }

        if self.is_codeword(&decoded_message) {
            Some(decoded_message)
        } else {
            None
        }
    }

    pub fn new(channel: &'a BinaryChannel, parity_check: &'a ParityCheckMatrix) -> Self {
        Self {
            channel,
            parity_check,
            transposer: Transposer::new(parity_check),
        }
    }

    fn bit_node_update(
        &self,
        intrinsec_likelyhoods: &[f64],
        extrinsec_likelyhoods: &SparseMatrix,
    ) -> Vec<f64> {
        self.transposer
            .transpose(extrinsec_likelyhoods)
            .rows_iter()
            .map(|row| row.map(|(val, _)| val).sum())
            .zip(intrinsec_likelyhoods)
            .map(|(ext, int): (f64, &f64)| ext + *int)
            .collect()
    }

    fn check_node_update(
        &self,
        extrinsec_likelyhoods: &SparseMatrix,
        total_likelyhoods: &[f64],
    ) -> SparseMatrix<'a> {
        let updated_values = extrinsec_likelyhoods
            .rows_iter()
            .map(|row| {
                row.map(|(val, pos): (f64, usize)| {
                    ((val - total_likelyhoods[pos]) / 2.0)
                })
                .product::<f64>()
                .tanh()
            })
            .map(|x: f64| -2.0 * x.atanh())
            .collect();
        SparseMatrix::from_parity_check(self.parity_check, updated_values)
    }

    // NOTE : Maybe implement parity_check.is_codeword(...).
    // It could be faster if the first 1 is near the beggining of message.
    fn is_codeword(&self, message: &[GF2]) -> bool {
        self.parity_check.dot(message).iter().all(|x| x == &GF2::B0)
    }
}

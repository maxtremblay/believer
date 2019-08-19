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
    pub fn decode(&self, message: &[GF2], max_iters: u32) -> Option<Vec<GF2>> {
        let mut extrinsec_likelyhoods = 
            SparseMatrix::from_parity_check(self.parity_check, vec![0.0; self.parity_check.len()]);
        let intrinsec_likelyhoods = self.channel.message_likelyhood(&message);
        let mut total_likelyhoods = intrinsec_likelyhoods.clone();
        let mut iter = 0;
        let mut decoded_message = message.to_vec();

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
        let likelyhood_diff_products: Vec<f64> = extrinsec_likelyhoods
            .rows_iter()
            .map(|row| {
                row.map(|(val, pos): (&f64, &usize)| {
                    ((*val - total_likelyhoods[*pos]) / 2.0)
                })
                .product::<f64>()
                .tanh()
            })
            .collect();
        
        let updated_values = self.parity_check
            .positions_iter()
            .map(|(row, col)| {
                extrinsec_likelyhoods
                    .get(row, col)
                    .map(|val| {
                        ((val - total_likelyhoods[col]) / 2.0).tanh()
                    })
                    .map(|denominator| {
                        -2.0 * (likelyhood_diff_products[row] / denominator).atanh()
                    })
                    .unwrap_or(0.0)
            })
            .collect();

        SparseMatrix::from_parity_check(self.parity_check, updated_values)
    }

    // NOTE : Maybe implement parity_check.is_codeword(...).
    // It could be faster if the first 1 is near the beggining of message.
    fn is_codeword(&self, message: &[GF2]) -> bool {
        self.parity_check.dot(message).iter().all(|x| x == &GF2::B0)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::channel::BinarySymmetricChannel;

    #[test]
    fn general_usage() {
        // Tests with a 3 bits repetition code over a 
        // bsc with error probability of 0.2.s
        let channel = BinarySymmetricChannel::new(0.2);
        let parity_check = ParityCheckMatrix::new(vec![
            vec![0, 1],
            vec![1, 2],
        ]);
        let decoder = Decoder::new(&channel, &parity_check);

        let decoded_message = decoder.decode(&vec![GF2::B0, GF2::B0, GF2::B1], 10);
        assert_eq!(decoded_message, Some(vec![GF2::B0; 3]));
    }
}

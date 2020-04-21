//! Note: Prototype idea
use crate::ParityCheckMatrix;
use crate::Pauli;

pub struct GF4Stabilizers {
    x_checks: ParityCheckMatrix,
    z_checks: ParityCheckMatrix,
}

impl GF4Stabilizers {

    pub fn empty_with_n_bits(n_bits: usize) -> Self {
        Self {
            x_checks: ParityCheckMatrix::with_n_bits(n_bits),
            z_checks: ParityCheckMatrix::with_n_bits(n_bits)
        }
    }


    pub fn from_dense_paulis(stabilizers: Vec<Vec<Pauli>>, n_qubits: usize) -> Self {
        let mut x_checks = Vec::with_capacity(stabilizers.len());
        let mut z_checks = Vec::with_capacity(stabilizers.len());
        stabilizers.iter().for_each(|stab| {
            let mut x_check = Vec::new();
            let mut z_check = Vec::new();
            stab.iter().enumerate().for_each(|(index, p)| {
                let gf4_representation = p.as_gf4();
                if gf4_representation.0 == 1 {
                    x_check.push(index);
                }
                if gf4_representation.1 == 1 {
                    z_check.push(index);
                }
            });
            x_checks.push(x_check);
            z_checks.push(z_check);
            
        });
        println!("x checks:{:?}", x_checks);
        println!("z checks:{:?}",z_checks);
        Self {
            x_checks: ParityCheckMatrix::with_n_bits(n_qubits).with_checks(x_checks),
            z_checks: ParityCheckMatrix::with_n_bits(n_qubits).with_checks(z_checks)
        }
    }

    pub fn from_sparse_paulis(stabilizers: Vec<Vec<(Pauli, usize)>>, n_qubits: usize) -> Self {
        let mut x_checks = Vec::with_capacity(stabilizers.len());
        let mut z_checks = Vec::with_capacity(stabilizers.len());
        stabilizers.iter().for_each(|stab| {
            let mut x_check = Vec::new();
            let mut z_check = Vec::new();
            stab.iter().for_each(|(pauli, position)| {
                let gf4_representation = pauli.as_gf4();
                if gf4_representation.0 == 1 {
                    x_check.push(*position);
                }
                if gf4_representation.1 == 1 {
                    z_check.push(*position);
                }
            });
            x_checks.push(x_check);
            z_checks.push(z_check);
            
        });
        Self {
            x_checks: ParityCheckMatrix::with_n_bits(n_qubits).with_checks(x_checks),
            z_checks: ParityCheckMatrix::with_n_bits(n_qubits).with_checks(z_checks)
        }
    }

    pub fn from_parity_check_matrices(
        x_checks: ParityCheckMatrix,
        z_checks: ParityCheckMatrix,
    ) -> Self {
        if x_checks.get_n_bits() != z_checks.get_n_bits() {
            panic!("different number of bits in each parity check matrix")
        }
        Self { x_checks, z_checks }
    }

    pub fn x_checks(&self) -> &ParityCheckMatrix {
        &self.x_checks
    }

    pub fn z_checks(&self) -> &ParityCheckMatrix {
        &self.z_checks
    }

    pub fn keep(&self, qubits: &[usize]) -> Self {
        Self {
            x_checks: self.x_checks.keep(qubits),
            z_checks: self.z_checks.keep(qubits),
        }
    }

    pub fn merge(&self) -> ParityCheckMatrix {
        self.x_checks.get_horizontal_concat_with(&self.z_checks)
    }

    pub fn n_qubits(&self) -> usize {
        self.x_checks().get_n_bits()
    }

    pub fn n_stabilizers(&self) -> usize {
        std::cmp::max(
            self.x_checks().get_n_checks(),
            self.z_checks().get_n_checks(),
        )
    }

    pub fn without(&self, qubits: &[usize]) -> Self {
        Self {
            x_checks: self.x_checks.without(qubits),
            z_checks: self.z_checks.without(qubits),
        }
    }
}

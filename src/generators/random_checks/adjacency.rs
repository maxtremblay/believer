use std::collections::{BTreeMap, BTreeSet};

// A tool for random check generator that prevent cycle smaller than the minimal girth during code
// construction.
//
// It keeps track of the adjacencies of every bit. Two bits are adjacent if creating a new check
// containing both of them would create a cycle smaller than the minimal girth.
pub(super) struct Adjacency {
    adjacencies: Vec<BTreeSet<usize>>,
    recursion_depth: usize,
}

impl Adjacency {
    // ***** Construction *****

    pub(super) fn new() -> Self {
        Self {
            adjacencies: Vec::new(),
            recursion_depth: 0,
        }
    }

    pub(super) fn with_n_bits(n_bits: usize) -> Self {
        Self {
            adjacencies: Self::initialize_adjacencies_with_n_bits(n_bits),
            recursion_depth: 0,
        }
    }

    fn initialize_adjacencies_with_n_bits(n_bits: usize) -> Vec<BTreeSet<usize>> {
        (0..n_bits)
            .map(|b| {
                let mut set = BTreeSet::new();
                set.insert(b);
                set
            })
            .collect()
    }

    // ***** Set recursion depth *****

    pub(super) fn set_recursion_depth(&mut self, depth: usize) {
        self.recursion_depth = depth;
    }

    // ***** Get adjacent bits *****

    pub(super) fn get_bits_adjacent_to(&self, bit: usize) -> Vec<usize> {
        if self.is_out_of_bound(bit) {
            Vec::new()
        } else {
            self.get_bits_adjacent_to_inbound_bit(bit)
        }
    }

    fn is_out_of_bound(&self, bit: usize) -> bool {
        bit >= self.get_n_bits()
    }

    fn get_n_bits(&self) -> usize {
        self.adjacencies.len()
    }

    fn get_bits_adjacent_to_inbound_bit(&self, bit: usize) -> Vec<usize> {
        if self.recursion_depth == 0 {
            vec![bit]
        } else {
            let mut getter = self.initialize_adjacent_bits_getter(bit);
            getter.get_adjacent_bits_with_depth(self.recursion_depth)
        }
    }

    fn initialize_adjacent_bits_getter(&self, source_bit: usize) -> AdjacentBitsGetter {
        AdjacentBitsGetter {
            source_bit,
            adjacent_bits: BTreeMap::new(),
            adjacencies: &self.adjacencies,
        }
    }

    // ***** Update adjacent bits *****

    pub(super) fn update_from_check(&mut self, check: &[usize]) {
        check
            .iter()
            .for_each(|bit| self.set_bits_in_check_adjacent_to(*bit, check));
    }

    fn set_bits_in_check_adjacent_to(&mut self, source_bit: usize, check: &[usize]) {
        check.iter().for_each(|adjacent_bit| {
            self.adjacencies[source_bit].insert(*adjacent_bit);
        });
    }
}

struct AdjacentBitsGetter<'a> {
    source_bit: usize,
    adjacent_bits: BTreeMap<usize, usize>,
    adjacencies: &'a [BTreeSet<usize>],
}

impl<'a> AdjacentBitsGetter<'a> {
    fn get_adjacent_bits_with_depth(&mut self, depth: usize) -> Vec<usize> {
        if self.adjacent_bits.is_empty() {
            self.recursively_set_adjacent_bits(self.source_bit, depth);
        }
        self.get_adjacent_bits_as_vec()
    }

    fn recursively_set_adjacent_bits(&mut self, bit: usize, depth: usize) {
        self.adjacent_bits.insert(bit, depth);
        if depth > 0 {
            self.adjacencies[bit].iter().for_each(|b| {
                if self.has_not_been_set_with_depth(b, depth - 1) {
                    self.recursively_set_adjacent_bits(*b, depth - 1);
                }
            })
        }
    }

    fn has_not_been_set_with_depth(&self, bit: &usize, depth: usize) -> bool {
        self.adjacent_bits
            .get(bit)
            .map(|last_set_depth| *last_set_depth <= depth)
            .unwrap_or(true)
    }

    fn get_adjacent_bits_as_vec(&self) -> Vec<usize> {
        self.adjacent_bits.keys().cloned().collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn without_recursion_every_bit_is_adjancent_only_to_itself() {
        let mut adjacency = Adjacency::with_n_bits(5);
        adjacency.update_from_check(&[0, 1, 2, 3, 4]);
        for bit in 0..5 {
            assert_eq!(adjacency.get_bits_adjacent_to(bit), vec![bit]);
        }
    }

    #[test]
    fn without_check_every_bit_is_adjacent_only_to_itself() {
        let mut adjacency = Adjacency::with_n_bits(5);
        adjacency.set_recursion_depth(10);
        for bit in 0..5 {
            assert_eq!(adjacency.get_bits_adjacent_to(bit), vec![bit]);
        }
    }

    #[test]
    fn with_one_full_check_and_one_level_of_recursions_every_bit_has_all_adjacent_bits() {
        let mut adjacency = Adjacency::with_n_bits(5);
        adjacency.set_recursion_depth(1);
        adjacency.update_from_check(&[0, 1, 2, 3, 4]);
        for bit in 0..5 {
            assert_eq!(adjacency.get_bits_adjacent_to(bit), vec![0, 1, 2, 3, 4]);
        }
    }

    #[test]
    fn with_few_checks_and_two_level_of_recursions_every_bit_has_some_adjacent_bits() {
        let mut adjacency = Adjacency::with_n_bits(6);
        adjacency.set_recursion_depth(2);

        adjacency.update_from_check(&[0, 1]);
        adjacency.update_from_check(&[2, 3]);
        adjacency.update_from_check(&[1, 2, 4]);
        adjacency.update_from_check(&[3, 5]);

        assert_eq!(adjacency.get_bits_adjacent_to(0), vec![0, 1, 2, 4]);
        assert_eq!(adjacency.get_bits_adjacent_to(1), vec![0, 1, 2, 3, 4]);
        assert_eq!(adjacency.get_bits_adjacent_to(2), vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(adjacency.get_bits_adjacent_to(3), vec![1, 2, 3, 4, 5]);
        assert_eq!(adjacency.get_bits_adjacent_to(4), vec![0, 1, 2, 3, 4]);
        assert_eq!(adjacency.get_bits_adjacent_to(5), vec![2, 3, 5]);
    }
}

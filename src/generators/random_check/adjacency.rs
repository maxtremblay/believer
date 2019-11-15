use std::collections::{BTreeMap, BTreeSet};

pub(super) struct Adjacency {
    adjacencies: Vec<BTreeSet<usize>>,
    recursion_depth: usize,
}

impl Adjacency {

    // ***** Construction *****

    pub(super) fn new() -> Self {
        Self {
            recursion_depth: 0,
            adjacencies: Vec::new(),
        }
    }

    pub(super) fn initialize_adjacencies(&mut self, n_bits: usize) {
        self.adjacencies = (0..n_bits)
            .map(|b| {
                let mut set = BTreeSet::new();
                set.insert(b);
                set
            })
            .collect();
    }

    pub(super) fn set_recursion_depth_from_girth(&mut self, girth: usize) {
        self.recursion_depth = girth / 2;
    }

    // ***** Get adjacent bits *****

    pub(super) fn get_bits_adjacent_to(&self, bit: usize) -> Vec<usize> {
        if self.is_out_of_bound(bit) {
            Vec::new()
        } else {
            self.get_bits_adjacent_to_inbound_bit(bit)
        }
    }

    fn is_out_of_bound(&self, bit:usize) -> bool {
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
        check.iter().for_each(|bit| self.set_bits_in_check_adjacent_to(*bit, check));
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
    adjacencies: &'a [BTreeSet<usize>]
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
        self.adjacent_bits[bit] <= depth
    }

    fn get_adjacent_bits_as_vec(&self) -> Vec<usize> {
        self.adjacent_bits.keys().cloned().collect()
    }
}

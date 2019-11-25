use super::DecodingResult;

/// An interface for simulation result. 
#[derive(Clone, Copy)]
pub struct SimulationResult {
    n_successes: u64,
    n_failures: u64,
}

impl SimulationResult {
    // ***** Construction *****

    // Creates a new empty `SimulationResult`.
    pub(super) fn new() -> Self {
        Self { n_successes: 0, n_failures: 0 }
    }

    // Creates a new `SimulationResult` from the number of successes and failures.
    pub(super) fn with_n_successes_and_failures(n_successes: u64, n_failures: u64) -> Self {
        Self { n_successes, n_failures }
    }

    // Creates the worse `SimulationResult`. That is, a simulation with failure rate 1.
    pub(super) fn worse_result() -> Self {
        Self { n_successes: 0, n_failures: 1 }
    }

    // ***** Updaters *****

    pub(super) fn add_decoding_result<D: DecodingResult>(&mut self, result: D) {
        if result.is_success() {
            self.n_successes += 1;
        } else {
            self.n_failures += 1;
        }
    }

    // ***** Checkers *****

    pub(super) fn has_not_at_least_one_success_and_one_failure(&self) -> bool {
        self.n_successes == 0 || self.n_failures == 0
    }

    // ***** Getters *****

    pub(super) fn combine_with(&self, other: SimulationResult) -> Self {
        Self {
            n_successes: self.n_successes + other.n_successes,
            n_failures: self.n_failures + other.n_failures,
        }
    }

    /// Get the effective failure rate of `self` for a given code `dimension`. 
    ///
    /// This is the equivalent failure rate per bit if `dimension` similar bits without error
    /// correction where used. 
    ///
    ///
    /// # Example 
    ///
    /// ```
    /// use believer::SimulationResult;
    /// let result = SimulationResult::from_n_successes_and_failures(9, 16);
    /// assert_eq!(result.get_effective_failure_rate_for_code_dimension(2), 0.4);
    /// ``` 
    pub fn get_effective_failure_rate_for_code_dimension(&self, dimension: u32) -> f64 {
        1.0 - self.get_effective_success_rate_for_code_dimension(dimension)
    }

    /// Get the effective success rate of `self` for a given code `dimension`. 
    ///
    /// This is the equivalent success rate per bit if `dimension` similar bits without error
    /// correction where used. 
    ///
    ///
    /// # Example 
    ///
    /// ```
    /// use believer::SimulationResult;
    /// let result = SimulationResult::from_n_successes_and_failures(9, 16);
    /// assert_eq!(result.get_effective_success_rate_for_code_dimension(2), 0.6);
    /// ```
    pub fn get_effective_success_rate_for_code_dimension(&self, dimension: u32) -> f64 {
        (self.get_success_rate() as f64).powf(1.0 / dimension as f64)
    }

    /// Get the failure rate of `self`.
    /// 
    /// # Example 
    /// 
    /// ```
    /// use believer::SimulationResult;
    /// let result = SimulationResult::from_n_successes_and_failures(9, 16);
    /// assert_eq!(result.get_failure_rate(), 0.64);
    /// ```
    pub fn get_failure_rate(&self) -> f64 {
        self.n_failures as f64 / self.get_n_iterations() as f64
    }

    /// Get the success rate of `self`.
    /// 
    /// # Example 
    /// 
    /// ```
    /// use believer::SimulationResult;
    /// let result = SimulationResult::from_n_successes_and_failures(9, 16);
    /// assert_eq!(result.get_failure_rate(), 0.36);
    /// ```
    pub fn get_success_rate(&self) -> f64 {
        self.n_successes as f64 / self.get_n_iterations() as f64
    }

    pub fn get_n_iterations(&self) -> u64 {
        self.n_failures + self.n_successes
    }

    pub fn get_n_failures(&self) -> u64 {
        self.n_failures
    }
    
    pub fn get_n_successes(&self) -> u64 {
        self.n_successes
    }
}
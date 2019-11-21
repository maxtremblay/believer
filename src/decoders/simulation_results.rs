/// An interface for simulation result. 
pub struct SimulationResult {
    n_successes: u64,
    n_failures: u64,
}

impl SimulationResult {
    // ***** Construction *****

    /// Creates a new `SimulationResult` from the number of successes and failures.
    pub fn with_n_successes_and_failures(n_successes: u64, n_failures: u64) -> Self {
        Self { n_successes, n_failures }
    }

    /// Creates the worse `SimulationResult`. That is, a simulation with failure rate 1.
    pub fn worse_result() -> Self {
        Self { n_successes: 0, n_failures: 1 }
    }

    // ***** Getters *****

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
        self.n_failures as f64 / self.get_total_n_iterations() as f64
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
        self.n_successes as f64 / self.get_total_n_iterations() as f64
    }

    fn get_total_n_iterations(&self) -> u64 {
        self.n_failures + self.n_successes
    }
}
use crate::arithmetic::{observation, power_scale, rescale, Sum};
use crate::Error;

/// Allocation-free univariate mean and centered-moment accumulator.
///
/// Uses compensated Welford updates and Chan/Pébay pairwise merges. Stores the
/// count, a mean with a low residual, and a compensated centered sum `M2` in O(1)
/// space. In addition to the [crate contracts](crate), differences from the mean,
/// merge differences, and the **unnormalized** `M2` must remain finite. Updates
/// can therefore fail even if a normalized variance would fit; use batch queries
/// for their wider moment range. Underflow in individual moment updates is
/// ordinary floating-point underflow and may lose tiny contributions.
///
/// Rejected updates leave the complete state unchanged. Count addition is
/// checked; ratios round counts above `2^53` to `f64`. Merge order affects rounding,
/// so neither exact associativity nor bitwise agreement with batch is promised.
/// Definitions: [GSL running statistics](https://www.gnu.org/software/gsl/doc/html/rstat.html).
/// Pairwise formulas: [Pébay, SAND2008-6212, §2](https://www.osti.gov/servlets/purl/1028931).
///
/// ```
/// use gsl_oxide_statistics::RunningStatistics;
/// let mut left = RunningStatistics::new();
/// for x in [1., 2.] { left.push(x)?; }
/// let mut right = RunningStatistics::default();
/// for x in [3., 4., 5.] { right.push(x)?; }
/// left.merge(&right)?;
/// assert_eq!(left.count(), 5);
/// assert_eq!(left.mean()?, 3.);
/// assert_eq!(left.sample_variance()?, 2.5);
/// left.push(6.)?;
/// assert_eq!(right.count(), 3);
/// # Ok::<(), gsl_oxide_statistics::Error>(())
/// ```
#[derive(Clone, Debug, Default)]
pub struct RunningStatistics {
    count: usize,
    mean: Sum,
    m2: Sum,
}

impl RunningStatistics {
    /// Creates an empty accumulator, equivalent to [`Default::default`].
    /// Numerical queries return [`Error::EmptyInput`] until a value is added.
    /// See the [GSL initialization reference](https://www.gnu.org/software/gsl/doc/html/rstat.html#initializing-the-accumulator).
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of successfully added observations, including merged observations.
    /// Always exact as `usize`; initially zero. [GSL reference](https://www.gnu.org/software/gsl/doc/html/rstat.html#adding-data-to-the-accumulator).
    pub fn count(&self) -> usize {
        self.count
    }

    /// Restores the empty state for reuse; cannot fail.
    /// [GSL reset reference](https://www.gnu.org/software/gsl/doc/html/rstat.html#initializing-the-accumulator).
    pub fn reset(&mut self) {
        *self = Self::new();
    }

    /// Adds one finite `value`, updating the mean and centered sum of squares.
    ///
    /// Returns [`Error::NonFiniteData`] for NaN/infinity, [`Error::CountOverflow`]
    /// for count overflow, or [`Error::NumericalRange`] for the internal range
    /// limits on [`RunningStatistics`]. On error nothing changes. A first finite
    /// value establishes its own mean and positive-zero population variance.
    /// [GSL add reference](https://www.gnu.org/software/gsl/doc/html/rstat.html#adding-data-to-the-accumulator).
    pub fn push(&mut self, value: f64) -> Result<(), Error> {
        observation(value)?;
        let count = self.count.checked_add(1).ok_or(Error::CountOverflow)?;
        if self.count == 0 {
            *self = Self {
                count,
                mean: Sum {
                    value,
                    correction: 0.0,
                },
                m2: Sum::default(),
            };
            return Ok(());
        }
        let mut delta = Sum {
            value,
            correction: 0.0,
        };
        delta.add(-self.mean.value)?;
        delta.add(-self.mean.correction)?;
        let mut candidate = self.clone();
        // Retain the subtraction residual as well as the mean-update residual:
        // e.g. [1e16, 1, -1e16] must not lose the middle observation in delta.
        candidate.mean.add(delta.value / count as f64)?;
        candidate.mean.add(delta.correction / count as f64)?;
        candidate.mean.total()?;
        let scale = power_scale(delta.total()?.abs());
        let d = delta.value / scale + delta.correction / scale;
        // Welford: delta^2 * n/(n+1), nonnegative without a roundoff clamp.
        let increment = rescale(d * d * (self.count as f64 / count as f64), scale, scale)?;
        candidate.m2.add(increment)?;
        candidate.m2.total()?;
        candidate.count = count;
        *self = candidate;
        Ok(())
    }

    /// Merges all observations in `other`, which is left unchanged.
    ///
    /// Empty states are identities. Uses count-weighted means and the pairwise
    /// correction `delta^2 * n_a*n_b/(n_a+n_b)` for `M2`; supports later pushes.
    /// Returns [`Error::CountOverflow`] or [`Error::NumericalRange`] under the
    /// struct's documented limits, without changing `self` on failure.
    /// [Pébay §2](https://www.osti.gov/servlets/purl/1028931);
    /// [GSL estimator definitions](https://www.gnu.org/software/gsl/doc/html/rstat.html#current-statistics).
    pub fn merge(&mut self, other: &Self) -> Result<(), Error> {
        let count = self
            .count
            .checked_add(other.count)
            .ok_or(Error::CountOverflow)?;
        if other.count == 0 {
            return Ok(());
        }
        if self.count == 0 {
            *self = other.clone();
            return Ok(());
        }
        // Anchor the mean in the larger block to keep the interpolation small.
        let (large, small) = if self.count >= other.count {
            (&*self, other)
        } else {
            (other, &*self)
        };
        let mut delta = small.mean;
        delta.add(-large.mean.value)?;
        delta.add(-large.mean.correction)?;
        let mut mean = large.mean;
        let ratio = small.count as f64 / count as f64;
        mean.add(delta.value * ratio)?;
        mean.add(delta.correction * ratio)?;
        mean.total()?;
        let scale = power_scale(delta.total()?.abs());
        let d = delta.value / scale + delta.correction / scale;
        let coefficient = (large.count as f64 / count as f64) * small.count as f64;
        let cross = rescale(d * d * coefficient, scale, scale)?;
        let mut m2 = large.m2;
        m2.add(small.m2.value)?;
        m2.add(small.m2.correction)?;
        m2.add(cross)?;
        m2.total()?;
        *self = Self { count, mean, m2 };
        Ok(())
    }

    /// Current finite arithmetic mean, requiring one observation.
    /// Returns [`Error::EmptyInput`] if empty. [GSL mean](https://www.gnu.org/software/gsl/doc/html/rstat.html#current-statistics).
    pub fn mean(&self) -> Result<f64, Error> {
        if self.count == 0 {
            Err(Error::EmptyInput)
        } else {
            self.mean.total()
        }
    }

    /// Current population variance `M2 / n` about the computed mean.
    /// Requires one observation; a singleton/constant returns positive zero.
    /// Returns [`Error::EmptyInput`] if empty; see the struct's moment limits.
    /// [GSL moment definitions](https://www.gnu.org/software/gsl/doc/html/rstat.html#current-statistics).
    pub fn population_variance(&self) -> Result<f64, Error> {
        self.mean()?;
        Ok(self.m2.total()? / self.count as f64)
    }

    /// Current sample variance `M2 / (n - 1)`, matching GSL's running variance.
    /// Requires two observations; returns [`Error::EmptyInput`] if empty or
    /// [`Error::InsufficientData`] for a singleton. Constants return positive zero.
    /// See the struct's moment limits and [GSL variance](https://www.gnu.org/software/gsl/doc/html/rstat.html#current-statistics).
    pub fn sample_variance(&self) -> Result<f64, Error> {
        self.mean()?;
        if self.count < 2 {
            return Err(Error::InsufficientData);
        }
        Ok(self.m2.total()? / (self.count - 1) as f64)
    }

    /// Finite square root of [`Self::population_variance`], requiring one value.
    /// Has the same empty-input error and range limits; singleton/constant gives
    /// positive zero. [GSL deviation definition](https://www.gnu.org/software/gsl/doc/html/rstat.html#current-statistics).
    pub fn population_standard_deviation(&self) -> Result<f64, Error> {
        Ok(self.population_variance()?.sqrt())
    }

    /// Finite square root of [`Self::sample_variance`], requiring two values.
    /// Has the same empty/insufficient-data errors and range limits; constants
    /// give positive zero. [GSL deviation](https://www.gnu.org/software/gsl/doc/html/rstat.html#current-statistics).
    pub fn sample_standard_deviation(&self) -> Result<f64, Error> {
        Ok(self.sample_variance()?.sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Returns the count and internal floating-point bit patterns of accumulator `s` for atomicity checks.
    fn snapshot(s: &RunningStatistics) -> (usize, [u64; 4]) {
        (
            s.count,
            [
                s.mean.value.to_bits(),
                s.mean.correction.to_bits(),
                s.m2.value.to_bits(),
                s.m2.correction.to_bits(),
            ],
        )
    }

    /// Checks count overflow is atomic.
    #[test]
    fn count_overflow_is_atomic() {
        let mut state = RunningStatistics {
            count: usize::MAX - 1,
            mean: Sum::default(),
            m2: Sum::default(),
        };
        // A rounded f64 count ratio must not affect exact usize counting.
        let mut merged = state.clone();
        let mut singleton = RunningStatistics::new();
        singleton.push(0.0).unwrap();
        merged.merge(&singleton).unwrap();
        assert_eq!(merged.count(), usize::MAX);
        state.push(0.0).unwrap();
        assert_eq!(state.count(), usize::MAX);
        assert_eq!(state.mean(), Ok(0.0));
        assert_eq!(state.sample_variance(), Ok(0.0));
        let before = snapshot(&state);
        assert_eq!(state.push(0.0), Err(Error::CountOverflow));
        assert_eq!(state.merge(&singleton), Err(Error::CountOverflow));
        assert_eq!(snapshot(&state), before);
        state.merge(&RunningStatistics::new()).unwrap();
        assert_eq!(snapshot(&state), before);
    }

    /// Checks internal moment overflow is atomic.
    #[test]
    fn internal_moment_overflow_is_atomic() {
        let mut state = RunningStatistics {
            count: 2,
            mean: Sum::default(),
            m2: Sum {
                value: f64::MAX,
                correction: 0.0,
            },
        };
        let before = snapshot(&state);
        assert_eq!(state.push(1e154), Err(Error::NumericalRange));
        assert_eq!(snapshot(&state), before);
        let other = state.clone();
        assert_eq!(state.merge(&other), Err(Error::NumericalRange));
        assert_eq!(snapshot(&state), before);
    }
}

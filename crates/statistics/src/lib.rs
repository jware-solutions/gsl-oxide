//! Independent, allocation-free statistics on borrowed `f64` slices and streams.
//!
//! Variances use the computed sample mean: population variance divides the
//! centered sum of squares by `n`, sample variance by `n - 1`. These are distinct
//! from variance about a known, fixed population mean. Standard deviations are
//! square roots of the corresponding variances. Covariance uses `n - 1`;
//! correlation is the centered, normalized Pearson coefficient.
//!
//! Every successful result is finite. Empty samples, insufficient degrees of
//! freedom, and participating NaN/infinite observations are errors. Constants
//! have positive-zero variance/deviation; zero-variance correlation is undefined.
//! Ordinary IEEE 754 underflow is allowed (including a result rounding to zero);
//! there is no near-constant cutoff. Scaling and compensated centered sums avoid
//! raw-moment subtraction, but do not promise correct rounding or bitwise
//! agreement across platforms or accumulation orders. Extremely disparate sample
//! magnitudes can lose contributions during floating-point scaling.
//!
//! Batch moments are scaled before division, so an unrepresentable *unnormalized*
//! sum of squares need not cause failure. Variance must fit in `f64`, including
//! for a standard-deviation query: that query can fail even when the mathematical
//! deviation is finite. Correlation operates on normalized centered data and does
//! not require representable variances. Streaming has the additional internal
//! range limits documented on [`RunningStatistics`]. Counts use checked `usize`
//! arithmetic; ratios use `f64`, so counts above `2^53` can incur integer-conversion
//! roundoff. All operations use O(1) auxiliary storage and never mutate samples.
//!
//! Weights are finite, nonnegative real reliability weights. Zero (including
//! `-0.0`) removes an observation **before** validating its value. Positive
//! weights are divided by their maximum; each resulting ratio must be a normal
//! positive `f64` (at least [`f64::MIN_POSITIVE`]), otherwise [`Error::WeightRange`]
//! is returned. This explicit range limit avoids silently discarding small
//! weights. Raw weight sums, squared weights, and weight/value products need not
//! be representable. [`weighted_reliability_variance`] uses `W - W2/W`, never the
//! frequency-weight denominator `W - 1`.
//!
//! Definitions: [GSL statistics](https://www.gnu.org/software/gsl/doc/html/statistics.html)
//! and [running statistics](https://www.gnu.org/software/gsl/doc/html/rstat.html).
//!
//! ```
//! use gsl_oxide_statistics::{mean, population_variance, sample_variance};
//! let data = [1., 2., 3., 4., 5.];
//! assert_eq!(mean(&data)?, 3.);
//! assert_eq!(population_variance(&data)?, 2.);
//! assert_eq!(sample_variance(&data)?, 2.5);
//! # Ok::<(), gsl_oxide_statistics::Error>(())
//! ```
//!
//! ```
//! use gsl_oxide_statistics::{weighted_mean, weighted_reliability_variance};
//! assert_eq!(weighted_mean(&[1., 4.], &[1., 3.])?, 3.25);
//! assert!((weighted_reliability_variance(&[1., 4.], &[1., 3.])? - 4.5).abs() < 1e-14);
//! # Ok::<(), gsl_oxide_statistics::Error>(())
//! ```

#![forbid(unsafe_code)]

mod arithmetic;
mod running;
mod slice;
mod weighted;

pub use running::RunningStatistics;
pub use slice::{
    correlation, mean, population_standard_deviation, population_variance, sample_covariance,
    sample_standard_deviation, sample_variance,
};
pub use weighted::{weighted_mean, weighted_population_variance, weighted_reliability_variance};

/// Failure of a checked statistics operation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    /// No observations were supplied to an unweighted numerical query.
    EmptyInput,
    /// Fewer than two observations (or positive weights) for a sample estimator.
    InsufficientData,
    /// Paired data or data and weights have different lengths.
    LengthMismatch,
    /// A participating observation is NaN or infinite.
    NonFiniteData,
    /// A weight is negative, NaN, or infinite.
    InvalidWeight,
    /// No positive weights, including an empty weighted sample.
    ZeroTotalWeight,
    /// A positive weight divided by the maximum weight is subnormal or zero.
    WeightRange,
    /// Correlation is undefined because at least one input is constant.
    ZeroVariance,
    /// A required result/state overflows or arithmetic leaves its valid range.
    NumericalRange,
    /// Adding observation counts would overflow `usize`.
    CountOverflow,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::EmptyInput => "empty sample",
            Self::InsufficientData => "insufficient observations for the estimator",
            Self::LengthMismatch => "sample lengths differ",
            Self::NonFiniteData => "participating observation is not finite",
            Self::InvalidWeight => "weight must be finite and nonnegative",
            Self::ZeroTotalWeight => "no positive weights",
            Self::WeightRange => "normalized positive weight is below the supported normal range",
            Self::ZeroVariance => "correlation requires two nonconstant inputs",
            Self::NumericalRange => "required statistics arithmetic is outside the supported range",
            Self::CountOverflow => "observation count overflows usize",
        })
    }
}

impl std::error::Error for Error {}

//! Independent safe Rust quadrature for finite intervals, breakpoints, infinite
//! intervals, and integrable endpoint singularities.
//!
//! All functions return `(integral, estimated_absolute_error)`. Estimates are
//! numerical diagnostics, **not rigorous bounds for arbitrary callbacks**. Finite
//! sampling can miss narrow features; supply useful breakpoints when known.
//!
//! Tolerances must be finite and nonnegative. At least one must be positive;
//! pure relative requests require `rel_tol >= 64 * f64::EPSILON`. Success means
//! `error <= max(abs_tol, rel_tol * abs(integral))`, with finite results and a
//! nonnegative estimate. Cancellation generally requires an absolute tolerance.
//! Unattainable accuracy returns an error, including requests below the measured
//! roundoff floor. Intermediate arithmetic must remain representable as `f64`.
//!
//! `max_evals` is a hard cap on calls to the original `FnMut` closure, including
//! initialization, all segments, both infinite tails, and singular-method mesh
//! refinements. A zero budget is valid but cannot integrate a nonempty interval.
//! Invalid arguments are rejected before any call. Valid zero-length finite
//! intervals return `(0.0, 0.0)` without calls. Reversal negates only the integral.
//! Any evaluated NaN or infinity is an error. Callbacks should describe a fixed
//! mathematical function even though they may keep counters or other state.
//!
//! Finite integration uses adaptive Gauss 7 / Kronrod 15 subdivision. Infinite
//! integration uses rational substitutions; endpoint integration uses tanh–sinh.
//! See the [GSL behavior reference](https://www.gnu.org/software/gsl/doc/html/integration.html)
//! and [DLMF quadrature mathematics](https://dlmf.nist.gov/3.5).
//!
//! ```
//! use gsl_oxide_integration::integrate;
//! let (area, error) = integrate(f64::sin, 0.0, std::f64::consts::PI,
//!                               1e-12, 1e-12, 10_000)?;
//! assert!((area - 2.0).abs() < 1e-12);
//! assert!(error <= 2e-12);
//! # Ok::<(), gsl_oxide_integration::Error>(())
//! ```

#![forbid(unsafe_code)]

mod adaptive;
mod infinite;
mod rule;
mod singular;

pub use infinite::{integrate_infinite, integrate_lower_infinite, integrate_upper_infinite};
pub use singular::integrate_singular;

/// Failure of a checked integration operation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    /// Invalid endpoints, breakpoints, or tolerances.
    InvalidInput,
    /// An evaluated original callback returned NaN or infinity.
    NonFiniteIntegrand,
    /// Another original callback call would exceed `max_evals`.
    EvaluationLimit,
    /// Roundoff, stagnation, or endpoint resolution prevents requested accuracy.
    Roundoff,
    /// A required interval or its quadrature nodes cannot be represented distinctly.
    Subdivision,
    /// A transformation, sum, or weighted result overflowed or broke down.
    NumericalBreakdown,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => "invalid integration arguments",
            Self::NonFiniteIntegrand => "integrand returned a non-finite value",
            Self::EvaluationLimit => "integrand evaluation budget exhausted",
            Self::Roundoff => "roundoff or stagnation prevents convergence",
            Self::Subdivision => "required subdivision or nodes cannot be represented",
            Self::NumericalBreakdown => "non-finite internal integration arithmetic",
        })
    }
}

impl std::error::Error for Error {}

/// Integrates a function over finite bounds using adaptive Gauss–Kronrod rules.
///
/// Returns `(integral, estimated_absolute_error)` under the [crate contracts](crate).
/// Smooth functions and ordinary local difficulties are supported; use
/// [`integrate_with_breakpoints`] for known interior discontinuities or peaks,
/// and [`integrate_singular`] for the documented endpoint singularities.
pub fn integrate<F>(
    f: F,
    a: f64,
    b: f64,
    abs_tol: f64,
    rel_tol: f64,
    max_evals: usize,
) -> Result<(f64, f64), Error>
where
    F: FnMut(f64) -> f64,
{
    integrate_with_breakpoints(f, a, b, &[], abs_tol, rel_tol, max_evals)
}

/// Integrates finite bounds with caller-supplied interior breakpoints.
///
/// Breakpoints must be finite, strictly increasing, and strictly inside
/// `(min(a,b), max(a,b))`, regardless of integration direction. Endpoints,
/// duplicates, and out-of-range points are invalid. A zero-length interval
/// accepts only an empty slice. The complete slice is validated before calls.
/// All segments share **one** callback budget and total error criterion, including
/// cancellation between segments. An empty slice is equivalent to [`integrate`].
/// Returns `(integral, estimated_absolute_error)`; see the [crate contracts](crate).
///
/// ```
/// use gsl_oxide_integration::integrate_with_breakpoints;
/// let (area, _) = integrate_with_breakpoints(
///     |x| if x < 0.3 { 1.0 } else { 2.0 }, 0.0, 1.0, &[0.3],
///     1e-12, 0.0, 1000)?;
/// assert!((area - 1.7).abs() < 1e-12);
/// # Ok::<(), gsl_oxide_integration::Error>(())
/// ```
pub fn integrate_with_breakpoints<F>(
    f: F,
    a: f64,
    b: f64,
    breakpoints: &[f64],
    abs_tol: f64,
    rel_tol: f64,
    max_evals: usize,
) -> Result<(f64, f64), Error>
where
    F: FnMut(f64) -> f64,
{
    validate_finite(a, b, abs_tol, rel_tol)?;
    let (lo, hi, sign) = ordered(a, b);
    let mut previous = lo;
    for &point in breakpoints {
        if !point.is_finite() || point <= previous || point >= hi {
            return Err(Error::InvalidInput);
        }
        previous = point;
    }
    if a == b {
        return Ok((0.0, 0.0));
    }
    let mut budget = Budget::new(f, max_evals);
    let segments = breakpoints
        .iter()
        .copied()
        .chain([hi])
        .scan(lo, |left, right| {
            let segment = (*left, right, 0);
            *left = right;
            Some(segment)
        });
    let (value, error) = adaptive::run(&mut |x, _| budget.eval(x), segments, abs_tol, rel_tol)?;
    Ok((sign * value, error))
}

fn validate_tolerances(abs: f64, rel: f64) -> Result<(), Error> {
    if !abs.is_finite()
        || !rel.is_finite()
        || abs < 0.0
        || rel < 0.0
        || (abs == 0.0 && rel < 64.0 * f64::EPSILON)
    {
        Err(Error::InvalidInput)
    } else {
        Ok(())
    }
}

fn validate_finite(a: f64, b: f64, abs: f64, rel: f64) -> Result<(), Error> {
    validate_tolerances(abs, rel)?;
    if !a.is_finite() || !b.is_finite() {
        Err(Error::InvalidInput)
    } else {
        Ok(())
    }
}

fn ordered(a: f64, b: f64) -> (f64, f64, f64) {
    if a <= b {
        (a, b, 1.0)
    } else {
        (b, a, -1.0)
    }
}

fn finite(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::NumericalBreakdown)
    }
}

fn target(abs: f64, rel: f64, value: f64) -> f64 {
    // Saturation keeps the target finite even for huge valid tolerances.
    abs.max((rel * value.abs()).min(f64::MAX))
}

struct Budget<F> {
    f: F,
    used: usize,
    limit: usize,
}

impl<F: FnMut(f64) -> f64> Budget<F> {
    fn new(f: F, limit: usize) -> Self {
        Self { f, used: 0, limit }
    }

    fn eval(&mut self, x: f64) -> Result<f64, Error> {
        finite(x)?;
        if self.used >= self.limit {
            return Err(Error::EvaluationLimit);
        }
        self.used += 1; // used < limit <= usize::MAX, so this cannot overflow.
        let y = (self.f)(x);
        if y.is_finite() {
            Ok(y)
        } else {
            Err(Error::NonFiniteIntegrand)
        }
    }
}

// Neumaier compensation retains small residuals in cancellation-heavy sums.
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}

impl Sum {
    fn add(&mut self, x: f64) -> Result<(), Error> {
        let next = finite(self.value + x)?;
        let residual = if self.value.abs() >= x.abs() {
            (self.value - next) + x
        } else {
            (x - next) + self.value
        };
        self.correction = finite(self.correction + residual)?;
        self.value = next;
        Ok(())
    }

    fn total(&self) -> Result<f64, Error> {
        finite(self.value + self.correction)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budget_counter_cannot_overflow() {
        let mut calls = 0;
        let mut budget = Budget {
            f: |_| {
                calls += 1;
                1.0
            },
            used: usize::MAX - 1,
            limit: usize::MAX,
        };
        assert_eq!(budget.eval(0.0), Ok(1.0));
        assert_eq!(budget.used, usize::MAX);
        assert_eq!(budget.eval(0.0), Err(Error::EvaluationLimit));
        assert_eq!(calls, 1);
    }

    #[test]
    fn compensated_residual_and_nonfinite_accumulation() {
        let mut sum = Sum::default();
        for x in [1e16, 1.0, -1e16] {
            sum.add(x).unwrap();
        }
        assert_eq!(sum.total(), Ok(1.0));
        let mut sum = Sum::default();
        sum.add(f64::MAX).unwrap();
        assert_eq!(sum.add(f64::MAX), Err(Error::NumericalBreakdown));
    }
}

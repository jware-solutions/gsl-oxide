//! Random distributions and sampling utilities corresponding to GSL 2.8 `randist`.
//!
//! Sampling takes a caller-owned [`RngCore`] (rand_core 0.6 / rand 0.8); this crate
//! does not seed or implement a generator. Names ending in `_density` denote
//! continuous densities, `_pdf` discrete probabilities, and `_log_density` or
//! `_log_pdf` their natural logarithms. Scalar distributions are re-exported
//! here, alongside multivariate distributions and generic selection utilities.
//!
//! Parameters must be finite, with positive scales/shapes/degrees of freedom.
//! Probabilities are in `[0,1]` unless an endpoint makes the distribution
//! undefined. Scalar density arguments may be infinite (density zero), but not NaN.
//! Densities may be infinite at a support boundary; logarithms may be infinite.
//! Sampling overflow and unrepresentable intermediate arithmetic return
//! [`Error::NumericalRange`]; ordinary underflow may round to zero. Discrete
//! samples use `u32` and report overflow rather than saturating. Invalid slice
//! dimensions return an error. See individual functions for parameter order,
//! support, and matrix layout. Results are not bitwise GSL RNG sequences.
//!
//! Algorithms and tables adapted from GNU GSL 2.8, GPL-3.0-or-later:
//! James Theiler, Brian Gough, Charles Karney, Jochen Voss, David Morrison,
//! Gavin E. Crooks, and other GSL contributors. The Landau approximation
//! originates in CERN DENLAN/RANLAN by K. S. Kölbig.
//!
//! ```
//! use gsl_oxide_randist::{normal_density, binomial_pdf, DiscreteDistribution};
//! assert!((normal_density(0.0, 1.0)? - 0.3989422804014327).abs() < 1e-15);
//! assert!((binomial_pdf(0, 0.5, 3)? - 0.125).abs() < 1e-15);
//! let categories = DiscreteDistribution::new(&[1.0, 3.0])?;
//! assert_eq!(categories.pdf(1), 0.75);
//! # Ok::<(), gsl_oxide_randist::Error>(())
//! ```

#![forbid(unsafe_code)]
#![warn(clippy::missing_docs_in_private_items)]

/// Scalar continuous distributions and their density functions.
mod continuous;
/// Integer-valued distributions, probability functions, and alias tables.
mod discrete;
/// CERN approximations for standard Landau sampling and density evaluation.
mod landau;
/// Dirichlet, multivariate normal, and Wishart distributions and estimators.
mod multivariate;
/// Normal-distribution samplers, tail samplers, and densities.
mod normal;
/// Generic selection utilities and uniform direction samplers.
mod sampling;

pub use continuous::*;
pub use discrete::*;
pub use landau::{landau, landau_density};
pub use multivariate::*;
pub use normal::*;
pub use rand_core::RngCore;
pub use sampling::*;

/// Failure of a checked distribution or sampling operation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    /// A distribution parameter or observation is outside its domain.
    Domain,
    /// An input slice is empty where a nonempty input is required.
    EmptyInput,
    /// Input dimensions disagree or a matrix is not square.
    DimensionMismatch,
    /// A requested sample is larger than the source without replacement.
    SampleTooLarge,
    /// A required calculation or sampled value is unrepresentable.
    NumericalRange,
}

impl std::fmt::Display for Error {
    /// Writes this error's human-readable description to the formatter `f`.
    ///
    /// Returns the formatter's result without changing the error.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Domain => "parameter or observation outside distribution domain",
            Self::EmptyInput => "empty input",
            Self::DimensionMismatch => "incompatible dimensions",
            Self::SampleTooLarge => "sample larger than source",
            Self::NumericalRange => "distribution arithmetic outside representable range",
        })
    }
}
impl std::error::Error for Error {}

/// Natural logarithm of `sqrt(2*pi)`, used in density normalization.
pub(crate) const LOG_SQRT_2_PI: f64 = 0.9189385332046727;

/// Checks that parameter `x` is finite and strictly positive.
///
/// Returns `Ok(())` on success, or `Error::Domain` for zero, negative values,
/// infinity, or NaN.
pub(crate) fn positive(x: f64) -> Result<(), Error> {
    if x > 0.0 && x.is_finite() {
        Ok(())
    } else {
        Err(Error::Domain)
    }
}
/// Checks that probability `p` lies in the closed interval `[0, 1]`.
///
/// Returns `Ok(())` for valid probabilities, including both endpoints, and
/// `Error::Domain` for all other values, including NaN.
pub(crate) fn probability(p: f64) -> Result<(), Error> {
    if (0.0..=1.0).contains(&p) {
        Ok(())
    } else {
        Err(Error::Domain)
    }
}
/// Checks that scalar observation `x` is not NaN.
///
/// Returns `Ok(())` for finite values and either infinity. Individual density
/// functions determine whether `x` is inside their support.
pub(crate) fn observation(x: f64) -> Result<(), Error> {
    if x.is_nan() {
        Err(Error::Domain)
    } else {
        Ok(())
    }
}
/// Returns computed value `x` unchanged when it is finite.
///
/// NaN and infinities return `Error::NumericalRange`; underflowed zeros remain
/// valid. Sampling and matrix calculations use this to check their results.
pub(crate) fn finite(x: f64) -> Result<f64, Error> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(Error::NumericalRange)
    }
}
/// Returns a computed probability, density, or logarithm `x` unless it is NaN.
///
/// NaN returns `Error::NumericalRange`. Infinities are retained for boundary
/// densities and logarithms; callers establish the required sign and support.
pub(crate) fn density(x: f64) -> Result<f64, Error> {
    if x.is_nan() {
        Err(Error::NumericalRange)
    } else {
        Ok(x)
    }
}
/// Converts a finite value `x` in `[0, u32::MAX]` into a count.
///
/// Fractional parts are truncated toward zero. Values outside this range,
/// including NaN and infinity, return `Error::NumericalRange`.
pub(crate) fn count(x: f64) -> Result<u32, Error> {
    if x >= 0.0 && x <= f64::from(u32::MAX) && x.is_finite() {
        Ok(x as u32)
    } else {
        Err(Error::NumericalRange)
    }
}
/// Uses random bits from `r` to return a uniform value on a 53-bit grid.
///
/// The result lies in `[0, 1)`: zero is possible, while one is excluded.
pub(crate) fn uniform<R: RngCore + ?Sized>(r: &mut R) -> f64 {
    (r.next_u64() >> 11) as f64 * (1.0 / 9007199254740992.0)
}
/// Uses random bits from `r` to return a uniform midpoint of a 52-bit grid.
///
/// The result lies strictly in `(0, 1)`, so logarithms and inverse transforms
/// can use it without encountering either endpoint.
pub(crate) fn open<R: RngCore + ?Sized>(r: &mut R) -> f64 {
    // Midpoints of 2^52 equal intervals: strictly between zero and one.
    ((r.next_u64() >> 12) as f64 + 0.5) * (1.0 / 4503599627370496.0)
}
/// Evaluates `ln(Gamma(x))` for a caller-validated positive finite `x`.
///
/// Returns the logarithmic component from the special-function implementation,
/// mapping its failures to `Error::NumericalRange`.
pub(crate) fn log_gamma(x: f64) -> Result<f64, Error> {
    gsl_oxide_special::gamma::ln_abs_gamma_with_sign(x)
        .map(|p| p.0)
        .map_err(|_| Error::NumericalRange)
}
/// Evaluates `ln(Beta(a, b))` for caller-validated positive finite shapes.
///
/// Returns the logarithm without exponentiating; special-function failures
/// are mapped to `Error::NumericalRange`.
pub(crate) fn log_beta(a: f64, b: f64) -> Result<f64, Error> {
    gsl_oxide_special::beta::ln_beta(a, b).map_err(|_| Error::NumericalRange)
}
/// Returns `x * ln(y)`, using the continuous convention `0 * ln(0) = 0`.
///
/// Callers supply nonnegative `y`. A zero coefficient `x` always returns zero;
/// other coefficients retain IEEE logarithmic boundary behavior.
pub(crate) fn xlog(x: f64, y: f64) -> f64 {
    if x == 0.0 {
        0.0
    } else {
        x * y.ln()
    }
}
/// Approximates the additive correction in Stirling's log-gamma formula.
///
/// For positive `x`, returns `ln(Gamma(x)) - ((x-1/2)ln(x)-x+ln(sqrt(2*pi)))`
/// through the `x^-9` term. Callers use `x >= 16` for the approximation's accuracy.
pub(crate) fn stirling(x: f64) -> f64 {
    let t = 1.0 / x;
    let t2 = t * t;
    t * (1.0 / 12.0
        + t2 * (-1.0 / 360.0 + t2 * (1.0 / 1260.0 + t2 * (-1.0 / 1680.0 + t2 / 1188.0))))
}
/// Returns `x - a - a * ln(x/a)` for positive `a` and `x`.
///
/// Uses a Taylor series when `x` is close to `a` to preserve small differences
/// in gamma and Poisson log probabilities. The caller validates both inputs.
pub(crate) fn deviance(a: f64, x: f64) -> f64 {
    let d = (x - a) / a;
    if d.abs() < 0.125 {
        let mut power = d * d;
        let mut s = 0.0;
        for k in 2..=32 {
            s += power / f64::from(k);
            power *= -d;
        }
        a * s
    } else {
        x - a - a * (x / a).ln()
    }
}

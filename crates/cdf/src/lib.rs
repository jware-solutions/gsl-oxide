//! Cumulative distributions and quantiles corresponding to all 86 functions in
//! GNU GSL 2.8's `gsl_cdf.h`.
//!
//! A `_cdf` returns `Pr(X <= x)`; a `_survival` returns `Pr(X > x)` directly,
//! retaining small upper tails. A `_quantile` inverts the CDF, and an
//! `_inverse_survival` inverts the survival function without forming `1 - p`.
//! Discrete counts use `u32`; parameter order follows GSL. Normal distributions
//! have mean zero; standard normal distributions also have unit variance.
//!
//! All functions return checked results. Parameters must be finite, with
//! strictly positive scales, shapes, and degrees of freedom. Observations can
//! be infinite but cannot be NaN. Probabilities are in `[0, 1]`; inverse
//! endpoints return the support endpoints, including infinities. Underflow
//! rounds to zero; unrepresentable finite quantiles return [`Error::Overflow`].
//! Beta and gamma shapes are limited to `1e8` (also applying to
//! distributions expressed through these functions). These limits return
//! [`Error::UnsupportedRange`]. Numerical iterations are bounded and return
//! [`Error::NonConvergence`] on failure. Accuracy is ordinary `f64` accuracy,
//! not guaranteed correct rounding; large balanced beta shapes can lose about
//! `1e-7` in absolute probability. Hypergeometric tails sum at most one million
//! recurrence terms, returning `NonConvergence` if that bound is exhausted.
//!
//! Algorithms and reference tests adapted from GNU GSL 2.8, GPL-3.0-or-later:
//! Jason H. Stover, Brian Gough, Przemyslaw Sliwa, Giulio Bottazzi, and other
//! GSL contributors. Normal quantiles use Wichura's AS 241 approximation.
//!
//! ```
//! use gsl_oxide_cdf::{normal_cdf, normal_quantile, normal_survival};
//! assert_eq!(normal_cdf(0.0, 2.0)?, 0.5);
//! assert_eq!(normal_quantile(0.5, 2.0)?, 0.0);
//! assert!(normal_survival(20.0, 1.0)? > 0.0);
//! # Ok::<(), gsl_oxide_cdf::Error>(())
//! ```

#![forbid(unsafe_code)]

mod continuous;
mod discrete;
mod inverse;
mod normal;
mod shaped;

pub use continuous::*;
pub use discrete::*;
pub use gsl_oxide_special::Error;
pub use normal::{
    normal_cdf, normal_inverse_survival, normal_quantile, normal_survival, standard_normal_cdf,
    standard_normal_inverse_survival, standard_normal_quantile, standard_normal_survival,
};
pub use shaped::*;

/// Checks that parameter `x` is finite and strictly positive, returning a domain error otherwise.
fn positive(x: f64) -> Result<(), Error> {
    if x.is_finite() && x > 0.0 {
        Ok(())
    } else {
        Err(Error::Domain)
    }
}

/// Checks that probability `p` lies in `[0,1]`, rejecting NaN and out-of-range values.
fn probability(p: f64) -> Result<(), Error> {
    if (0.0..=1.0).contains(&p) {
        Ok(())
    } else {
        Err(Error::Domain)
    }
}

/// Checks observation `x` for NaN; finite values and infinities are accepted.
fn observation(x: f64) -> Result<(), Error> {
    if x.is_nan() {
        Err(Error::Domain)
    } else {
        Ok(())
    }
}

/// Returns quantile `x`, allowing infinite endpoints for `p` equal to zero or one and reporting interior overflow.
fn finite_quantile(x: f64, p: f64) -> Result<f64, Error> {
    if x.is_nan() {
        Err(Error::NonConvergence)
    } else if x.is_infinite() && p > 0.0 && p < 1.0 {
        Err(Error::Overflow)
    } else {
        Ok(x)
    }
}

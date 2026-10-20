//! Independent, safe Rust implementations of selected GSL special functions.

#![forbid(unsafe_code)]

pub mod bessel;
pub mod beta;
pub mod error_functions;
pub mod gamma;
pub mod incomplete_beta;
pub mod incomplete_gamma;

/// Failure of a checked special-function evaluation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Error {
    /// An argument is NaN, outside the real domain, or a pole.
    Domain,
    /// A documented finite argument or order limit was exceeded.
    UnsupportedRange,
    /// The result cannot be represented as a finite `f64`.
    Overflow,
    /// An iteration limit or a numerical breakdown prevented convergence.
    NonConvergence,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Domain => "argument outside the function's real domain",
            Self::UnsupportedRange => "argument exceeds the supported numerical range",
            Self::Overflow => "result overflows f64",
            Self::NonConvergence => "numerical iteration did not converge",
        })
    }
}

impl std::error::Error for Error {}

pub(crate) const LN_SQRT_2_PI: f64 = 0.918_938_533_204_672_8;
pub(crate) const EULER: f64 = 0.577_215_664_901_532_9;

pub(crate) fn finite_exp(log: f64) -> Result<f64, Error> {
    let value = log.exp();
    if value.is_infinite() {
        Err(Error::Overflow)
    } else if value.is_nan() {
        Err(Error::NonConvergence)
    } else {
        Ok(value)
    }
}

// log(1+u)-u, retaining the quadratic term near zero.
pub(crate) fn log1pmx(u: f64) -> f64 {
    if u.abs() >= 0.125 {
        return u.ln_1p() - u;
    }
    let mut term = -u * u / 2.0;
    let mut sum = term;
    for k in 3..32 {
        term *= -u * f64::from(k - 1) / f64::from(k);
        sum += term;
    }
    sum
}

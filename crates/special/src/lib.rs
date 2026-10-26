//! Safe Rust implementations of selected GSL special functions.
//!
//! Error functions, gamma/beta, regularized incomplete gamma/beta, and
//! integer-order cylindrical Bessel functions, real exponential integrals,
//! Legendre polynomials and generalized Laguerre polynomials form the subset.
//! Each module documents its supported domains and numerical limits.
//!
//! ```
//! use gsl_oxide_special::error_functions::erf;
//! assert!((erf(1.0) - 0.8427007929497149).abs() < 2.0e-15);
//! ```
//!
//! ```
//! use gsl_oxide_special::{bessel, gamma, incomplete_gamma};
//! let (log_abs, sign) = gamma::ln_abs_gamma_with_sign(-0.5)?;
//! assert!(log_abs > 1.0 && sign == -1);
//! assert!(incomplete_gamma::upper_regularized(1.0, 100.0)? > 0.0);
//! assert!(bessel::i0_scaled(1000.0)? > 0.01);
//! # Ok::<(), gsl_oxide_special::Error>(())
//! ```
//!
//! ```
//! use gsl_oxide_special::{exponential_integrals as expint, laguerre, legendre};
//! assert!((expint::e1(1.0)? - 0.21938393439552027).abs() < 2e-15);
//! assert_eq!(expint::e2(0.0)?, 1.0);
//! assert!((expint::ei(-1.0)? + expint::e1(1.0)?).abs() < 2e-15);
//! // Scaled integrals remain useful when the unscaled results underflow/overflow.
//! assert!(expint::e1_scaled(1000.0)? > 0.0009);
//! assert!(expint::e2_scaled(1000.0)? > 0.0009);
//! assert!(expint::ei_scaled(1000.0)? > 0.001);
//! assert_eq!(legendre::polynomial(3, 0.5)?, -0.4375);
//! assert_eq!(laguerre::generalized_polynomial(2, 0.5, 1.0)?, -0.125);
//! assert_eq!(laguerre::generalized_polynomial(2, 0.0, 0.0)?, 1.0);
//! # Ok::<(), gsl_oxide_special::Error>(())
//! ```

#![forbid(unsafe_code)]

pub mod bessel;
pub mod beta;
pub mod error_functions;
pub mod exponential_integrals;
pub mod gamma;
pub mod incomplete_beta;
pub mod incomplete_gamma;
pub mod laguerre;
pub mod legendre;

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

//! Real exponential integrals, following the
//! [GSL convention](https://www.gnu.org/software/gsl/doc/html/specfunc.html#exponential-integrals).
//!
//! For `x>0`, `E_n(x)=integral_1^infinity exp(-xt)/t^n dt` for n=1,2.
//! On the negative axis we return the real continuation: `E1(x)=-Ei(-x)`
//! and `E2(x)=exp(-x)-x*E1(x)`. The positive-argument integral above does
//! not converge there. `Ei(x)=PV integral_-infinity^x exp(t)/t dt`, with
//! a principal value needed when x is positive.
//!
//! All finite real arguments are supported, except the E1/Ei poles at
//! either signed zero (`Domain`). E2 and its scaled form return 1 at
//! either zero. NaN returns `Domain`; either infinity returns
//! `UnsupportedRange`. Unrepresentable results return `Overflow`;
//! underflow returns zero with the result's sign. Series and continued
//! fractions are limited to 512 terms (`NonConvergence` on exhaustion or
//! breakdown). Roundoff is not a correct-rounding guarantee; relative
//! accuracy is inappropriate near zeros, especially the negative E2 zero.

use crate::{finite_exp, incomplete_gamma::guard, Error, EULER};

const LIMIT: usize = 512;
const ASYMPTOTIC: f64 = 50.0;

fn validate(x: f64, pole: bool) -> Result<(), Error> {
    if x.is_nan() || (pole && x == 0.0) {
        Err(Error::Domain)
    } else if !x.is_finite() {
        Err(Error::UnsupportedRange)
    } else {
        Ok(())
    }
}

fn add(sum: &mut f64, correction: &mut f64, term: f64) {
    let adjusted = term - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

// DLMF 6.6.1/2. The alternating series is only used for |x|<=1;
// the positive series for Ei has no large-term cancellation.
fn series(x: f64, limit: usize) -> Result<f64, Error> {
    let mut power = 1.0;
    let mut sum = 0.0;
    let mut correction = 0.0;
    for k in 1..=limit {
        power *= x / k as f64;
        let term = power / k as f64;
        add(&mut sum, &mut correction, term);
        if term.abs() <= 0.25 * f64::EPSILON * sum.abs() || term == 0.0 {
            return Ok((EULER + x.abs().ln()) + sum);
        }
    }
    Err(Error::NonConvergence)
}

// Integrate the Taylor series of exp(t)/t about Ei's zero (DLMF 6.13.1).
// The split root retains distances smaller than an ulp of the root.
fn near_root(x: f64, limit: usize) -> Result<f64, Error> {
    const ROOT_HI: f64 = 0.372_507_410_781_366_6;
    const ROOT_LO: f64 = 1.314_018_341_438_602_9e-17;
    let delta = (x - ROOT_HI) - ROOT_LO;
    let mut exponential = ROOT_HI.exp();
    let mut coefficient = exponential / ROOT_HI;
    let mut power = delta;
    let mut sum = coefficient * power;
    let mut correction = 0.0;
    for k in 1..=limit {
        exponential /= k as f64;
        coefficient = (exponential - coefficient) / ROOT_HI;
        power *= delta;
        let term = coefficient * power / (k + 1) as f64;
        add(&mut sum, &mut correction, term);
        if term.abs() <= 0.25 * f64::EPSILON * sum.abs() {
            return Ok(sum);
        }
    }
    Err(Error::NonConvergence)
}

// Scaled asymptotics: DLMF 6.12.1/2 and 8.20.1. Starting at 50
// permits termination at f64 precision well before the least term.
// sign=-1 gives E_n(x), sign=+1 gives Ei(x) (n=1) or -E2(-x) (n=2).
fn asymptotic(x: f64, n: u32, sign: f64, limit: usize) -> Result<f64, Error> {
    let mut term = 1.0;
    let mut sum = 1.0;
    let mut correction = 0.0;
    for k in 0..limit {
        let next = term * (sign * (k as f64 + f64::from(n)) / x);
        if next.abs() > term.abs() {
            return Err(Error::NonConvergence);
        }
        add(&mut sum, &mut correction, next);
        if next.abs() <= 0.25 * f64::EPSILON * sum.abs() {
            return Ok(sum / x);
        }
        term = next;
    }
    Err(Error::NonConvergence)
}

// Contracted continued fraction, DLMF 8.19.17 / 8.9.2, modified Lentz.
// Evaluate E2 directly, avoiding exp(-x)-x*E1(x) cancellation at large x.
fn fraction(n: u32, x: f64, limit: usize) -> Result<f64, Error> {
    let mut b = x + f64::from(n);
    let mut c = 1e300;
    let mut d = 1.0 / b;
    let mut value = d;
    for k in 1..=limit {
        let a = -(k as f64) * (k as f64 + f64::from(n) - 1.0);
        b += 2.0;
        d = 1.0 / guard(b + a * d);
        c = guard(b + a / c);
        let change = c * d;
        value *= change;
        if !value.is_finite() || value <= 0.0 {
            return Err(Error::NonConvergence);
        }
        if (change - 1.0).abs() <= 2.0 * f64::EPSILON {
            return Ok(value);
        }
    }
    Err(Error::NonConvergence)
}

fn positive_en_scaled(n: u32, x: f64) -> Result<f64, Error> {
    if x <= 1.0 {
        let e1 = -series(-x, LIMIT)?;
        return Ok(if n == 1 {
            x.exp() * e1
        } else {
            x.exp() * (-x).mul_add(e1, (-x).exp())
        });
    }
    if x >= ASYMPTOTIC {
        asymptotic(x, n, -1.0, LIMIT)
    } else {
        fraction(n, x, LIMIT)
    }
}

fn positive_ei_scaled(x: f64) -> Result<f64, Error> {
    if x >= ASYMPTOTIC {
        asymptotic(x, 1, 1.0, LIMIT)
    } else {
        let value = if (0.25..=0.5).contains(&x) {
            near_root(x, LIMIT)?
        } else {
            series(x, LIMIT)?
        };
        Ok((-x).exp() * value)
    }
}

// Split the exponential so exp(exponent) cannot overflow or prematurely
// underflow when its product with the scaled result is representable.
fn unscale(value: f64, exponent: f64) -> Result<f64, Error> {
    if exponent.abs() <= 700.0 {
        return Ok(value * exponent.exp());
    }
    if exponent.abs() < 1400.0 {
        let half = (0.5 * exponent).exp();
        let result = (value * half) * half;
        return if result.is_finite() {
            Ok(result)
        } else {
            Err(Error::Overflow)
        };
    }
    Ok(finite_exp(exponent + value.abs().ln())?.copysign(value))
}

/// Returns the real exponential integral `E1(x)`, unnormalized.
///
/// For positive x this is `integral_1^infinity exp(-xt)/t dt`; for negative
/// x it is the real continuation `-Ei(-x)`. All finite nonzero x are
/// supported; either signed zero and NaN return [`Error::Domain`].
/// Infinities, overflow, underflow and bounded convergence follow the
/// [module contract](self). See the [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#exponential-integrals).
pub fn e1(x: f64) -> Result<f64, Error> {
    validate(x, true)?;
    unscale(e1_scaled(x)?, -x)
}

/// Returns `exp(x)*E1(x)`, evaluated directly to retain large-argument tails.
///
/// The normalization, finite nonzero domain, poles and errors are those
/// of [`e1`]. See the [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#exponential-integrals).
pub fn e1_scaled(x: f64) -> Result<f64, Error> {
    validate(x, true)?;
    if x > 0.0 {
        positive_en_scaled(1, x)
    } else {
        Ok(-positive_ei_scaled(-x)?)
    }
}

/// Returns the real, unnormalized second exponential integral `E2(x)`.
///
/// For positive x it is `integral_1^infinity exp(-xt)/t^2 dt`; negative x
/// uses the real continuation `exp(-x)-x*E1(x)`. Returns its continuous
/// limit 1 at either signed zero. All finite x are supported. NaN returns
/// [`Error::Domain`]; other errors and underflow follow the [module contract](self).
/// See the [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#exponential-integrals).
pub fn e2(x: f64) -> Result<f64, Error> {
    validate(x, false)?;
    unscale(e2_scaled(x)?, -x)
}

/// Returns `exp(x)*E2(x)` directly, including when E2 overflows or underflows.
///
/// Has [`e2`]'s normalization, all-finite domain, exact value 1 at either
/// zero, and exceptional behavior. See the [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#exponential-integrals).
pub fn e2_scaled(x: f64) -> Result<f64, Error> {
    validate(x, false)?;
    if x == 0.0 {
        Ok(1.0)
    } else if x > 0.0 {
        positive_en_scaled(2, x)
    } else if -x >= ASYMPTOTIC {
        Ok(-asymptotic(-x, 2, 1.0, LIMIT)?)
    } else {
        Ok(x.mul_add(positive_ei_scaled(-x)?, 1.0))
    }
}

/// Returns the real `Ei(x)=PV integral_-infinity^x exp(t)/t dt`.
///
/// The principal value applies to positive x; negative x needs no PV.
/// Supports all finite nonzero x; either signed zero and NaN return
/// [`Error::Domain`]. Uses a root-centered series near its positive zero.
/// Infinities, overflow, underflow and convergence follow the [module contract](self).
/// See the [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#exponential-integrals).
pub fn ei(x: f64) -> Result<f64, Error> {
    validate(x, true)?;
    unscale(ei_scaled(x)?, x)
}

/// Returns `exp(-x)*Ei(x)` directly, retaining useful large-argument values.
///
/// Has [`ei`]'s real convention, all-finite nonzero domain, poles and errors.
/// See the [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#exponential-integrals).
pub fn ei_scaled(x: f64) -> Result<f64, Error> {
    validate(x, true)?;
    if x > 0.0 {
        positive_ei_scaled(x)
    } else {
        Ok(-positive_en_scaled(1, -x)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounded_iteration_failures() {
        assert_eq!(series(2.0, 1), Err(Error::NonConvergence));
        assert_eq!(fraction(1, 2.0, 1), Err(Error::NonConvergence));
        assert_eq!(fraction(2, 2.0, 1), Err(Error::NonConvergence));
        assert_eq!(asymptotic(50.0, 1, 1.0, 1), Err(Error::NonConvergence));
        assert_eq!(asymptotic(1.0, 2, 1.0, LIMIT), Err(Error::NonConvergence));
        assert_eq!(near_root(0.4, 1), Err(Error::NonConvergence));
    }
}

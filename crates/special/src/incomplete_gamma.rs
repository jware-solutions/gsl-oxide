//! Regularized gamma tails; [GSL definitions](https://www.gnu.org/software/gsl/doc/html/specfunc.html#incomplete-gamma-functions).

use crate::{
    gamma::{correction, ln_gamma1p, positive_log},
    log1pmx, Error, LN_SQRT_2_PI,
};

const LIMIT: usize = 32_768;
const MAX_SHAPE: f64 = 1_000_000.0;

// log(x^a exp(-x)/Gamma(a)); deviance form prevents catastrophic
// cancellation at the large-shape transition. DLMF 8.11.6, 5.11.1.
fn prefactor(a: f64, x: f64) -> f64 {
    if a >= 16.0 {
        let u = (x - a) / a;
        let deviance = if u.abs() < 0.5 {
            a * log1pmx(u)
        } else {
            a * (x.ln() - a.ln()) + a - x
        };
        deviance + 0.5 * a.ln() - LN_SQRT_2_PI - correction(a)
    } else {
        a * x.ln() - x - positive_log(a)
    }
}

fn tails(a: f64, x: f64, limit: usize) -> Result<(f64, f64), Error> {
    if !a.is_finite() || a <= 0.0 || x.is_nan() || x < 0.0 {
        return Err(Error::Domain);
    }
    if a > MAX_SHAPE {
        return Err(Error::UnsupportedRange);
    }
    if x == 0.0 {
        return Ok((0.0, 1.0));
    }
    if x == f64::INFINITY {
        return Ok((1.0, 0.0));
    }
    if a == 1.0 {
        return Ok((-(-x).exp_m1(), (-x).exp()));
    }
    // Small-shape/small-x upper tail: integrate the alternating exponential
    // series (DLMF 8.7.3), and use expm1 for 1-x^a/Gamma(a+1).
    // This retains Q even when P rounds to one.
    if a <= 0.25 && x <= a + 1.0 {
        let log = a * x.ln() - ln_gamma1p(a);
        let mut term = x;
        let mut sum = a * term / (a + 1.0);
        for n in 2..=limit {
            term *= -x / n as f64;
            let delta = a * term / (a + n as f64);
            sum += delta;
            if delta.abs() <= 2.0 * f64::EPSILON * sum.abs() {
                let q = -log.exp_m1() + log.exp() * sum;
                let p = if q > 0.5 {
                    log.exp() * (1.0 - sum)
                } else {
                    1.0 - q
                };
                return checked(p, q);
            }
        }
        return Err(Error::NonConvergence);
    }
    let log_pref = prefactor(a, x);
    if x < a + 1.0 {
        // DLMF 8.7.1: positive lower-tail series; decreasing ratios give
        // a geometric bound on the remaining terms, including x near a.
        let mut term = 1.0;
        let mut sum = 1.0;
        for n in 1..=limit {
            term *= x / (a + n as f64);
            sum += term;
            let ratio = x / (a + n as f64 + 1.0);
            if ratio < 1.0 && term * ratio / (1.0 - ratio) <= 2.0 * f64::EPSILON * sum {
                let log_p = if a < 0.25 {
                    a * x.ln() - x - ln_gamma1p(a) + sum.ln()
                } else {
                    log_pref - a.ln() + sum.ln()
                };
                return checked(log_p.exp(), -log_p.exp_m1());
            }
        }
    } else {
        // DLMF 8.9.2, modified Lentz evaluation. The selected region has
        // positive initial denominator. Guard tiny denominators throughout.
        let mut b = x + 1.0 - a;
        let mut c = 1.0e300;
        let mut d = 1.0 / b;
        let mut h = d;
        for n in 1..=limit {
            let an = -(n as f64) * (n as f64 - a);
            b += 2.0;
            d = guard(b + an * d);
            c = guard(b + an / c);
            d = 1.0 / d;
            let delta = c * d;
            h *= delta;
            if !h.is_finite() || h <= 0.0 {
                return Err(Error::NonConvergence);
            }
            if (delta - 1.0).abs() <= 4.0 * f64::EPSILON {
                let log_q = log_pref + h.ln();
                return checked(-log_q.exp_m1(), log_q.exp());
            }
        }
    }
    Err(Error::NonConvergence)
}

pub(crate) fn guard(x: f64) -> f64 {
    if x.abs() < 1.0e-300 {
        1.0e-300_f64.copysign(x)
    } else {
        x
    }
}

fn checked(p: f64, q: f64) -> Result<(f64, f64), Error> {
    if (0.0..=1.0).contains(&p) && (0.0..=1.0).contains(&q) {
        Ok((p, q))
    } else {
        Err(Error::NonConvergence)
    }
}

/// Returns `P(a,x)=integral_0^x t^(a-1)exp(-t)dt / Gamma(a)`.
///
/// Requires finite `0<a<=1e6` and `x>=0`, including positive infinity.
/// NaN, invalid shapes and negative x return [`Error::Domain`]; larger
/// shapes return [`Error::UnsupportedRange`]. Returns exactly 0 at either
/// signed zero and 1 at infinity. Tails may underflow to positive zero.
/// Bounded series/fraction failure returns [`Error::NonConvergence`].
/// Large transitions use a cancellation-free logarithmic prefactor and up
/// to 32768 terms.
/// See [DLMF §8.25](https://dlmf.nist.gov/8.25).
pub fn lower_regularized(a: f64, x: f64) -> Result<f64, Error> {
    Ok(tails(a, x, LIMIT)?.0)
}

/// Returns `Q(a,x)=integral_x^infinity t^(a-1)exp(-t)dt / Gamma(a)`.
///
/// Has [`lower_regularized`]'s domain, restrictions and errors, with exact values
/// 1 at zero and 0 at infinity. Computes the upper small tail directly; underflow
/// returns positive zero. See [DLMF §8.25](https://dlmf.nist.gov/8.25).
///
/// ```
/// use gsl_oxide_special::incomplete_gamma::{lower_regularized, upper_regularized};
/// assert!(upper_regularized(1.0, 100.0)? > 0.0);
/// assert_eq!(lower_regularized(1.0, 100.0)?, 1.0);
/// # Ok::<(), gsl_oxide_special::Error>(())
/// ```
pub fn upper_regularized(a: f64, x: f64) -> Result<f64, Error> {
    Ok(tails(a, x, LIMIT)?.1)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exhaustion_is_an_error() {
        for (a, x) in [(10.0, 9.0), (10.0, 12.0), (0.01, 0.5)] {
            assert_eq!(tails(a, x, 1), Err(Error::NonConvergence));
            assert!(tails(a, x, LIMIT).is_ok());
        }
    }
}

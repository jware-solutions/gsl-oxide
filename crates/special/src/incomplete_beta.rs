//! Regularized beta tails; [GSL definitions](https://www.gnu.org/software/gsl/doc/html/specfunc.html#incomplete-beta-function).

use crate::{
    beta::{ln_beta, log_a_beta},
    incomplete_gamma::guard,
    Error,
};
const LIMIT: usize = 16_384;

fn fraction(a: f64, b: f64, x: f64, limit: usize) -> Result<f64, Error> {
    // DLMF 8.17.22 continued fraction, two partial numerators per step.
    // Modified Lentz evaluation avoids unbounded numerator/denominator pairs.
    let mut c = 1.0;
    let mut d = 1.0 / guard(1.0 - (a + b) * x / (a + 1.0));
    let mut h = d;
    for m in 1..=limit {
        let m = m as f64;
        let even = m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m));
        let odd = -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0));
        let mut delta = 0.0;
        for numerator in [even, odd] {
            d = 1.0 / guard(1.0 + numerator * d);
            c = guard(1.0 + numerator / c);
            delta = d * c;
            h *= delta;
        }
        if !h.is_finite() || h <= 0.0 {
            return Err(Error::NonConvergence);
        }
        if (delta - 1.0).abs() <= 4.0 * f64::EPSILON {
            return Ok(h);
        }
    }
    Err(Error::NonConvergence)
}

fn tails(a: f64, b: f64, x: f64, limit: usize) -> Result<(f64, f64), Error> {
    if !a.is_finite() || !b.is_finite() || a <= 0.0 || b <= 0.0 || !(0.0..=1.0).contains(&x) {
        return Err(Error::Domain);
    }
    if a > 1e6 || b > 1e6 {
        return Err(Error::UnsupportedRange);
    }
    if x == 0.0 {
        return Ok((0.0, 1.0));
    }
    if x == 1.0 {
        return Ok((1.0, 0.0));
    }
    if a == 1.0 && b == 1.0 {
        return Ok((x, 1.0 - x));
    }
    if b == 1.0 {
        let log = a * x.ln();
        return Ok((x.powf(a), -log.exp_m1()));
    }
    if a == 1.0 {
        let log = b * (-x).ln_1p();
        return Ok((-log.exp_m1(), log.exp()));
    }
    if a == b && x == 0.5 {
        return Ok((0.5, 0.5));
    }
    // Preserve logs from the ORIGINAL x. Reconstructing x after symmetry
    // would discard small supplied x when 1-x rounds to one.
    let swap = x > (a + 1.0) / (a + b + 2.0);
    let (aa, bb, xx) = if swap { (b, a, 1.0 - x) } else { (a, b, x) };
    if aa <= 0.25 {
        // DLMF 8.17.7: x^a/(aB) [1+a sum (1-b)_n x^n/(n!(a+n))].
        // Region selection bounds b*x by about 1.25, avoiding cancellation
        // of large series terms. Evaluate 1-I via expm1, preserving O(a)
        // tails when the small shape makes I nearly one even in this region.
        let log_x = if swap { (-x).ln_1p() } else { x.ln() };
        let log = aa * log_x - log_a_beta(aa, bb);
        let mut term = 1.0;
        let mut sum = 0.0;
        for n in 1..=limit {
            term *= (n as f64 - bb) * xx / n as f64;
            let delta = term / (aa + n as f64);
            sum += delta;
            if delta == 0.0 || delta.abs() <= 2.0 * f64::EPSILON * sum.abs() {
                let upper = -log.exp_m1() - log.exp() * aa * sum;
                let lower = if upper > 0.5 {
                    log.exp() * (1.0 + aa * sum)
                } else {
                    1.0 - upper
                };
                if !(0.0..=1.0).contains(&upper) {
                    return Err(Error::NonConvergence);
                }
                return Ok(if swap { (upper, lower) } else { (lower, upper) });
            }
        }
        return Err(Error::NonConvergence);
    }
    let front = a * x.ln() + b * (-x).ln_1p() - ln_beta(a, b)?;
    let log = front + fraction(aa, bb, xx, limit)?.ln() - aa.ln();
    let small = log.exp();
    let large = -log.exp_m1();
    if !(0.0..=1.0).contains(&small) || !(0.0..=1.0).contains(&large) {
        return Err(Error::NonConvergence);
    }
    Ok(if swap { (large, small) } else { (small, large) })
}

/// Returns `I_x(a,b)=integral_0^x t^(a-1)(1-t)^(b-1)dt / B(a,b)`.
///
/// Requires finite `0<a,b<=1e6` and `0<=x<=1`. NaN, infinities and invalid
/// domains return [`Error::Domain`]; larger shapes return
/// [`Error::UnsupportedRange`]. Returns exact 0/1 at x=0/1 (after shape
/// validation); underflow returns positive zero. Symmetry and a bounded
/// continued fraction select a stable tail. Failure or numerical breakdown
/// returns [`Error::NonConvergence`]. See [DLMF §8.17](https://dlmf.nist.gov/8.17).
/// Log-factor cancellation can cost about 1e-9 absolute error at the largest
/// balanced shapes; full f64 relative accuracy is not promised there.
///
/// ```
/// use gsl_oxide_special::incomplete_beta::{regularized, complement};
/// assert!((regularized(2.0, 3.0, 0.5)? - 0.6875).abs() < 1e-13);
/// assert!(complement(0.5, 8.0, 0.999)? > 0.0);
/// # Ok::<(), gsl_oxide_special::Error>(())
/// ```
pub fn regularized(a: f64, b: f64, x: f64) -> Result<f64, Error> {
    Ok(tails(a, b, x, LIMIT)?.0)
}

/// Returns `1-I_x(a,b)`, evaluating the upper small tail directly.
///
/// Has [`regularized`]'s domain, restrictions and errors, with exact 1/0
/// at x=0/1. Subnormal tails are retained when representable and otherwise
/// underflow to positive zero. See [DLMF 8.17.4](https://dlmf.nist.gov/8.17.E4).
pub fn complement(a: f64, b: f64, x: f64) -> Result<f64, Error> {
    Ok(tails(a, b, x, LIMIT)?.1)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exhaustion_is_an_error() {
        assert_eq!(tails(30.0, 40.0, 0.4, 1), Err(Error::NonConvergence));
        assert!(tails(30.0, 40.0, 0.4, LIMIT).is_ok());
        assert_eq!(tails(0.1, 0.2, 0.3, 1), Err(Error::NonConvergence));
        assert!(tails(0.1, 0.2, 0.3, LIMIT).is_ok());
    }
}

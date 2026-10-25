//! Generalized Laguerre polynomials; [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#laguerre-functions).

use crate::Error;

const MAX_DEGREE: u32 = 1000;
const MAX_ALPHA: f64 = 1000.0;
const MAX_X: f64 = 10000.0;
const SCALE: f64 = 1.340_780_792_994_259_7e154; // 2^512, exact
const THRESHOLD: f64 = 3.273_390_607_896_142e150; // 2^500, exact

fn restore(mut value: f64, scales: u32) -> Result<f64, Error> {
    for _ in 0..scales {
        value *= SCALE;
        if !value.is_finite() {
            return Err(Error::Overflow);
        }
    }
    Ok(value)
}

/// Returns the unnormalized generalized Laguerre polynomial `L_n^alpha(x)`.
///
/// Uses GSL's normalization `(alpha+1)_n/n! * 1F1(-n,alpha+1,x)`;
/// alpha=0 gives the ordinary Laguerre polynomial. Requires finite
/// `-1<alpha<=1000`, finite `|x|<=10000`, and `n<=1000`.
/// Invalid alpha and non-finite parameters return [`Error::Domain`];
/// exceeded finite limits or degree return [`Error::UnsupportedRange`].
/// Validation precedes even the exact n=0 and n=1 shortcuts.
///
/// Uses at most n steps of the three-term recurrence (DLMF 18.9.1),
/// with a difference recurrence for `|x|<=1`, power-of-two scaling to
/// avoid intermediate overflow, and the binomial product at either signed
/// zero. Unrepresentable results return
/// [`Error::Overflow`], numerical breakdown [`Error::NonConvergence`];
/// ordinary underflow is allowed. Roundoff grows with degree and relative
/// accuracy deteriorates near roots; use absolute error allowances there.
/// See the [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#laguerre-functions).
pub fn generalized(n: u32, alpha: f64, x: f64) -> Result<f64, Error> {
    if !alpha.is_finite() || alpha <= -1.0 || !x.is_finite() {
        return Err(Error::Domain);
    }
    if n > MAX_DEGREE || alpha > MAX_ALPHA || x.abs() > MAX_X {
        return Err(Error::UnsupportedRange);
    }
    if n == 0 {
        return Ok(1.0);
    }
    if n == 1 {
        return Ok((alpha + 1.0) - x);
    }
    if x == 0.0 {
        // DLMF 18.5.12 at x=0. No subtractive cancellation when alpha
        // approaches -1; factors after the first are close to one.
        let mut value = alpha + 1.0;
        for k in 2..=n {
            value *= 1.0 + alpha / f64::from(k);
            if !value.is_finite() {
                return Err(Error::Overflow);
            }
        }
        return Ok(value);
    }
    let small_x = x.abs() <= 1.0;
    let mut previous = if small_x { alpha - x } else { 1.0 };
    let mut current = (alpha + 1.0) - x;
    let mut scales = 0;
    for k in 1..n {
        if current.abs().max(previous.abs()) > THRESHOLD {
            previous /= SCALE;
            current /= SCALE;
            scales += 1;
        }
        let k = f64::from(k);
        let next = if small_x {
            // d_k=L_k-L_(k-1). This avoids subtracting two nearly equal
            // O(1) terms at each step when n*x is small (DLMF 18.9.1).
            let difference = (k + alpha).mul_add(previous, -x * current) / (k + 1.0);
            previous = difference;
            current + difference
        } else {
            let a = (2.0 * k + alpha + 1.0 - x) / (k + 1.0);
            let b = (k + alpha) / (k + 1.0);
            let next = a.mul_add(current, -b * previous);
            previous = current;
            next
        };
        if !next.is_finite() {
            return Err(Error::NonConvergence);
        }
        current = next;
    }
    restore(current, scales)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn restoration_checks_overflow_and_signed_zero() {
        assert_eq!(restore(1.0, 2), Err(Error::Overflow));
        assert_eq!(restore(0.25, 2), Ok(2.0_f64.powi(1022)));
        assert_eq!(restore(-0.0, 2).unwrap().to_bits(), (-0.0_f64).to_bits());
    }
}

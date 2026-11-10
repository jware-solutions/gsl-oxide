//! Positive-shape beta functions; [GSL definitions](https://www.gnu.org/software/gsl/doc/html/specfunc.html#gamma-and-beta-functions).

use crate::{
    finite_exp,
    gamma::{correction, correction_from_inverse, log_ratio, positive_log},
    Error, LN_SQRT_2_PI,
};

// log(a B(a,b)) for 0<a<=1/4, retaining O(a) when this approaches zero.
/// Returns `ln(a * B(a,b))` for small positive `a` and positive `b`, preserving small logarithmic corrections.
pub(crate) fn log_a_beta(a: f64, mut b: f64) -> f64 {
    let mut shift = 0.0;
    while b < 16.0 {
        let ratio = a / b;
        shift += if ratio.is_finite() {
            ratio.ln_1p()
        } else {
            a.ln() - b.ln()
        };
        b += 1.0;
    }
    crate::gamma::ln_gamma1p(a) - crate::gamma::log_ratio(b, a) + shift
}

/// Returns `ln(B(a,b))`, where `B(a,b)=Gamma(a)Gamma(b)/Gamma(a+b)`.
///
/// Requires positive finite shapes; NaN, infinities and nonpositive shapes
/// return [`Error::Domain`]. A log outside f64 returns [`Error::Overflow`].
/// Uses beta recurrence and combined Stirling terms, retaining log-ratios
/// for unequal shapes without forming an overflowing `a+b` or subtracting
/// large log gammas. See [DLMF §5.12](https://dlmf.nist.gov/5.12).
/// The log remains available when [`beta`] underflows; no zeros or poles
/// occur in this domain. Rounding grows with the magnitude of the log.
pub fn ln_beta(mut a: f64, mut b: f64) -> Result<f64, Error> {
    if !a.is_finite() || !b.is_finite() || a <= 0.0 || b <= 0.0 {
        return Err(Error::Domain);
    }
    if a > b {
        std::mem::swap(&mut a, &mut b);
    }
    if a == 1.0 {
        return Ok(-b.ln());
    }
    if b == 1.0 {
        return Ok(-a.ln());
    }
    if a <= 0.25 {
        return Ok(log_a_beta(a, b) - a.ln());
    }
    if a < 16.0 && b >= 16.0 {
        return Ok(positive_log(a) - log_ratio(b, a));
    }
    let mut shift = 0.0;
    // B(a,b)=(1+b/a) B(a+1,b); at most 32 steps total.
    while a < 16.0 {
        shift += if b / a < f64::MAX {
            (b / a).ln_1p()
        } else {
            b.ln() - a.ln()
        };
        a += 1.0;
    }
    while b < 16.0 {
        shift += (a / b).ln_1p();
        b += 1.0;
    }
    if a > b {
        std::mem::swap(&mut a, &mut b);
    }
    let r = a / b;
    let l = r.ln_1p();
    let log_sum = b.ln() + l;
    let inv_sum = (1.0 / b) / (1.0 + r);
    let log_r = if r >= 0.5 { r.ln() } else { a.ln() - b.ln() };
    let value = (a - 0.5) * (log_r - l) - (b - 0.5) * l - 0.5 * log_sum
        + LN_SQRT_2_PI
        + correction(a)
        + correction(b)
        - correction_from_inverse(inv_sum)
        + shift;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Overflow)
    }
}

/// Returns the positive beta integral `integral_0^1 t^(a-1)(1-t)^(b-1) dt`.
///
/// Uses [`ln_beta`]'s positive finite domain and errors. Overflow returns
/// [`Error::Overflow`]. Representable subnormals are retained; rounded
/// underflow returns positive zero. Gamma factors are never formed.
/// See [DLMF 5.12.1](https://dlmf.nist.gov/5.12.E1).
pub fn beta(a: f64, b: f64) -> Result<f64, Error> {
    let log = ln_beta(a, b)?;
    if a == 1.0 || b == 1.0 {
        let value = 1.0 / if a == 1.0 { b } else { a };
        return if value.is_finite() {
            Ok(value)
        } else {
            Err(Error::Overflow)
        };
    }
    finite_exp(log)
}

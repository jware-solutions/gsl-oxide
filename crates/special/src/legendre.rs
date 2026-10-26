//! Ordinary Legendre polynomials; [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#legendre-polynomials).

use crate::Error;

const MAX_DEGREE: u32 = 1000;

/// Returns the ordinary, unnormalized Legendre polynomial `P_n(x)`.
///
/// Requires finite `-1<=x<=1` and `n<=1000`. Invalid x (including NaN
/// and infinities) returns [`Error::Domain`], even for n=0; larger n
/// returns [`Error::UnsupportedRange`]. Endpoints are exact:
/// `P_n(1)=1`, `P_n(-1)=(-1)^n`. Odd degrees retain x's signed zero.
///
/// Uses at most n recurrence steps (DLMF 18.9.1), with a difference
/// recurrence near the endpoints to reduce cancellation. Roundoff grows
/// with degree; use absolute error near roots. Subnormal inputs use the
/// leading Taylor term, multiplying by x only at the end. A numerical breakdown
/// returns [`Error::NonConvergence`]. See the
/// [GSL reference](https://www.gnu.org/software/gsl/doc/html/specfunc.html#legendre-polynomials).
pub fn polynomial(n: u32, x: f64) -> Result<f64, Error> {
    if !x.is_finite() || !(-1.0..=1.0).contains(&x) {
        return Err(Error::Domain);
    }
    if n > MAX_DEGREE {
        return Err(Error::UnsupportedRange);
    }
    if n == 0 {
        return Ok(1.0);
    }
    if n == 1 || (x == 0.0 && n % 2 == 1) {
        return Ok(x);
    }
    let t = x.abs();
    if t < f64::MIN_POSITIVE {
        // DLMF 18.5.7: P_(2m)(0)=(-1)^m binomial(2m,m)/4^m,
        // P_(2m+1)'(0)=(2m+1)P_(2m)(0). For subnormal x the omitted
        // quadratic/cubic terms are negligible. Multiply by x only once:
        // a degree recurrence on subnormals would round at every step.
        let mut central = 1.0;
        for k in 1..=n / 2 {
            central *= -f64::from(2 * k - 1) / f64::from(2 * k);
        }
        return Ok(if n % 2 == 0 {
            central
        } else {
            (f64::from(n) * central) * x
        });
    }
    let mut previous = 1.0;
    let mut current = t;
    if t > 0.5 {
        // Reinsch-style rearrangement: d_k=P_k-P_(k-1),
        // d_(k+1)=((2k+1)(t-1)P_k+k*d_k)/(k+1).
        let mut difference = t - 1.0;
        for k in 1..n {
            let k = f64::from(k);
            difference = ((2.0 * k + 1.0) * (t - 1.0)).mul_add(current, k * difference) / (k + 1.0);
            current += difference;
        }
    } else {
        for k in 1..n {
            let k = f64::from(k);
            let next = ((2.0 * k + 1.0) * t).mul_add(current, -k * previous) / (k + 1.0);
            previous = current;
            current = next;
        }
    }
    if !current.is_finite() {
        return Err(Error::NonConvergence);
    }
    Ok(if x.is_sign_negative() && n % 2 == 1 {
        -current
    } else {
        current
    })
}

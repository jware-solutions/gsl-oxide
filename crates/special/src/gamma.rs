//! Gamma and signed log gamma; [GSL definitions](https://www.gnu.org/software/gsl/doc/html/specfunc.html#gamma-and-beta-functions).

use crate::{finite_exp, Error, EULER, LN_SQRT_2_PI};

// B_(2k)/(2k(2k-1)), k=1..8, from the rational Bernoulli numbers.
// DLMF 5.11.1 and its first-omitted-term bound (5.11(ii)). At x>=16
// the omitted term is < 1e-21. These are formula coefficients, not a fit.
const STIRLING: [f64; 8] = [
    1.0 / 12.0,
    -1.0 / 360.0,
    1.0 / 1260.0,
    -1.0 / 1680.0,
    1.0 / 1188.0,
    -691.0 / 360360.0,
    1.0 / 156.0,
    -3617.0 / 122400.0,
];

pub(crate) fn correction(x: f64) -> f64 {
    correction_from_inverse(1.0 / x)
}

pub(crate) fn correction_from_inverse(inv: f64) -> f64 {
    let square = inv * inv;
    STIRLING.iter().rev().fold(0.0, |s, c| s * square + c) * inv
}

// Euler--Maclaurin evaluation of zeta(k), k>=2 (DLMF 25.2.9), only
// for the local log-gamma Taylor series. Bernoulli/factorial coefficients
// are obtained from the same rational Bernoulli numbers above.
fn zeta(k: i32) -> f64 {
    let s = f64::from(k);
    let mut sum = (1..16).map(|n| f64::from(n).powi(-k)).sum::<f64>();
    sum += 16.0_f64.powi(1 - k) / (s - 1.0) + 0.5 * 16.0_f64.powi(-k);
    let mut rising = s;
    let mut factorial = 1.0;
    for (j, c) in STIRLING.iter().enumerate() {
        let m = 2 * j + 2;
        if j > 0 {
            rising *= (s + (m - 3) as f64) * (s + (m - 2) as f64);
        }
        // B_m/m! = STIRLING[j]/(m-2)!
        if j > 0 {
            factorial *= ((m - 3) * (m - 2)) as f64;
        }
        sum += c / factorial * rising * 16.0_f64.powi(-k - m as i32 + 1);
    }
    sum
}

// log Gamma(b+a)-log Gamma(b), b>=16, 0<a<=16.
// Algebraic difference of DLMF 5.11.1, preserving O(a) terms even when
// b+a rounds to b. Used by both complete and incomplete beta.
pub(crate) fn log_ratio(b: f64, a: f64) -> f64 {
    let r = a / b;
    let l = r.ln_1p();
    let mut difference = 0.0;
    let inv = 1.0 / b;
    let mut power = inv;
    for (j, c) in STIRLING.iter().enumerate() {
        difference += c * power * (-((2 * j + 1) as f64) * l).exp_m1();
        power *= inv * inv;
    }
    a * b.ln() + a * ((a - 0.5) / b) + (b + a - 0.5) * crate::log1pmx(r) + difference
}

pub(crate) fn ln_gamma1p(z: f64) -> f64 {
    // DLMF 5.7.3: log Gamma(1+z)=-Euler*z+sum (-z)^k*zeta(k)/k.
    // |z|<=1/4; the remainder after 32 terms is < 2e-21.
    let mut power = z * z;
    let mut sum = -EULER * z;
    for k in 2..=32 {
        sum += power * zeta(k) / f64::from(k);
        power *= -z;
    }
    sum
}

pub(crate) fn positive_log(mut x: f64) -> f64 {
    if (0.75..=1.25).contains(&x) {
        return ln_gamma1p(x - 1.0);
    }
    if (1.75..=2.25).contains(&x) {
        let z = x - 2.0;
        return z.ln_1p() + ln_gamma1p(z);
    }
    if x < 0.25 {
        return ln_gamma1p(x) - x.ln();
    }
    let mut shift = 0.0;
    // At most sixteen recurrence steps.
    while x < 16.0 {
        shift += x.ln();
        x += 1.0;
    }
    (x - 0.5) * x.ln() - x + LN_SQRT_2_PI + correction(x) - shift
}

/// Returns `(ln(abs(Gamma(x))), sign)`, with sign `+1` or `-1`.
///
/// All finite nonpole real inputs are supported. Zero (either sign), negative
/// integers, NaN and either infinity return [`Error::Domain`]. Logarithmic
/// overflow returns [`Error::Overflow`]; a finite log survives gamma overflow
/// and underflow. Uses recurrence, reflection (DLMF 5.5.3), and
/// [Stirling's expansion](https://dlmf.nist.gov/5.11).
/// Ordinary f64 accuracy is intended, not correct rounding; small log values
/// outside the local Taylor regions can have about 1e-14 absolute roundoff.
///
/// ```
/// let (log, sign) = gsl_oxide_special::gamma::ln_abs_gamma_with_sign(-0.5)?;
/// assert_eq!(sign, -1);
/// assert!((log - (2.0 * std::f64::consts::PI.sqrt()).ln()).abs() < 1e-14);
/// # Ok::<(), gsl_oxide_special::Error>(())
/// ```
pub fn ln_abs_gamma_with_sign(x: f64) -> Result<(f64, i32), Error> {
    if !x.is_finite() || (x <= 0.0 && x == x.trunc()) {
        return Err(Error::Domain);
    }
    let (log, sign) = if x > 0.0 {
        (positive_log(x), 1)
    } else if x > -0.25 {
        (ln_gamma1p(x) - (-x).ln(), -1)
    } else {
        // Reduce before multiplying by pi. All representable negative
        // nonintegers have |x|<2^52; nearest-integer reduction is exact.
        let nearest = x.round();
        let r = x - nearest;
        let sine = (std::f64::consts::PI * r).sin();
        let parity = if nearest % 2.0 == 0.0 { 1 } else { -1 };
        let sign = if sine > 0.0 { parity } else { -parity };
        // Gamma(1-x)=(-x)*Gamma(-x), avoiding loss of x when 1-x rounds.
        (
            std::f64::consts::PI.ln() - sine.abs().ln() - (-x).ln() - positive_log(-x),
            sign,
        )
    };
    if log.is_finite() {
        Ok((log, sign))
    } else {
        Err(Error::Overflow)
    }
}

/// Returns `Gamma(x) = integral_0^infinity t^(x-1) exp(-t) dt` for positive
/// `x`, with real continuation to negative nonintegers.
///
/// Has the same argument and pole contract as [`ln_abs_gamma_with_sign`].
/// Unrepresentable magnitudes return [`Error::Overflow`]; underflow returns
/// signed zero.
/// The signed logarithm remains useful outside the representable gamma range.
/// See [DLMF §5.11](https://dlmf.nist.gov/5.11).
pub fn gamma(x: f64) -> Result<f64, Error> {
    let (log, sign) = ln_abs_gamma_with_sign(x)?;
    Ok(finite_exp(log)?.copysign(f64::from(sign)))
}

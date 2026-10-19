//! Gaussian error functions.

/// Returns the Gaussian error function of the real argument `x`:
/// `erf(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t*t) dt`.
///
/// All finite `f64` inputs are supported; results lie in `[-1, 1]`.
/// NaN returns NaN, infinities return their signed unit limits, and either
/// signed zero is preserved. Ordinary floating-point semantics suffice, so
/// the value is returned directly rather than through `Result`.
///
/// A convergent series avoids subtractive cancellation. Results have ordinary
/// `f64` rounding error; correct rounding is not guaranteed. For `|x| >= 6`,
/// the result is the signed unit limit because the remaining tail is smaller
/// than half the spacing immediately below 1.
///
/// Reference: [GSL error function](https://www.gnu.org/software/gsl/doc/html/specfunc.html#error-function).
///
/// ```
/// use gsl_oxide_special::error_functions::erf;
///
/// assert!((erf(1.0) - 0.8427007929497149).abs() < 2.0e-15);
/// ```
pub fn erf(x: f64) -> f64 {
    if x.is_nan() || x == 0.0 {
        return x;
    }

    let magnitude = x.abs();
    if magnitude >= 6.0 {
        // For x > 0, the Gaussian tail is bounded by exp(-x*x)/(x*sqrt(pi)).
        // Check before squaring, including for infinite or very large inputs.
        return 1.0_f64.copysign(x);
    }

    // Original evaluation of the mathematical series in DLMF 7.6.2:
    // https://dlmf.nist.gov/7.6.E2
    // Coefficients follow the recurrence from erf'(x) = 2*exp(-x*x)/sqrt(pi);
    // no third-party code, fitted coefficients or coefficient tables are used.
    let squared = magnitude * magnitude;
    let mut term = 1.0;
    let mut sum = 1.0;
    let mut correction = 0.0;
    let mut denominator = 3.0;
    // With squared < 36, 160 terms put the remainder far below f64 precision.
    for _ in 0..160 {
        term *= 2.0 * squared / denominator;
        let adjusted = term - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
        if term <= f64::EPSILON * sum {
            break;
        }
        denominator += 2.0;
    }

    // Multiply the argument last so subnormal inputs retain their scale.
    let value = magnitude * (std::f64::consts::FRAC_2_SQRT_PI * (-squared).exp() * sum);
    value.min(1.0).copysign(x)
}

//! Integer-order cylindrical Bessel functions, using the
//! [GSL definitions](https://www.gnu.org/software/gsl/doc/html/specfunc.html#bessel-functions).
//!
//! `J_n(x) = sum_k (-1)^k (x/2)^(n+2k)/(k!(n+k)!)`; I uses the
//! same series with positive terms. Y is the real irregular solution
//! `lim_(v->n) [cos(pi*v)J_v(x)-J_(-v)(x)]/sin(pi*v)`.
//! `K_n(x) = integral_0^infinity exp(-x*cosh(t))*cosh(n*t) dt`.
//! Scaled I/K multiply these definitions by `exp(-|x|)` / `exp(x)`.
//!
//! All entry points return `Result<f64, Error>`. Supported orders are
//! `0..=1000`; finite arguments satisfy `|x|<=1e6` for J/I and `0<x<=1e6`
//! for Y/K. NaN and invalid real domains return `Domain`; other infinities and
//! exceeded finite/order limits return `UnsupportedRange`. Y and K reject
//! zero and negative infinity (their limits at zero are negative/positive
//! infinity respectively).
//! J/I use integer parity on negative arguments: order zero is 1 at zero,
//! other orders are zero, with a negative zero for odd orders at -0.
//! Overflow is checked, underflow returns signed zero for J/I and positive
//! zero for K. Bounded kernel failure returns `NonConvergence`.
//!
//! Small series, local Taylor propagation of the Bessel ODE, and large-x
//! expansions provide the base orders. Recurrence follows
//! [DLMF §10.74](https://dlmf.nist.gov/10.74): J backward when order dominates,
//! I backward, Y/K forward, with normalization and rescaling. Accuracy is
//! ordinary f64, with absolute rather than relative accuracy near J/Y zeros;
//! accumulated recurrence roundoff grows with order (roughly 1e-11 relative
//! allowances are appropriate at order 1000 away from zeros).

use crate::{finite_exp, gamma::positive_log, Error, EULER};
use std::f64::consts::{FRAC_1_SQRT_2, PI};

const MAX_ORDER: u32 = 1000;
const MAX_X: f64 = 1e6;

/// Validates order `n` and argument `x`; `positive` selects the strictly positive domain. Returns the argument magnitude or an error.
fn validate(n: u32, x: f64, positive: bool) -> Result<f64, Error> {
    if x.is_nan() || (positive && x <= 0.0) {
        return Err(Error::Domain);
    }
    if n > MAX_ORDER || x.abs() > MAX_X {
        return Err(Error::UnsupportedRange);
    }
    Ok(x.abs())
}

/// Applies integer-order parity to `value` for argument `x` and order `n`.
fn parity(n: u32, x: f64, value: f64) -> f64 {
    if n % 2 == 1 && x.is_sign_negative() {
        -value
    } else {
        value
    }
}

// DLMF 10.2.2/10.25.2. For x<=1 this series is well conditioned for
// every supported order. Logarithmic leading factors preserve tiny results.
/// Evaluates the small-argument series for order `n` and argument `x`; `modified` selects I instead of J.
fn small(n: u32, x: f64, modified: bool) -> f64 {
    if x == 0.0 {
        return if n == 0 { 1.0 } else { 0.0 };
    }
    if n == 1 && x < 1e-150 {
        return x / 2.0;
    }
    let mut term = 1.0;
    let mut sum = 1.0;
    let factor = if modified { x * x / 4.0 } else { -x * x / 4.0 };
    for k in 1..=64 {
        term *= factor / (f64::from(k) * f64::from(n + k));
        sum += term;
        if term.abs() <= f64::EPSILON * sum.abs() {
            break;
        }
    }
    if n == 0 {
        sum
    } else if n == 1 {
        x / 2.0 * sum
    } else {
        (f64::from(n) * (x.ln() - std::f64::consts::LN_2) - positive_log(f64::from(n) + 1.0)
            + sum.ln())
        .exp()
    }
}

// DLMF 10.8.2 / 10.31.2 and derivatives Y1=-Y0', K1=-K0'.
// Harmonic-number coefficients follow directly from those formulas.
/// Returns the small-argument order-zero and order-one pair at `x`; `modified` selects K instead of Y.
fn singular_small(x: f64, modified: bool) -> (f64, f64) {
    let l = x.ln() - std::f64::consts::LN_2 + EULER;
    let f0 = small(0, x, modified);
    let f1 = small(1, x, modified);
    if x < 1e-150 {
        return if modified {
            (-l, 1.0 / x)
        } else {
            (2.0 / PI * l, -(2.0 / PI) / x)
        };
    }
    let mut term = 1.0;
    let mut harmonic = 0.0;
    let mut sum = 0.0;
    let mut derivative = 0.0;
    for k in 1..=64 {
        let k = f64::from(k);
        harmonic += 1.0 / k;
        term *= (if modified { 1.0 } else { -1.0 }) * (x * x / 4.0) / (k * k);
        let delta = harmonic * term;
        sum += delta;
        derivative += 2.0 * k * delta / x;
        if delta.abs() <= f64::EPSILON * sum.abs() {
            break;
        }
    }
    if modified {
        (-l * f0 + sum, f0 / x + l * f1 - derivative)
    } else {
        (
            2.0 / PI * (l * f0 - sum),
            2.0 / PI * (-f0 / x + l * f1 + derivative),
        )
    }
}

// Original Taylor stepping of x²y''+xy'+x²y=0, DLMF 10.2.1.
// Expanding around a positive center gives the coefficient recurrence below.
// Steps <=1/4 from x>=1 have analytic radius >=1, so 32 terms give a
// negligible truncation error. At most 124 steps, avoiding long alternating
// series cancellation. Both J0 and Y0 are stable oscillatory ODE solutions.
/// Advances a Bessel solution and its derivative from `center` to `target`, returning the order-zero and order-one values.
fn march(mut center: f64, target: f64, mut value: f64, mut derivative: f64) -> (f64, f64) {
    for _ in 0..124 {
        if center >= target {
            break;
        }
        let h = (target - center).min(0.25);
        let mut c = [0.0; 34];
        c[0] = value;
        c[1] = derivative;
        for k in 0..32 {
            let kf = k as f64;
            let previous = if k >= 1 { c[k - 1] } else { 0.0 };
            let previous2 = if k >= 2 { c[k - 2] } else { 0.0 };
            c[k + 2] = -(center * (kf + 1.0) * (2.0 * kf + 1.0) * c[k + 1]
                + (kf * kf + center * center) * c[k]
                + 2.0 * center * previous
                + previous2)
                / (center * center * (kf + 2.0) * (kf + 1.0));
        }
        value = c.iter().rev().fold(0.0, |s, coef| s * h + coef);
        derivative = (1..c.len()).rev().fold(0.0, |s, k| s * h + k as f64 * c[k]);
        center += h;
    }
    (value, -derivative)
}

// Hankel expansion DLMF 10.17.1, 10.17.3--4. Coefficients are products
// (4n²-(2k-1)²)/(8k), generated here. Only base orders and x>=32 are used;
// stop at small terms or before the least term starts increasing.
/// Evaluates the large-`x` asymptotic expansion for order `n`; `second` selects Y instead of J.
fn oscillatory_asym(n: u32, x: f64, second: bool) -> f64 {
    let mut even = 1.0;
    let mut odd = 0.0;
    let mut term = 1.0_f64;
    for k in 1..=128 {
        let kf = f64::from(k);
        let next = term * (4.0 * f64::from(n * n) - (2.0 * kf - 1.0).powi(2)) / (8.0 * kf * x);
        if next.abs() > term.abs() {
            break;
        }
        term = next;
        let signed = if k % 4 == 1 || k % 4 == 0 {
            term
        } else {
            -term
        };
        if k % 2 == 0 {
            even += signed;
        } else {
            odd += signed;
        }
        if term.abs() < f64::EPSILON * 0.25 {
            break;
        }
    }
    // Form sin/cos(x-pi/4-n*pi/2) without subtracting a small phase from x.
    let (sx, cx) = x.sin_cos();
    let (sin, cos) = if n == 0 {
        ((sx - cx) * FRAC_1_SQRT_2, (cx + sx) * FRAC_1_SQRT_2)
    } else {
        (-(sx + cx) * FRAC_1_SQRT_2, (sx - cx) * FRAC_1_SQRT_2)
    };
    (2.0 / (PI * x)).sqrt()
        * if second {
            sin * even + cos * odd
        } else {
            cos * even - sin * odd
        }
}

/// Returns the order-zero and order-one J or Y values at positive `x`, selected by `second`.
fn jy_base(x: f64, second: bool) -> (f64, f64) {
    if x <= 1.0 {
        return if second {
            singular_small(x, false)
        } else {
            (small(0, x, false), small(1, x, false))
        };
    }
    if x >= 32.0 {
        return (
            oscillatory_asym(0, x, second),
            oscillatory_asym(1, x, second),
        );
    }
    let (v0, v1) = if second {
        singular_small(1.0, false)
    } else {
        (small(0, 1.0, false), small(1, 1.0, false))
    };
    march(1.0, x, v0, -v1)
}

// DLMF 10.40.1--2, same product coefficients as the Hankel expansion.
/// Returns the exponentially scaled large-`x` I or K asymptotic for order `n`; `second` selects K.
fn modified_asym(n: u32, x: f64, second: bool) -> f64 {
    let mut term = 1.0_f64;
    let mut sum = 1.0;
    for k in 1..=128 {
        let kf = f64::from(k);
        let next = term * (4.0 * f64::from(n * n) - (2.0 * kf - 1.0).powi(2)) / (8.0 * kf * x)
            * if second { 1.0 } else { -1.0 };
        if next.abs() > term.abs() {
            break;
        }
        term = next;
        sum += term;
        if term.abs() < f64::EPSILON * 0.25 {
            break;
        }
    }
    sum * if second {
        (PI / (2.0 * x)).sqrt()
    } else {
        1.0 / (2.0 * PI * x).sqrt()
    }
}

/// Returns exponentially scaled I of base order `n` at nonnegative `x`.
fn i_base(n: u32, x: f64) -> f64 {
    if x >= 32.0 {
        return modified_asym(n, x, false);
    }
    // Positive series; at x<=32 all terms are safely representable, and
    // the term after 128 steps is <1e-110 of the sum. Thus the fixed limit
    // exceeds the convergence work throughout this base-order region.
    let mut term = if n == 0 { 1.0 } else { x / 2.0 };
    let mut sum = term;
    for k in 1..=128 {
        term *= (x * x / 4.0) / (f64::from(k) * f64::from(n + k));
        sum += term;
        if term <= f64::EPSILON * sum {
            break;
        }
    }
    (-x).exp() * sum
}

// K_n(x)e^x = integral_0^infinity exp(-x(cosh(t)-1))*cosh(n*t)dt,
// DLMF 10.32.9, n=0,1. The even analytic integrand permits trapezoidal
// refinement. Truncate at x(cosh(t)-1)>=48; beyond that the integral is
// < 1e-19 of the base value for 1<=x<=32. This is a private K kernel,
// not a general quadrature API. At most 7 refinements and 4096 samples each.
/// Computes scaled K0 and K1 at `x` by at most `refinements` quadrature refinements; reports nonconvergence.
fn k_integral(x: f64, refinements: usize) -> Result<(f64, f64), Error> {
    let mut previous = (0.0_f64, 0.0_f64);
    for level in 0..refinements {
        let h = 0.25 / 2.0_f64.powi(level as i32);
        let mut s0 = 0.5;
        let mut s1 = 0.5;
        let mut ended = false;
        for k in 1..=4096 {
            let t = k as f64 * h;
            let c = t.cosh();
            // cosh(t)-1 = 2sinh(t/2)^2 avoids near-zero cancellation.
            let exponent = x * 2.0 * (t / 2.0).sinh().powi(2);
            let v = (-exponent).exp();
            s0 += v;
            s1 += v * c;
            if exponent >= 48.0 {
                ended = true;
                break;
            }
        }
        if !ended {
            return Err(Error::NonConvergence);
        }
        let value = (h * s0, h * s1);
        if level > 0
            && (value.0 - previous.0).abs() <= 8.0 * f64::EPSILON * value.0
            && (value.1 - previous.1).abs() <= 8.0 * f64::EPSILON * value.1
        {
            return Ok(value);
        }
        previous = value;
    }
    Err(Error::NonConvergence)
}

/// Returns exponentially scaled K0 and K1 at positive `x`, or an error if quadrature fails.
fn k_base(x: f64) -> Result<(f64, f64), Error> {
    if x <= 1.0 {
        let (v0, v1) = singular_small(x, true);
        return Ok((v0 * x.exp(), v1 * x.exp()));
    }
    if x >= 32.0 {
        return Ok((modified_asym(0, x, true), modified_asym(1, x, true)));
    }
    k_integral(x, 7)
}

// Miller recurrence, normalized with the larger of J0/J1.
// Rescale the saved target along with both recurrence values.
/// Evaluates J of order `n` at `x` using a backward recurrence starting at `start` and normalizing against base orders.
fn j_backward(n: u32, x: f64, start: u32) -> f64 {
    let mut next = 0.0;
    let mut current = 1.0;
    let mut saved = 0.0;
    let mut order_one = 0.0;
    for k in (1..=start).rev() {
        if k == n {
            saved = current;
        }
        let previous = 2.0 * f64::from(k) / x * current - next;
        next = current;
        current = previous;
        if current.abs() > 1e200 {
            current *= 1e-200;
            next *= 1e-200;
            saved *= 1e-200;
        }
        if k == 1 {
            order_one = next;
        }
    }
    let (j0, j1) = jy_base(x, false);
    if j0.abs() >= j1.abs() {
        saved / current * j0
    } else {
        saved / order_one * j1
    }
}

/// Returns J of order `n` at `x`, allowing `refinements` recurrence checks before reporting nonconvergence.
fn j_value(n: u32, x: f64, refinements: usize) -> Result<f64, Error> {
    if x <= 1.0 {
        return Ok(small(n, x, false));
    }
    let (v0, v1) = jy_base(x, false);
    if n == 0 {
        return Ok(v0);
    }
    if n == 1 {
        return Ok(v1);
    }
    if f64::from(n) > x {
        let start = n + 64 + (40.0 * f64::from(n)).sqrt() as u32;
        let mut previous = j_backward(n, x, start);
        for extra in 1..=refinements {
            let value = j_backward(n, x, start + 32 * extra as u32);
            if value.is_finite() && (value - previous).abs() <= 32.0 * f64::EPSILON * value.abs() {
                return Ok(value);
            }
            previous = value;
        }
        return Err(Error::NonConvergence);
    }
    let (mut a, mut b) = (v0, v1);
    for k in 1..n {
        let c = 2.0 * f64::from(k) / x * b - a;
        a = b;
        b = c;
    }
    Ok(b)
}

/// Returns the logarithm of exponentially scaled I of order `n` at `x`, or a convergence error.
fn i_log(n: u32, x: f64) -> Result<f64, Error> {
    i_log_kernel(n, x, 8)
}

/// Evaluates the logarithm of scaled I for `n, x`, limiting backward-recurrence convergence checks to `refinements`.
fn i_log_kernel(n: u32, x: f64, refinements: usize) -> Result<f64, Error> {
    if n <= 1 {
        return Ok(i_base(n, x).ln());
    }
    if x <= 1.0 {
        // Avoid losing a representable unscaled value if its scaled value
        // underflows; in this region exp(x) is bounded and n is finite.
        let mut term = 1.0;
        let mut sum = 1.0;
        for k in 1..=64 {
            term *= (x * x / 4.0) / (f64::from(k) * f64::from(n + k));
            sum += term;
            if term <= f64::EPSILON * sum {
                break;
            }
        }
        return Ok(f64::from(n) * (x.ln() - std::f64::consts::LN_2)
            - positive_log(f64::from(n) + 1.0)
            + sum.ln()
            - x);
    }
    let start = n + 64 + (80.0 * x).sqrt().ceil() as u32;
    let mut previous = i_backward_log(n, x, start);
    for extra in 1..=refinements {
        let value = i_backward_log(n, x, start + 32 * extra as u32);
        if (value - previous).abs() <= 16.0 * f64::EPSILON * value.abs().max(1.0) {
            return Ok(value);
        }
        previous = value;
    }
    Err(Error::NonConvergence)
}

// Logarithmic Miller ratios avoid losing I_n when exp(-x)I_n underflows
// but I_n itself remains representable (e.g. n=1000, x=400).
/// Returns log scaled I for order `n` at `x`, using backward recurrence beginning at `start`.
fn i_backward_log(n: u32, x: f64, start: u32) -> f64 {
    let mut ratio = 0.0;
    let mut log = i_base(0, x).ln();
    let mut correction = 0.0;
    for k in (1..=start).rev() {
        ratio = 1.0 / (2.0 * f64::from(k) / x + ratio);
        if k <= n {
            let term = ratio.ln() - correction;
            let next = log + term;
            correction = (next - log) - term;
            log = next;
        }
    }
    log
}

// Dominant-solution forward recurrence, with logarithmic scale tracking.
/// Advances a dominant Bessel recurrence to order `n` at `x` from `first, second`; `modified` selects K. Returns the scaled value and its logarithmic scale.
fn forward_log(
    n: u32,
    x: f64,
    first: f64,
    second: f64,
    modified: bool,
) -> Result<(f64, f64), Error> {
    if n == 0 {
        return Ok((first, 0.0));
    }
    if n == 1 {
        return Ok((second, 0.0));
    }
    let (mut a, mut b) = (first, second);
    let mut log_scale = 0.0;
    for k in 1..n {
        // Scale BEFORE multiplying by 2k/x, including very small x.
        let factor_log = (2.0 * f64::from(k)).ln() - x.ln();
        let magnitude = a.abs().max(b.abs());
        if magnitude.ln() + factor_log > 460.0 {
            let scale = magnitude;
            a /= scale;
            b /= scale;
            log_scale += scale.ln();
            if factor_log > 460.0 {
                // Once the dominant small-x solution exceeds f64 it cannot
                // return to range; K unscaling may recover it, so retain logs.
                let log = b.abs().ln() + log_scale + factor_log;
                let mut result = log;
                for j in k + 1..n {
                    result += (2.0 * f64::from(j)).ln() - x.ln();
                }
                return Ok((if b < 0.0 { -1.0 } else { 1.0 }, result));
            }
        }
        let c = (2.0 * f64::from(k) / x) * b + if modified { a } else { -a };
        a = b;
        b = c;
        if !b.is_finite() {
            return Err(Error::Overflow);
        }
    }
    Ok((b, log_scale))
}

/// Returns cylindrical `J_0(x)` (regular solution of Bessel's equation).
/// See the [module contract](self) for domain, zero, limits and errors.
pub fn j0(x: f64) -> Result<f64, Error> {
    jn(0, x)
}
/// Returns cylindrical `J_1(x)`; odd parity, including signed zero.
/// See the [module contract](self) for domain, limits and errors.
pub fn j1(x: f64) -> Result<f64, Error> {
    jn(1, x)
}
/// Returns cylindrical `J_n(x)`, the solution regular at zero.
/// Nonnegative orders `n<=1000`, `|x|<=1e6`; see the [module contract](self).
pub fn jn(n: u32, x: f64) -> Result<f64, Error> {
    let a = validate(n, x, false)?;
    Ok(parity(n, x, j_value(n, a, 8)?))
}
/// Returns cylindrical `Y_0(x)` (irregular Bessel solution).
/// Positive domain and singularity/error behavior: [module contract](self).
pub fn y0(x: f64) -> Result<f64, Error> {
    yn(0, x)
}
/// Returns cylindrical `Y_1(x)`; see the [module contract](self).
pub fn y1(x: f64) -> Result<f64, Error> {
    yn(1, x)
}
/// Returns cylindrical `Y_n(x)` for `n<=1000`, `0<x<=1e6`.
/// Zero/negative/NaN inputs are domain errors; see the [module contract](self).
pub fn yn(n: u32, x: f64) -> Result<f64, Error> {
    validate(n, x, true)?;
    let (a, b) = jy_base(x, true);
    if n == 0 {
        return Ok(a);
    }
    if !b.is_finite() {
        return Err(Error::Overflow);
    }
    let (value, scale) = forward_log(n, x, a, b, false)?;
    Ok(finite_exp(value.abs().ln() + scale)?.copysign(value))
}
/// Returns `I_0(x)`, the regular solution of the modified Bessel equation.
/// See the [module contract](self); overflow is an error.
pub fn i0(x: f64) -> Result<f64, Error> {
    i_n(0, x)
}
/// Returns `I_1(x)`; odd parity. See the [module contract](self).
pub fn i1(x: f64) -> Result<f64, Error> {
    i_n(1, x)
}
/// Returns `I_n(x)` for `n<=1000`, `|x|<=1e6`.
/// Overflow is an error; underflow retains parity. See the [module contract](self).
pub fn i_n(n: u32, x: f64) -> Result<f64, Error> {
    i_evaluate(n, x, false)
}
/// Returns `exp(-|x|) I_0(x)` directly; see the [module contract](self).
pub fn i0_scaled(x: f64) -> Result<f64, Error> {
    i_n_scaled(0, x)
}
/// Returns `exp(-|x|) I_1(x)` directly; see the [module contract](self).
pub fn i1_scaled(x: f64) -> Result<f64, Error> {
    i_n_scaled(1, x)
}
/// Returns `exp(-|x|) I_n(x)` for `n<=1000`, `|x|<=1e6`.
/// Avoids unscaled overflow. Same parity, zero and errors as [`i_n`].
///
/// ```
/// use gsl_oxide_special::bessel::{i0, i0_scaled};
/// assert!(i0(1000.0).is_err());
/// assert!(i0_scaled(1000.0)? > 0.01);
/// # Ok::<(), gsl_oxide_special::Error>(())
/// ```
pub fn i_n_scaled(n: u32, x: f64) -> Result<f64, Error> {
    i_evaluate(n, x, true)
}
/// Evaluates I of order `n` at `x`; `scaled` selects exponential scaling and range or convergence failures are errors.
fn i_evaluate(n: u32, x: f64, scaled: bool) -> Result<f64, Error> {
    let a = validate(n, x, false)?;
    if a == 0.0 {
        return Ok(parity(n, x, if n == 0 { 1.0 } else { 0.0 }));
    }
    if a <= 1.0 && n <= 1 {
        let value = small(n, a, true) * if scaled { (-a).exp() } else { 1.0 };
        return Ok(parity(n, x, value));
    }
    let log = i_log(n, a)? + if scaled { 0.0 } else { a };
    Ok(parity(n, x, finite_exp(log)?))
}
/// Returns the decaying modified Bessel solution `K_0(x)`.
/// See the [module contract](self) for positive domain, singularity and errors.
pub fn k0(x: f64) -> Result<f64, Error> {
    kn(0, x)
}
/// Returns `K_1(x)`; see the [module contract](self).
pub fn k1(x: f64) -> Result<f64, Error> {
    kn(1, x)
}
/// Returns `K_n(x)` for `n<=1000`, `0<x<=1e6`.
/// Underflow returns positive zero; overflow is checked. See the [module contract](self).
pub fn kn(n: u32, x: f64) -> Result<f64, Error> {
    k_evaluate(n, x, false)
}
/// Returns `exp(x) K_0(x)` directly; see the [module contract](self).
pub fn k0_scaled(x: f64) -> Result<f64, Error> {
    kn_scaled(0, x)
}
/// Returns `exp(x) K_1(x)` directly; see the [module contract](self).
pub fn k1_scaled(x: f64) -> Result<f64, Error> {
    kn_scaled(1, x)
}
/// Returns `exp(x) K_n(x)` for `n<=1000`, `0<x<=1e6`.
/// Avoids unscaled underflow; otherwise shares [`kn`]'s domain and errors.
pub fn kn_scaled(n: u32, x: f64) -> Result<f64, Error> {
    k_evaluate(n, x, true)
}
/// Evaluates K of order `n` at positive `x`; `scaled` selects exponential scaling and range or convergence failures are errors.
fn k_evaluate(n: u32, x: f64, scaled: bool) -> Result<f64, Error> {
    validate(n, x, true)?;
    let (a, b) = k_base(x)?;
    if n > 0 && !b.is_finite() {
        return Err(Error::Overflow);
    }
    let (value, scale) = forward_log(n, x, a, b, true)?;
    finite_exp(value.ln() + scale - if scaled { 0.0 } else { x })
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Checks integral exhaustion.
    #[test]
    fn integral_exhaustion() {
        assert_eq!(k_integral(2.0, 1), Err(Error::NonConvergence));
        assert!(k_integral(2.0, 7).is_ok());
    }
    /// Checks miller exhaustion.
    #[test]
    fn miller_exhaustion() {
        assert_eq!(i_log_kernel(20, 15.0, 0), Err(Error::NonConvergence));
        assert!(i_log_kernel(20, 15.0, 8).is_ok());
        assert_eq!(j_value(20, 15.0, 0), Err(Error::NonConvergence));
        assert!(j_value(20, 15.0, 8).is_ok());
    }
}

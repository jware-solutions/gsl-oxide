use crate::{
    finite, ordered,
    rule::{geometry, weighted},
    target, validate_finite, Budget, Error, Sum,
};
use std::f64::consts::PI;

/// Integrates finite bounds with logarithmic or algebraic endpoint singularities.
///
/// Uses the original tanh–sinh transformation and successively halved trapezoidal
/// meshes: `x=mid+half*tanh(pi/2*sinh(u))`. Intended for functions smooth inside
/// the interval, with smooth endpoint factors multiplying logarithms or powers
/// `distance^alpha`, **alpha >= -1/2**. Stronger singularities, interior
/// discontinuities, and oscillatory endpoint behavior are outside this method's
/// supported class. Smooth nonsingular functions are also accepted.
///
/// Returns `(integral, estimated_absolute_error)` under the [common contracts](crate).
/// Endpoints are never evaluated. Separate absolute tail estimates assume the
/// documented endpoint behavior has set in: decreasing transformed terms and a
/// conservative `4*distance*abs(f(x))` envelope. Mesh differences, tail estimates,
/// mapping uncertainty, and summation roundoff all contribute to the error.
/// Agreement of successive meshes alone cannot certify convergence.
///
/// Only the rounded `x` is passed to the callback. Near a nonzero endpoint its
/// distance loses relative accuracy; inverse square roots there often limit
/// attainable absolute accuracy to roughly `1e-7` times the integral's scale.
/// At a zero endpoint much smaller distances can be represented. The routine
/// retains an estimate for unresolved tails and returns [`Error::Roundoff`] if
/// resolution or 16 mesh halvings prevent convergence. Smaller budgets can fail
/// earlier with [`Error::EvaluationLimit`]. There is no fixed epsilon exclusion.
///
/// Mathematics: [Takahasi and Mori, 1974](https://ems.press/journals/prims/articles/2686),
/// [Mori's review](https://ems.press/journals/prims/articles/2317).
/// See also [Boost's discussion of endpoint conditioning](https://www.boost.org/doc/libs/1_84_0/libs/math/doc/html/math_toolkit/double_exponential/de_tanh_sinh.html).
/// No implementation code or tables from these sources are used.
///
/// ```
/// let (value, error) = gsl_oxide_integration::integrate_singular(
///     |x| 1.0 / x.sqrt(), 0.0, 1.0, 1e-10, 1e-10, 20_000)?;
/// assert!((value - 2.0).abs() < 2e-10);
/// assert!(error <= 2e-10);
/// # Ok::<(), gsl_oxide_integration::Error>(())
/// ```
pub fn integrate_singular<F: FnMut(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    abs_tol: f64,
    rel_tol: f64,
    max_evals: usize,
) -> Result<(f64, f64), Error> {
    validate_finite(a, b, abs_tol, rel_tol)?;
    if a == b {
        return Ok((0.0, 0.0));
    }
    let (a, b, sign) = ordered(a, b);
    let (mid, half) = geometry(a, b)?;
    let width = b - a;
    let scale_uncertainty = if width.is_finite() {
        (2.0 * half - width).abs() / width
    } else {
        0.0
    };
    let mut budget = Budget::new(f, max_evals);
    let mut previous = None;
    let mut previous_error = f64::MAX;
    let mut stagnant = 0;
    let mut h = 0.5;
    // Meshes are recomputed, keeping only local sums. All repeated samples count.
    for level in 0..=16 {
        let mut sum = Sum::default();
        let mut absolute = Sum::default();
        let mut uncertainty = Sum::default();
        let central = weighted(budget.eval(mid)?, half, PI * 0.5 * h)?;
        sum.add(central)?;
        absolute.add(central.abs())?;
        let preliminary = previous.unwrap_or(central);
        let tail_goal = target(abs_tol, rel_tol, preliminary) / 16.0;
        let mut tails = 0.0;
        for right in [false, true] {
            let tail = side(
                &mut budget,
                a,
                b,
                half,
                h,
                right,
                tail_goal,
                &mut sum,
                &mut absolute,
                &mut uncertainty,
            )?;
            tails = finite(tails + tail)?;
        }
        let value = sum.total()?;
        let roundoff = ((64.0 * f64::EPSILON + scale_uncertainty) * absolute.total()?)
            .max(f64::from_bits(1) * budget.used as f64);
        let error = if let Some(old) = previous {
            finite(4.0 * (value - old).abs() + tails + uncertainty.total()? + roundoff)?
        } else {
            f64::MAX
        };
        if level >= 2 && error <= target(abs_tol, rel_tol, value) {
            return Ok((sign * value, error));
        }
        if level > 2 && error >= previous_error * 0.9 {
            stagnant += 1;
        } else {
            stagnant = 0;
        }
        if stagnant >= 3 || roundoff > target(abs_tol, rel_tol, value) {
            return Err(Error::Roundoff);
        }
        previous = Some(value);
        previous_error = error;
        h *= 0.5;
    }
    Err(Error::Roundoff)
}

#[allow(clippy::too_many_arguments)]
fn side<F: FnMut(f64) -> f64>(
    budget: &mut Budget<F>,
    a: f64,
    b: f64,
    half: f64,
    h: f64,
    right: bool,
    goal: f64,
    sum: &mut Sum,
    absolute: &mut Sum,
    uncertainty: &mut Sum,
) -> Result<f64, Error> {
    let endpoint = if right { b } else { a };
    let mut last_x = a + half;
    let mut history = [f64::INFINITY; 3];
    let mut envelope = f64::INFINITY;
    // Even on the finest supported mesh this reaches u=8, beyond underflow
    // of exp(-pi*sinh(u)). Thus representational collapse is encountered before
    // the loop bound, rather than imposing an unresolved finite tail cutoff.
    for k in 1..=1_048_576 {
        let u = f64::from(k) * h;
        // q/(1+q) retains endpoint distance even when tanh rounds to +/-1.
        let q = (-PI * u.sinh()).exp();
        let distance = half * (2.0 * q / (1.0 + q));
        let x = if right { b - distance } else { a + distance };
        if !(a < x && x < b) || x == last_x {
            // Never skip a collapsed node as though it contributed zero. Retain
            // the last tail envelope; subsequent meshes can improve it, or fail.
            return finite(envelope).map_err(|_| Error::Roundoff);
        }
        last_x = x;
        let y = budget.eval(x)?;
        let weight = 2.0 * PI * u.cosh() * q / (1.0 + q).powi(2) * h;
        let term = weighted(y, half, weight)?;
        sum.add(term)?;
        absolute.add(term.abs())?;
        let represented_distance = (x - endpoint).abs();
        // First-order endpoint sensitivity for the documented powers/logs;
        // retain rounding uncertainty instead of assuming exact transformed x.
        let relative_position_error =
            ((represented_distance - distance).abs() + f64::EPSILON * x.abs() + f64::from_bits(1))
                / distance;
        uncertainty.add(finite(term.abs() * relative_position_error.min(1.0))?)?;
        history.rotate_left(1);
        history[2] = term.abs();
        let decreasing = history[2] <= history[1] && history[1] <= history[0];
        if k >= 3 && u >= 2.0 && decreasing {
            // For alpha >= -1/2 the remaining integral is asymptotically at most
            // 2*d*|f|. Factor 4 leaves room for smooth factors/logarithms.
            let endpoint_tail = finite((y.abs() * represented_distance) * 4.0)?;
            let ratio = if history[0] > 0.0 {
                (history[2] / history[0]).sqrt()
            } else {
                0.0
            };
            let mesh_tail = if ratio < 1.0 {
                finite(2.0 * history[2] / (1.0 - ratio))?
            } else {
                f64::MAX
            };
            envelope = endpoint_tail.max(mesh_tail);
            if envelope <= goal {
                return Ok(envelope);
            }
        }
    }
    Err(Error::Roundoff)
}

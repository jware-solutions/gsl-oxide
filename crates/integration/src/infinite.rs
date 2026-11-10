use crate::{adaptive, finite, validate_tolerances, Budget, Error};

/// Integrates `(a, +infinity)` for a finite lower bound.
///
/// Returns `(integral, estimated_absolute_error)` under the [crate contracts](crate).
/// Uses `x=a+(1-t)/t`, Jacobian `1/t^2`, with open quadrature nodes in `(0,1)`;
/// the original closure is never called at infinity. Assumes a convergent
/// improper integral with sufficiently regular transformed tails. Very slow
/// decay or oscillations may exhaust the budget. No arbitrary finite cutoff is
/// introduced. Overflow or loss of distinct transformed arguments is an error.
/// See [GSL's infinite substitutions](https://www.gnu.org/software/gsl/doc/html/integration.html#qagi-adaptive-integration-on-infinite-intervals).
///
/// ```
/// let (value, _) = gsl_oxide_integration::integrate_upper_infinite(
///     |x| (-x).exp(), 0.0, 1e-11, 1e-11, 20_000)?;
/// assert!((value - 1.0).abs() < 1e-11);
/// # Ok::<(), gsl_oxide_integration::Error>(())
/// ```
pub fn integrate_upper_infinite<F: FnMut(f64) -> f64>(
    f: F,
    a: f64,
    abs_tol: f64,
    rel_tol: f64,
    max_evals: usize,
) -> Result<(f64, f64), Error> {
    semi(f, a, 1.0, abs_tol, rel_tol, max_evals)
}

/// Integrates `(-infinity, b)` for a finite upper bound.
///
/// Uses `x=b-(1-t)/t` and Jacobian `1/t^2`. Returns
/// `(integral, estimated_absolute_error)` with [`integrate_upper_infinite`]'s
/// convergence assumptions, restrictions, and [common contracts](crate).
pub fn integrate_lower_infinite<F: FnMut(f64) -> f64>(
    f: F,
    b: f64,
    abs_tol: f64,
    rel_tol: f64,
    max_evals: usize,
) -> Result<(f64, f64), Error> {
    semi(f, b, -1.0, abs_tol, rel_tol, max_evals)
}

/// Evaluates the budgeted integrand at the point mapped from `t` and `boundary` in direction `sign`, including the Jacobian. Invalid mappings and nonfinite results return errors.
fn transformed<F: FnMut(f64) -> f64>(
    budget: &mut Budget<F>,
    t: f64,
    boundary: f64,
    sign: f64,
) -> Result<f64, Error> {
    if !(0.0 < t && t < 1.0) {
        return Err(Error::Subdivision);
    }
    let offset = (1.0 - t) / t;
    let x = finite(sign.mul_add(offset, boundary))?;
    if (sign > 0.0 && x <= boundary) || (sign < 0.0 && x >= boundary) {
        return Err(Error::Subdivision);
    }
    // Sequential division avoids an underflowing t*t or an overflowing 1/t^2.
    finite((budget.eval(x)? / t) / t)
}

/// Integrates `f` from `boundary` toward the infinity selected by `sign`, using tolerances `abs`, `rel` and callback limit `limit`; returns the estimate and absolute error.
fn semi<F: FnMut(f64) -> f64>(
    f: F,
    boundary: f64,
    sign: f64,
    abs: f64,
    rel: f64,
    limit: usize,
) -> Result<(f64, f64), Error> {
    validate_tolerances(abs, rel)?;
    if !boundary.is_finite() {
        return Err(Error::InvalidInput);
    }
    let mut budget = Budget::new(f, limit);
    adaptive::run(
        &mut |t, _| transformed(&mut budget, t, boundary, sign),
        [(0.0, 1.0, 0)],
        abs,
        rel,
    )
}

/// Integrates the ordinary improper integral over `(-infinity, +infinity)`.
///
/// Both tails must converge independently. Two separately estimated transformed
/// intervals split at zero share one adaptive controller, total error criterion,
/// and original callback budget. Opposing tail errors are added, never canceled;
/// no symmetric pairing or principal-value interpretation is used. General
/// divergence detection for arbitrary callbacks is not possible. Detected
/// stagnation, overflow, and work exhaustion return errors.
/// Returns `(integral, estimated_absolute_error)`; see
/// [`integrate_upper_infinite`] and the [common contracts](crate).
pub fn integrate_infinite<F: FnMut(f64) -> f64>(
    f: F,
    abs_tol: f64,
    rel_tol: f64,
    max_evals: usize,
) -> Result<(f64, f64), Error> {
    validate_tolerances(abs_tol, rel_tol)?;
    let mut budget = Budget::new(f, max_evals);
    adaptive::run(
        &mut |t, part| transformed(&mut budget, t, 0.0, if part == 0 { -1.0 } else { 1.0 }),
        [(0.0, 1.0, 0), (0.0, 1.0, 1)],
        abs_tol,
        rel_tol,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Checks transform overflow collapse and jacobian range.
    #[test]
    fn transform_overflow_collapse_and_jacobian_range() {
        let mut budget = Budget::new(|_| panic!("invalid transformed argument"), 10);
        for t in [0.0, 1.0, -1.0, f64::NAN] {
            assert_eq!(
                transformed(&mut budget, t, 0.0, 1.0),
                Err(Error::Subdivision)
            );
        }
        assert_eq!(
            transformed(&mut budget, f64::MIN_POSITIVE, f64::MAX, 1.0),
            Err(Error::NumericalBreakdown)
        );
        assert_eq!(
            transformed(&mut budget, 0.5, f64::MAX, 1.0),
            Err(Error::Subdivision)
        );
        assert_eq!(budget.used, 0);
        let mut budget = Budget::new(|_| 1.0, 1);
        assert_eq!(
            transformed(&mut budget, 1e-200, 0.0, 1.0),
            Err(Error::NumericalBreakdown)
        );
        assert_eq!(budget.used, 1);
        // Do not form 1/t^2: it would overflow before multiplication by a tiny
        // but finite original callback value could recover a finite result.
        let mut budget = Budget::new(
            |x| {
                assert!(x.is_finite());
                1e-300
            },
            1,
        );
        let value = transformed(&mut budget, 1e-200, 0.0, 1.0).unwrap();
        assert!((value / 1e100 - 1.0).abs() < 5e-16);
    }
}

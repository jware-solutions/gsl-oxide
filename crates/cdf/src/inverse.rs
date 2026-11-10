use crate::Error;

// Bisection in log(x) handles quantiles spanning the entire f64 exponent
// range, including subnormal lower tails. No probability complement is formed.
/// Inverts the monotone probability callback `evaluate` at `p` over positive values up to `max`. Set `upper` for a decreasing survival function; returns an overflow or convergence error when necessary.
pub(crate) fn positive_quantile(
    p: f64,
    upper: bool,
    max: f64,
    evaluate: impl Fn(f64) -> Result<f64, Error>,
) -> Result<f64, Error> {
    let mut lo = -746.0;
    let mut hi = max.ln();
    let at_max = evaluate(max)?;
    if (!upper && at_max < p) || (upper && at_max > p) {
        return Err(Error::Overflow);
    }
    for _ in 0..128 {
        let mid = lo + (hi - lo) * 0.5;
        let x = mid.exp();
        let value = evaluate(x)?;
        if value == p {
            return Ok(x);
        }
        if (value < p) != upper {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo <= 2.0 * f64::EPSILON * (1.0 + mid.abs()) {
            return Ok((lo + (hi - lo) * 0.5).exp().min(max));
        }
    }
    Err(Error::NonConvergence)
}

// Finite comparisons only: NaN/infinity use direct assertions, and signed zero
// uses to_bits() in the calling tests. Approximate equality cannot check a zero's sign.
pub(crate) fn assert_close(actual: f64, expected: f64, abs_tol: f64, rel_tol: f64, context: &str) {
    assert!(
        abs_tol.is_finite() && abs_tol >= 0.0 && rel_tol.is_finite() && rel_tol >= 0.0,
        "{context}: invalid tolerances: abs_tol={abs_tol:?}, rel_tol={rel_tol:?}; \
         actual={actual:?}, expected={expected:?}"
    );
    assert!(
        actual.is_finite() && expected.is_finite(),
        "{context}: comparison requires finite values: actual={actual:?}, expected={expected:?}"
    );
    let absolute_error = (actual - expected).abs();
    let bound = abs_tol.max(rel_tol * expected.abs());
    assert!(
        absolute_error <= bound,
        "{context}: actual={actual:?}, expected={expected:?}, \
         absolute error={absolute_error:e}, allowed bound={bound:e}"
    );
}

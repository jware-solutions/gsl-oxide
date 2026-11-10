// Match the workspace's finite absolute/relative comparison convention.
/// Asserts that two values agree within the supplied absolute and relative tolerances.
pub(crate) fn close(actual: f64, expected: f64, abs: f64, rel: f64) {
    assert!(abs.is_finite() && abs >= 0.0 && rel.is_finite() && rel >= 0.0);
    assert!(actual.is_finite() && expected.is_finite());
    let bound = abs.max(rel * expected.abs());
    assert!(
        (actual - expected).abs() <= bound,
        "actual={actual:.17e}, expected={expected:.17e}, difference={:.4e}, bound={bound:.4e}",
        (actual - expected).abs()
    );
}

/// Asserts that the supplied value is the positive IEEE-754 zero.
pub(crate) fn positive_zero(value: f64) {
    assert_eq!(value.to_bits(), 0.0_f64.to_bits());
}

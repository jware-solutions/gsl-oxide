use gsl_oxide_integration::Error;

pub(crate) fn assert_close(actual: f64, expected: f64, abs: f64, rel: f64) {
    assert!(abs.is_finite() && abs >= 0.0 && rel.is_finite() && rel >= 0.0);
    assert!(actual.is_finite() && expected.is_finite());
    let bound = abs.max(rel * expected.abs());
    assert!(
        (actual - expected).abs() <= bound,
        "actual={actual:.17e}, expected={expected:.17e}, difference={:.4e}, bound={bound:.4e}",
        (actual - expected).abs()
    );
}

pub(crate) fn check(
    result: Result<(f64, f64), Error>,
    expected: f64,
    abs: f64,
    rel: f64,
) -> (f64, f64) {
    let (value, error) = result.unwrap();
    assert_close(value, expected, abs, rel);
    assert!(error.is_finite() && error >= 0.0);
    assert!(
        error <= abs.max(rel * value.abs()),
        "stopping target: {value} +/- {error}"
    );
    (value, error)
}

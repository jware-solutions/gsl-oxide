#![allow(dead_code)]
pub(crate) const TEST_TOL0: f64 = 2.0 * f64::EPSILON;
pub(crate) const TEST_TOL1: f64 = 16.0 * f64::EPSILON;
pub(crate) const TEST_TOL2: f64 = 256.0 * f64::EPSILON;
pub(crate) const TEST_TOL3: f64 = 2048.0 * f64::EPSILON;
pub(crate) const TEST_TOL4: f64 = 16384.0 * f64::EPSILON;
pub(crate) const TEST_TOL5: f64 = 131072.0 * f64::EPSILON;
pub(crate) const TEST_TOL6: f64 = 1048576.0 * f64::EPSILON;

// Same relative-error convention as gsl_test_rel: absolute tolerance at zero.
// Check infinity's sign as well as its presence.
/// Asserts that `actual` matches `expected` within the given relative tolerance, using `context` in failures.
pub(crate) fn assert_relative(actual: f64, expected: f64, tolerance: f64, context: &str) {
    let pass = if expected.is_infinite() {
        actual == expected
    } else if expected == 0.0 {
        actual.is_finite() && actual.abs() <= tolerance
    } else {
        actual.is_finite() && ((actual - expected) / expected).abs() <= tolerance
    };
    assert!(
        pass,
        "{context}: observed {actual:.18e}, expected {expected:.18e}, tolerance {tolerance:.3e}"
    );
}
#[allow(unused_macros)]
macro_rules! check {
    ($function:ident, ($($arg:expr),*), $expected:expr, $tolerance:expr) => {
        crate::support::assert_relative($function($($arg),*).expect(stringify!($function($($arg),*))),
            $expected, $tolerance, stringify!($function($($arg),*)))
    };
}
#[allow(unused_imports)]
pub(crate) use check;

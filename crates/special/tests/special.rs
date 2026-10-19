//! Tests for the implemented special-function public API.

use gsl_oxide_special::error_functions::erf;

fn assert_close(x: f64, expected: f64, abs_tol: f64, rel_tol: f64) {
    let actual = erf(x);
    let tolerance = abs_tol.max(rel_tol * expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "erf({x:?}) = {actual:?}, expected {expected:?}, tolerance {tolerance:e}"
    );
}

#[test]
fn reference_values() {
    // Rounded values of the defining Gaussian integral:
    // https://www.gnu.org/software/gsl/doc/html/specfunc.html#error-function
    // Independently checked with high-precision evaluation of its Taylor series
    // (DLMF 7.6.1), rather than the positive-term series used by the implementation.
    // Allow 16 epsilons relative error for recurrence, summation and exp rounding.
    let cases = [
        (0.0, 0.0),
        (1.0e-10, 1.128_379_167_095_512_6e-10),
        (0.1, 0.112_462_916_018_284_9),
        (0.5, 0.520_499_877_813_046_5),
        (0.75, 0.711_155_633_653_515_1),
        (1.0, 0.842_700_792_949_714_9),
        (-1.0, -0.842_700_792_949_714_9),
        (1.5, 0.966_105_146_475_310_8),
        (2.0, 0.995_322_265_018_952_7),
        (3.0, 0.999_977_909_503_001_4),
        (4.0, 0.999_999_984_582_742_1),
        (5.0, 0.999_999_999_998_462_6),
        (5.5, 0.999_999_999_999_992_7),
    ];
    for (x, expected) in cases {
        assert_close(x, expected, 0.0, 16.0 * f64::EPSILON);
    }
}

#[test]
fn non_finite_inputs() {
    assert!(erf(f64::NAN).is_nan());
    assert!(erf(-f64::NAN).is_nan());
    assert_eq!(erf(f64::INFINITY), 1.0);
    assert_eq!(erf(f64::NEG_INFINITY), -1.0);
}

#[test]
fn preserves_signed_zero() {
    assert_eq!(erf(0.0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(erf(-0.0).to_bits(), (-0.0_f64).to_bits());
}

#[test]
fn odd_symmetry() {
    for x in [1.0e-300, 1.0e-10, 0.1, 0.5, 1.0, 2.5, 5.5, 6.0, f64::MAX] {
        let expected = -erf(x);
        let actual = erf(-x);
        assert!((actual - expected).abs() <= 2.0 * f64::EPSILON * expected.abs());
    }
}

#[test]
fn tiny_inputs_follow_the_derivative_at_zero() {
    // erf(x) = 2*x/sqrt(pi) + O(x^3); the relative correction is at most x^2/3.
    for x in [f64::MIN_POSITIVE, 1.0e-300, 1.0e-160, 1.0e-10] {
        let expected = std::f64::consts::FRAC_2_SQRT_PI * x;
        assert_close(x, expected, 0.0, 4.0 * f64::EPSILON);
        assert_close(-x, -expected, 0.0, 4.0 * f64::EPSILON);
    }
    // Subnormal rounding gets an absolute allowance of one subnormal step.
    let step = f64::from_bits(1);
    for x in [step, f64::from_bits(16), f64::MIN_POSITIVE / 2.0] {
        let expected = std::f64::consts::FRAC_2_SQRT_PI * x;
        assert_close(x, expected, step, 0.0);
        assert_close(-x, -expected, step, 0.0);
        assert!(erf(x) > 0.0);
        assert!(erf(-x) < 0.0);
    }
}

#[test]
fn large_inputs_and_saturation_boundary() {
    // The tail bound exp(-x^2)/(x*sqrt(pi)) is < half an ulp below 1 at x >= 6.
    let below = f64::from_bits(6.0_f64.to_bits() - 1);
    let above = f64::from_bits(6.0_f64.to_bits() + 1);
    for x in [below, 6.0, above, 10.0, 1.0e200, f64::MAX] {
        assert_close(x, 1.0, 16.0 * f64::EPSILON, 0.0);
        assert_close(-x, -1.0, 16.0 * f64::EPSILON, 0.0);
        assert!((0.0..=1.0).contains(&erf(x)));
        assert!((-1.0..=0.0).contains(&erf(-x)));
    }
}

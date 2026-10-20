//! Tests for the implemented special-function public API.

use gsl_oxide_special::error_functions::erf;

// Finite comparisons only: NaN/infinity use direct assertions, and signed zero
// uses to_bits() below. Approximate equality does not establish a zero's sign.
fn assert_close(actual: f64, expected: f64, abs_tol: f64, rel_tol: f64, context: &str) {
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

#[test]
fn comparison_accepts_absolute_relative_and_exact_bounds() {
    // Binary fractions make the accepted boundaries exact in these checks.
    assert_close(1.25, 1.0, 0.25, 0.125, "absolute bound");
    assert_close(0.75, 1.0, 0.25, 0.125, "absolute bound, other side");
    assert_close(9.0, 8.0, 0.25, 0.125, "relative bound");
    assert_close(-9.0, -8.0, 0.25, 0.125, "negative relative bound");
    assert_close(1.0, 1.0, 0.0, 0.0, "zero tolerances");
    assert_close(0.0, 0.0, 0.0, 0.0, "exact zero");
}

#[test]
#[should_panic(expected = "allowed bound")]
fn comparison_rejects_just_outside_absolute_bound() {
    let actual = f64::from_bits(1.25_f64.to_bits() + 1);
    assert_close(actual, 1.0, 0.25, 0.125, "outside absolute bound");
}

#[test]
#[should_panic(expected = "allowed bound")]
fn comparison_rejects_just_outside_relative_bound() {
    let actual = f64::from_bits(9.0_f64.to_bits() + 1);
    assert_close(actual, 8.0, 0.25, 0.125, "outside relative bound");
}

#[test]
fn comparison_rejects_invalid_tolerances() {
    for invalid in [-f64::EPSILON, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for (abs_tol, rel_tol) in [(invalid, 0.0), (0.0, invalid)] {
            assert!(
                std::panic::catch_unwind(|| {
                    assert_close(1.0, 1.0, abs_tol, rel_tol, "invalid tolerance");
                })
                .is_err(),
                "accepted abs_tol={abs_tol:?}, rel_tol={rel_tol:?}"
            );
        }
    }
}

#[test]
fn comparison_rejects_non_finite_values() {
    for non_finite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for (actual, expected) in [
            (non_finite, 1.0),
            (1.0, non_finite),
            (non_finite, non_finite),
        ] {
            assert!(
                std::panic::catch_unwind(|| {
                    assert_close(actual, expected, 1.0, 1.0, "non-finite comparison");
                })
                .is_err(),
                "accepted actual={actual:?}, expected={expected:?}"
            );
        }
    }
}

#[test]
fn reference_values() {
    // Retain Phase 1's rounded literals, documented there as Taylor-series
    // evaluations (DLMF 7.6.1): https://dlmf.nist.gov/7.6.E1
    // DLMF supplies the series, not a numerical table. Published numerical
    // anchors at 0.5 and 2.0 also appear in the mpmath erf/erfinv examples:
    // https://mpmath.readthedocs.io/en/latest/functions/expintegrals.html#erf
    // The GSL manual supplies the definition, not these numerical literals:
    // https://www.gnu.org/software/gsl/doc/html/specfunc.html#error-function
    // Keep the existing 16-epsilon relative allowance: rounded f64 references
    // plus recurrence, summation and exp rounding; no GSL bitwise match claimed.
    const REL_TOL: f64 = 16.0 * f64::EPSILON;
    // Each row is (input, expected, absolute tolerance, relative tolerance).
    let cases = [
        (0.0, 0.0, 0.0, 0.0),
        (1.0e-10, 1.128_379_167_095_512_6e-10, 0.0, REL_TOL),
        (0.1, 0.112_462_916_018_284_9, 0.0, REL_TOL),
        (0.5, 0.520_499_877_813_046_5, 0.0, REL_TOL),
        (0.75, 0.711_155_633_653_515_1, 0.0, REL_TOL),
        (1.0, 0.842_700_792_949_714_9, 0.0, REL_TOL),
        (1.5, 0.966_105_146_475_310_8, 0.0, REL_TOL),
        (2.0, 0.995_322_265_018_952_7, 0.0, REL_TOL),
        (3.0, 0.999_977_909_503_001_4, 0.0, REL_TOL),
        (4.0, 0.999_999_984_582_742_1, 0.0, REL_TOL),
        (5.0, 0.999_999_999_998_462_6, 0.0, REL_TOL),
        (5.5, 0.999_999_999_999_992_7, 0.0, REL_TOL),
    ];
    for (x, expected, abs_tol, rel_tol) in cases {
        // Negative references use mathematical oddness (DLMF 7.4.1), not erf(x):
        // https://dlmf.nist.gov/7.4.E1
        for (input, reference) in [(x, expected), (-x, -expected)] {
            assert_close(
                erf(input),
                reference,
                abs_tol,
                rel_tol,
                &format!("reference erf({input:?})"),
            );
        }
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
    // A relationship check, supplementary to the independent reference values.
    // Allow two epsilons relative rounding, or one subnormal step at tiny scale.
    for x in [
        f64::from_bits(1),
        f64::from_bits(16),
        f64::MIN_POSITIVE,
        1.0e-300,
        1.0e-10,
        0.1,
        0.5,
        1.0,
        2.5,
        5.5,
        6.0,
        f64::MAX,
    ] {
        assert_close(
            erf(-x),
            -erf(x),
            f64::from_bits(1),
            2.0 * f64::EPSILON,
            &format!("odd symmetry at x={x:?}"),
        );
    }
}

#[test]
fn tiny_inputs_follow_the_derivative_at_zero() {
    // DLMF 7.6.1: erf(x) = 2*x/sqrt(pi) + O(x^3). Here |x| <= 1e-10,
    // so the relative correction <= x^2/3 is far below the rounding allowance.
    for x in [f64::MIN_POSITIVE, 1.0e-300, 1.0e-160, 1.0e-10] {
        let expected = std::f64::consts::FRAC_2_SQRT_PI * x;
        for (input, reference) in [(x, expected), (-x, -expected)] {
            assert_close(
                erf(input),
                reference,
                0.0,
                4.0 * f64::EPSILON,
                &format!("derivative near zero: erf({input:?})"),
            );
        }
    }
}

#[test]
fn subnormals_and_normal_transition() {
    // The cubic correction in DLMF 7.6.1 is smaller than a subnormal step here.
    // One step covers rounding the slope and product; relative bounds can
    // underflow. A strict sign check also prevents zero from passing at x=step.
    let step = f64::from_bits(1);
    let normal_bits = f64::MIN_POSITIVE.to_bits();
    for x in [
        step,
        f64::from_bits(2),
        f64::from_bits(16),
        f64::MIN_POSITIVE / 2.0,
        f64::from_bits(normal_bits - 1),
        f64::MIN_POSITIVE,
        f64::from_bits(normal_bits + 1),
    ] {
        let expected = std::f64::consts::FRAC_2_SQRT_PI * x;
        for (input, reference) in [(x, expected), (-x, -expected)] {
            let actual = erf(input);
            assert_close(
                actual,
                reference,
                step,
                0.0,
                &format!("subnormal/normal transition: erf({input:?})"),
            );
            assert!(
                actual != 0.0 && actual.is_sign_negative() == input.is_sign_negative(),
                "erf({input:?}) = {actual:?}: expected a nonzero result with the input's sign"
            );
        }
    }
}

#[test]
fn saturation_boundary() {
    // The Gaussian tail bound exp(-x^2)/(x*sqrt(pi)) is < half an ulp below 1
    // even at the predecessor of 6. Keep the existing 16-epsilon allowance for
    // the series evaluation below the branch; the documented branch is exact.
    let below = f64::from_bits(6.0_f64.to_bits() - 1);
    let above = f64::from_bits(6.0_f64.to_bits() + 1);
    // Negation reverses ordering: -above < -6 < -below.
    for (input, expected) in [(below, 1.0), (-below, -1.0)] {
        let actual = erf(input);
        assert_close(
            actual,
            expected,
            16.0 * f64::EPSILON,
            0.0,
            &format!("below saturation: erf({input:?})"),
        );
        assert!(actual.abs() <= 1.0, "erf({input:?}) = {actual:?}");
    }
    for x in [6.0, above] {
        assert_eq!(erf(x), 1.0, "erf({x:?}) must saturate");
        assert_eq!(erf(-x), -1.0, "erf({:?}) must saturate", -x);
    }
}

#[test]
fn very_large_finite_inputs() {
    // Exact unit limits also require finite, bounded, correctly signed results.
    // Squaring 1e200 or MAX would overflow; the supported result must not.
    for x in [10.0, 1.0e200, f64::MAX] {
        assert_eq!(erf(x), 1.0, "erf({x:?}) must reach the positive limit");
        assert_eq!(erf(-x), -1.0, "erf({:?}) must reach the negative limit", -x);
    }
}

#[test]
fn finite_outputs_have_mathematical_range_and_sign() {
    for magnitude in [
        0.0,
        f64::from_bits(1),
        f64::MIN_POSITIVE / 2.0,
        f64::MIN_POSITIVE,
        1.0e-300,
        1.0e-10,
        0.1,
        0.75,
        1.0,
        2.5,
        5.5,
        6.0,
        1.0e200,
        f64::MAX,
    ] {
        for x in [magnitude, -magnitude] {
            let actual = erf(x);
            assert!(
                actual.is_finite() && (-1.0..=1.0).contains(&actual),
                "erf({x:?}) = {actual:?}: expected a finite value in [-1, 1]"
            );
            if x == 0.0 {
                assert_eq!(actual, 0.0, "erf({x:?}) must be zero");
            } else {
                assert!(
                    actual != 0.0 && actual.is_sign_negative() == x.is_sign_negative(),
                    "erf({x:?}) = {actual:?}: expected a nonzero result with the input's sign"
                );
            }
        }
    }
}

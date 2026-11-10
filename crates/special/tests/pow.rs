//! GSL 2.8 power regression cases, IEEE edge cases, and downstream trait usage.
use gsl_oxide_special::{error_functions::erf, Pow as _};
use std::f64::consts::PI;

mod support;
use support::assert_close;

const SMALL_POWERS: [fn(&f64) -> f64; 8] = [
    f64::pow_2,
    f64::pow_3,
    f64::pow_4,
    f64::pow_5,
    f64::pow_6,
    f64::pow_7,
    f64::pow_8,
    f64::pow_9,
];

/// Checks downstream expression uses public trait.
#[test]
fn downstream_expression_uses_public_trait() {
    const N_FAC2_DBLE: [f64; 7] = [1.0, 1.0, 2.0, 3.0, 8.0, 15.0, 48.0];
    let ns = 3;
    let n32 = 2_u32;
    let x = 4.0_f64;
    let actual = N_FAC2_DBLE[2 * (ns - 1) + 2] / 2.0_f64.pow_uint(n32 + 1)
        * (PI / x.pow_uint(2 * n32 + 1)).sqrt();
    assert_close(
        actual,
        3.0 * PI.sqrt() / 16.0,
        0.0,
        2e-15,
        "downstream expression",
    );
    assert_eq!(erf(0.0), 0.0);
}

/// Checks small powers have exact binary anchors.
#[test]
fn small_powers_have_exact_binary_anchors() {
    let expected = [
        2.25,
        -3.375,
        5.0625,
        -7.59375,
        11.390625,
        -17.0859375,
        25.62890625,
        -38.443359375,
    ];
    for (power, expected) in SMALL_POWERS.into_iter().zip(expected) {
        assert_eq!(power(&-1.5), expected);
    }
}

/// Checks integer powers match gsl test cases.
#[test]
fn integer_powers_match_gsl_test_cases() {
    // GSL 2.8 sys/test.c compares these cases with libm's real-exponent pow.
    let base = -314.0_f64 / 100.0;
    for n in -9..=9 {
        assert_close(
            base.pow_int(n),
            base.powf(f64::from(n)),
            0.0,
            1e-15,
            "GSL signed powers",
        );
    }
    for n in 0..=9 {
        assert_close(
            base.pow_uint(n),
            base.powf(f64::from(n)),
            0.0,
            1e-15,
            "GSL unsigned powers",
        );
    }
    for (power, n) in SMALL_POWERS.into_iter().zip(2..=9) {
        assert_close(
            power(&base),
            base.powf(f64::from(n)),
            0.0,
            1e-15,
            "GSL small powers",
        );
    }
    for n in [i32::MIN, i32::MAX] {
        assert_close(
            1.0000001_f64.pow_int(n),
            1.0000001_f64.powf(f64::from(n)),
            0.0,
            1e-6,
            "GSL extreme exponents",
        );
    }
}

/// Checks full integer ranges preserve exponent and parity.
#[test]
fn full_integer_ranges_preserve_exponent_and_parity() {
    for n in [i32::MIN, i32::MIN + 1, -1, 0, 1, i32::MAX] {
        assert_eq!(1.0_f64.pow_int(n), 1.0);
        assert_eq!((-1.0_f64).pow_int(n), if n % 2 == 0 { 1.0 } else { -1.0 });
    }
    for n in [0, 1, i32::MAX as u32, i32::MAX as u32 + 1, u32::MAX] {
        assert_eq!(1.0_f64.pow_uint(n), 1.0);
        assert_eq!((-1.0_f64).pow_uint(n), if n % 2 == 0 { 1.0 } else { -1.0 });
    }
    assert_eq!(2.0_f64.pow_int(i32::MIN), 0.0);
    assert_eq!(0.5_f64.pow_int(i32::MIN), f64::INFINITY);
    assert_eq!(2.0_f64.pow_uint(u32::MAX), f64::INFINITY);
    assert_eq!(0.5_f64.pow_uint(u32::MAX), 0.0);
}

/// Checks exponent zero and one handle all bases.
#[test]
fn exponent_zero_and_one_handle_all_bases() {
    for base in [
        0.0,
        -0.0,
        1.0,
        -2.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ] {
        assert_eq!(base.pow_uint(0), 1.0);
        assert_eq!(base.pow_int(0), 1.0);
        if base.is_nan() {
            assert!(base.pow_uint(1).is_nan());
            assert!(base.pow_int(1).is_nan());
        } else {
            assert_eq!(base.pow_uint(1).to_bits(), base.to_bits());
            assert_eq!(base.pow_int(1).to_bits(), base.to_bits());
        }
    }
}

/// Checks zeros infinities and nan follow ieee arithmetic.
#[test]
fn zeros_infinities_and_nan_follow_ieee_arithmetic() {
    for n in [2, 3, 8, 9] {
        let negative_zero = if n % 2 == 0 { 0.0_f64 } else { -0.0_f64 };
        let negative_infinity = if n % 2 == 0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
        for (base, expected) in [
            (0.0_f64, 0.0_f64),
            (-0.0, negative_zero),
            (f64::INFINITY, f64::INFINITY),
            (f64::NEG_INFINITY, negative_infinity),
        ] {
            assert_eq!(base.pow_uint(n as u32).to_bits(), expected.to_bits());
            assert_eq!(base.pow_int(n).to_bits(), expected.to_bits());
            assert_eq!(base.pow_int(-n).to_bits(), (1.0 / expected).to_bits());
        }
        assert!(f64::NAN.pow_uint(n as u32).is_nan());
        assert!(f64::NAN.pow_int(n).is_nan());
        assert!(f64::NAN.pow_int(-n).is_nan());
    }
    for (power, n) in SMALL_POWERS.into_iter().zip(2..=9) {
        let negative_zero = if n % 2 == 0 { 0.0_f64 } else { -0.0_f64 };
        let negative_infinity = if n % 2 == 0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
        assert_eq!(power(&0.0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(power(&-0.0).to_bits(), negative_zero.to_bits());
        assert_eq!(power(&f64::INFINITY), f64::INFINITY);
        assert_eq!(power(&f64::NEG_INFINITY), negative_infinity);
        assert!(power(&f64::NAN).is_nan());
    }
}

/// Checks overflow underflow and reciprocal first keep gsl behavior.
#[test]
fn overflow_underflow_and_reciprocal_first_keep_gsl_behavior() {
    assert_eq!(2.0_f64.pow_uint(1023), 2.0_f64.powi(1023));
    assert_eq!(2.0_f64.pow_uint(1024), f64::INFINITY);
    assert_eq!(2.0_f64.pow_int(-1074).to_bits(), 1);
    assert_eq!(2.0_f64.pow_int(-1075).to_bits(), 0);
    assert_eq!((-2.0_f64).pow_int(-1075).to_bits(), (-0.0_f64).to_bits());
    // Squaring the original base first would overflow and lose this subnormal.
    assert_eq!(2.0_f64.powi(512).pow_int(-2).to_bits(), 1_u64 << 50);
    assert_eq!(f64::from_bits(1).pow_uint(1).to_bits(), 1);
    assert_eq!(f64::from_bits(1).pow_int(-1), f64::INFINITY);
    for (power, n) in SMALL_POWERS.into_iter().zip(2..=9) {
        assert_eq!(power(&f64::MAX), f64::INFINITY);
        assert_eq!(
            power(&-f64::MAX),
            if n % 2 == 0 {
                f64::INFINITY
            } else {
                f64::NEG_INFINITY
            }
        );
    }
}

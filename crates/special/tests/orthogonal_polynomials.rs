//! Polynomial anchors from explicit formulas, with identities as supplements.
use gsl_oxide_special::{laguerre::generalized, legendre, Error};
mod support;
use support::assert_close;

#[test]
fn legendre_low_orders_endpoints_and_central_values() {
    // DLMF 18.5.16, exact polynomial coefficients and special values.
    for x in [-1.0_f64, -0.75, -0.1, -0.0, 0.0, 0.25, 0.5, 0.9, 1.0] {
        assert_eq!(legendre::p(0, x), Ok(1.0));
        assert_eq!(legendre::p(1, x), Ok(x));
        assert_close(
            legendre::p(2, x).unwrap(),
            (3.0 * x * x - 1.0) / 2.0,
            3e-16,
            2e-15,
            "P2",
        );
        assert_close(
            legendre::p(3, x).unwrap(),
            x * (5.0 * x * x - 3.0) / 2.0,
            3e-16,
            2e-15,
            "P3",
        );
        assert_close(
            legendre::p(4, x).unwrap(),
            (35.0 * x.powi(4) - 30.0 * x * x + 3.0) / 8.0,
            5e-16,
            3e-15,
            "P4",
        );
    }
    for n in [0, 1, 2, 3, 99, 100, 999, 1000] {
        assert_eq!(legendre::p(n, 1.0), Ok(1.0));
        assert_eq!(
            legendre::p(n, -1.0),
            Ok(if n % 2 == 0 { 1.0 } else { -1.0 })
        );
        if n % 2 == 1 {
            assert_eq!(legendre::p(n, 0.0).unwrap().to_bits(), 0.0_f64.to_bits());
            assert_eq!(
                legendre::p(n, -0.0).unwrap().to_bits(),
                (-0.0_f64).to_bits()
            );
        }
    }
    // P_(2m)(0)=(-1)^m binomial(2m,m)/4^m (DLMF 18.5.7).
    for (n, expected) in [
        (2, -0.5),
        (4, 0.375),
        (10, -63.0 / 256.0),
        (100, 0.07958923738717877),
        (1000, 0.0252250181783608),
    ] {
        assert_close(
            legendre::p(n, 0.0).unwrap(),
            expected,
            0.0,
            5e-15,
            "central binomial anchors",
        );
    }
}

// Independent finite expression P_n(1-2z)=sum (-n)_k(n+1)_k z^k/(k!)^2,
// DLMF 18.5.7 with alpha=beta=0. This is not the degree recurrence.
// Used only where the absolute sum of terms is small, so cancellation
// in this test expression cannot swamp the chosen absolute tolerances.
fn legendre_explicit(n: u32, x: f64) -> f64 {
    let z = (1.0 - x) / 2.0;
    let mut term = 1.0;
    let mut sum = 1.0;
    for k in 1..=n {
        term *= -f64::from(n - k + 1) * f64::from(n + k) * z / f64::from(k).powi(2);
        sum += term;
    }
    sum
}

#[test]
fn legendre_high_degrees_and_difficult_regions() {
    for n in [100, 999, 1000] {
        for x in [
            1.0 - 2.0_f64.powi(-30),
            f64::from_bits(1.0_f64.to_bits() - 1),
            0.999999,
            0.99999711,
        ] {
            let expected = legendre_explicit(n, x);
            assert_close(
                legendre::p(n, x).unwrap(),
                expected,
                8e-15,
                8e-14,
                "high-degree explicit anchor",
            );
            let parity = if n % 2 == 0 { 1.0 } else { -1.0 };
            assert_close(
                legendre::p(n, -x).unwrap(),
                parity * expected,
                8e-15,
                8e-14,
                "negative endpoint neighbor",
            );
        }
    }
    // Exact P2/P3 roots, absolute allowance includes input rounding.
    for (n, root) in [(2, (1.0_f64 / 3.0).sqrt()), (3, (3.0_f64 / 5.0).sqrt())] {
        for bits in [root.to_bits() - 1, root.to_bits(), root.to_bits() + 1] {
            assert_close(
                legendre::p(n, f64::from_bits(bits)).unwrap(),
                0.0,
                8e-16,
                0.0,
                "root neighbors",
            );
        }
    }
    // Regression: the degree recurrence formerly rounded the least
    // subnormal input to one quantum instead of the expected 25 quanta.
    for x in [
        f64::from_bits(1),
        f64::from_bits(100),
        f64::from_bits(f64::MIN_POSITIVE.to_bits() - 1),
        f64::MIN_POSITIVE,
        f64::from_bits(f64::MIN_POSITIVE.to_bits() + 1),
        1e-100,
    ] {
        // P999'(0)=999*P998(0), higher terms negligible here.
        // P998(0)=-1000/999*P1000(0), from the binomial closed form.
        let expected = -1000.0 * 0.0252250181783608 * x;
        assert_close(
            legendre::p(999, x).unwrap(),
            expected,
            f64::from_bits(1),
            8e-15,
            "tiny odd polynomial",
        );
        assert!(legendre::p(999, x).unwrap() < 0.0);
        assert_eq!(legendre::p(999, -x).unwrap(), -legendre::p(999, x).unwrap());
        assert_close(
            legendre::p(1000, x).unwrap(),
            0.0252250181783608,
            0.0,
            8e-15,
            "tiny even polynomial",
        );
    }
}

#[test]
fn legendre_domains_and_supplemental_identities() {
    for n in [0, 1, 1000] {
        for x in [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            1.01,
            -1.01,
            f64::from_bits(1.0_f64.to_bits() + 1),
            -f64::from_bits(1.0_f64.to_bits() + 1),
        ] {
            assert_eq!(legendre::p(n, x), Err(Error::Domain));
        }
    }
    for n in [1001, u32::MAX] {
        assert_eq!(legendre::p(n, 0.0), Err(Error::UnsupportedRange));
        assert_eq!(legendre::p(n, 1.0), Err(Error::UnsupportedRange));
    }
    for x in [0.1_f64, 0.5, 0.9, 0.99] {
        for n in [1, 2, 10, 99, 500, 999] {
            let pn = legendre::p(n, x).unwrap();
            assert_eq!(
                legendre::p(n, -x).unwrap(),
                if n % 2 == 0 { pn } else { -pn }
            );
            let k = f64::from(n);
            let residual = (k + 1.0) * legendre::p(n + 1, x).unwrap() - (2.0 * k + 1.0) * x * pn
                + k * legendre::p(n - 1, x).unwrap();
            assert_close(residual, 0.0, 6e-13, 0.0, "degree recurrence");
        }
    }
    for bits in [
        0.5_f64.to_bits() - 1,
        0.5_f64.to_bits(),
        0.5_f64.to_bits() + 1,
    ] {
        // Independent explicit expression is well conditioned at degree 10.
        let x = f64::from_bits(bits);
        assert_close(
            legendre::p(10, x).unwrap(),
            legendre_explicit(10, x),
            6e-14,
            0.0,
            "recurrence branch neighbors",
        );
    }
}

#[test]
fn laguerre_independent_values_and_low_degree_formulas() {
    // Published examples, no external oracle execution:
    // https://mpmath.readthedocs.io/en/latest/functions/orthogonal.html#laguerre-polynomials
    assert_close(
        generalized(5, 0.0, 0.25).unwrap(),
        0.037_263_997_395_833_33,
        2e-16,
        2e-15,
        "ordinary L5",
    );
    assert_eq!(generalized(2, 0.0, 10000.0), Ok(49980001.0));
    // DLMF 18.5.17_5: exact low-order formulas, including noninteger alpha.
    for alpha in [-0.75_f64, -0.5, 0.0, 0.5, 2.5, 1000.0] {
        for x in [-10000.0_f64, -10.0, -0.0, 0.0, 0.25, 1.0, 10.0, 10000.0] {
            assert_eq!(generalized(0, alpha, x), Ok(1.0));
            assert_eq!(generalized(1, alpha, x), Ok(alpha + 1.0 - x));
            let l2 = 0.5 * x * x - (alpha + 2.0) * x + 0.5 * (alpha + 1.0) * (alpha + 2.0);
            // Integer coefficients of 6*L3 avoid rounding each 1/6 term
            // before cancellation (all these inputs are dyadic rationals).
            let l3 = (((-x + 3.0 * (alpha + 3.0)) * x - 3.0 * (alpha + 2.0) * (alpha + 3.0)) * x
                + (alpha + 1.0) * (alpha + 2.0) * (alpha + 3.0))
                / 6.0;
            assert_close(
                generalized(2, alpha, x).unwrap(),
                l2,
                4e-15,
                3e-15,
                "generalized L2",
            );
            assert_close(
                generalized(3, alpha, x).unwrap(),
                l3,
                4e-15,
                3e-15,
                "generalized L3",
            );
        }
    }
    // At x=0, L_n^alpha(0)=binomial(n+alpha,n), DLMF 18.5.12.
    for n in [2, 10, 100, 999, 1000] {
        assert_eq!(generalized(n, 0.0, 0.0), Ok(1.0));
        assert_close(
            generalized(n, 1.0, -0.0).unwrap(),
            f64::from(n + 1),
            0.0,
            1e-14, // up to 1000 rounded factors in the zero-argument product
            "binomial alpha=1",
        );
        assert_close(
            generalized(n, 2.0, 0.0).unwrap(),
            f64::from(n + 1) * f64::from(n + 2) / 2.0,
            0.0,
            1e-14,
            "binomial alpha=2",
        );
    }
}

// DLMF 18.5.12, the explicit finite coefficient sum. The coefficients
// are constructed independently via products, not the production degree
// recurrence. Positive x is restricted to |n*x|<=2, where cancellation
// is mild; negative x has entirely positive terms. These are independently
// justified exact finite expressions (with ordinary f64 rounding).
fn laguerre_explicit(n: u32, alpha: f64, x: f64) -> f64 {
    let mut sum = 0.0;
    let mut power_over_factorial = 1.0;
    for j in 0..=n {
        let mut binomial = 1.0;
        for k in 1..=n - j {
            binomial *= (alpha + f64::from(j + k)) / f64::from(k);
        }
        sum += binomial * power_over_factorial;
        power_over_factorial *= -x / f64::from(j + 1);
    }
    sum
}

#[test]
fn laguerre_high_degrees_near_roots_and_alpha_boundary() {
    let near_minus_one = f64::from_bits((-1.0_f64).to_bits() - 1);
    for (n, alpha, x, absolute, relative) in [
        // Up to 1000 propagated recurrence steps: allow 3e-13 relative;
        // near roots use absolute error against the coefficient formula.
        (1000, 0.0, 0.001, 8e-15, 3e-13),
        // Near the first root, 1000 propagated steps lose relative accuracy;
        // allow 2e-14 absolute against the well-conditioned coefficient sum.
        (1000, 0.0, 0.001445, 2e-14, 3e-13),
        (1000, 0.5, 0.001, 0.0, 3e-13),
        (1000, -0.75, 0.001, 8e-15, 3e-13),
        (1000, near_minus_one, 0.001, 2e-15, 3e-13),
        (1000, near_minus_one, -0.001, 2e-15, 3e-13),
        (1000, near_minus_one, 0.0, 0.0, 3e-14),
        (1000, 2.5, -1.0, 0.0, 3e-13),
        (250, 0.0, -1000.0, 0.0, 3e-13), // power-of-two scaling, finite result
    ] {
        let expected = laguerre_explicit(n, alpha, x);
        assert_close(
            generalized(n, alpha, x).unwrap(),
            expected,
            absolute,
            relative,
            "high-degree finite expression",
        );
    }
    // Positive coefficients for negative x make this independent finite
    // sum safe even at the overflow transition. Keep finite final values
    // when multiplying unscaled recurrence intermediates would overflow.
    for n in [260, 280, 285, 290, 300] {
        let expected = laguerre_explicit(n, 0.0, -1000.0);
        if expected.is_finite() {
            assert_close(
                generalized(n, 0.0, -1000.0).unwrap(),
                expected,
                0.0,
                3e-13,
                "Laguerre overflow boundary",
            );
        } else {
            assert_eq!(generalized(n, 0.0, -1000.0), Err(Error::Overflow));
        }
    }
    // L2^alpha roots = alpha+2 +/- sqrt(alpha+2), DLMF 18.5.17_5.
    for alpha in [0.0_f64, 0.5, -0.75] {
        for root in [
            alpha + 2.0 - (alpha + 2.0).sqrt(),
            alpha + 2.0 + (alpha + 2.0).sqrt(),
        ] {
            for bits in [root.to_bits() - 1, root.to_bits(), root.to_bits() + 1] {
                assert_close(
                    generalized(2, alpha, f64::from_bits(bits)).unwrap(),
                    0.0,
                    4e-15,
                    0.0,
                    "Laguerre root neighbors",
                );
            }
        }
    }
    for x in [
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MIN_POSITIVE,
        -f64::MIN_POSITIVE,
    ] {
        assert_eq!(generalized(1000, 0.0, x), Ok(1.0));
        assert_close(
            generalized(1000, near_minus_one, x).unwrap(),
            (near_minus_one + 1.0) / 1000.0,
            0.0,
            3e-12,
            "near -1 with tiny x",
        );
    }
}

#[test]
fn laguerre_domains_limits_overflow_and_recurrence() {
    for n in [0, 1, 1000] {
        for alpha in [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            -1.0,
            -2.0,
            f64::from_bits((-1.0_f64).to_bits() + 1),
        ] {
            assert_eq!(generalized(n, alpha, 0.0), Err(Error::Domain));
        }
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(generalized(n, 0.0, x), Err(Error::Domain));
        }
        assert_eq!(
            generalized(n, f64::from_bits(1000.0_f64.to_bits() + 1), 0.0),
            Err(Error::UnsupportedRange)
        );
        for x in [
            f64::from_bits(10000.0_f64.to_bits() + 1),
            -f64::from_bits(10000.0_f64.to_bits() + 1),
        ] {
            assert_eq!(generalized(n, 0.0, x), Err(Error::UnsupportedRange));
        }
    }
    for n in [1001, u32::MAX] {
        assert_eq!(generalized(n, 0.0, 0.0), Err(Error::UnsupportedRange));
    }
    assert_eq!(generalized(1000, 1000.0, 0.0), Err(Error::Overflow));
    assert_eq!(generalized(1000, 0.0, -10000.0), Err(Error::Overflow));
    assert_eq!(generalized(1000, 0.0, 10000.0), Err(Error::Overflow));
    for magnitude in [1.0_f64, 10000.0] {
        for sign in [-1.0, 1.0] {
            // Adjacent inputs at the recurrence switch and finite range
            // boundary. Degree 2 has an independently exact formula.
            for bits in [
                magnitude.to_bits() - 1,
                magnitude.to_bits(),
                magnitude.to_bits() + 1,
            ] {
                let x = sign * f64::from_bits(bits);
                if x.abs() > 10000.0 {
                    assert_eq!(generalized(2, 0.0, x), Err(Error::UnsupportedRange));
                    continue;
                }
                let expected = (x * x - 4.0 * x + 2.0) / 2.0;
                assert_close(
                    generalized(2, 0.0, x).unwrap(),
                    expected,
                    5e-16,
                    3e-15,
                    "Laguerre branch neighbors",
                );
            }
        }
    }
    for alpha in [-0.75, 0.0, 0.5, 10.0] {
        for x in [-1.0, 0.01, 0.25, 1.0, 100.0] {
            for n in [1, 2, 10, 99, 999] {
                let k = f64::from(n);
                let a = (k + 1.0) * generalized(n + 1, alpha, x).unwrap();
                let b = (2.0 * k + alpha + 1.0 - x) * generalized(n, alpha, x).unwrap()
                    - (k + alpha) * generalized(n - 1, alpha, x).unwrap();
                assert_close(a, b, 2e-12, 3e-13, "Laguerre recurrence");
            }
        }
    }
}

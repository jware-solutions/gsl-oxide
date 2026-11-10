//! Bessel anchors, stable order recurrences, scaled forms and numerical limits.
use gsl_oxide_special::{bessel::*, Error};
mod support;
use support::assert_close;

/// Checks published anchors.
#[test]
fn published_anchors() {
    // Numerical examples, no oracle executed:
    // https://mpmath.readthedocs.io/en/latest/functions/bessel.html
    // Decimal anchors have 15+ digits; tolerate f64 series/ODE/recurrence
    // rounding with relative 8e-14 (absolute 2e-15 near oscillatory zeros).
    for (n, x, y) in [
        (2, 1000.0, -0.024777229528606),
        (4, 0.75, 0.000801070086542314),
        (1, std::f64::consts::PI, 0.28461534317975276),
        (0, 10000.0, -0.0070961603533888015),
    ] {
        assert_close(jn(n, x).unwrap(), y, 2e-15, 8e-14, "published J");
    }
    assert_close(
        y1(std::f64::consts::PI).unwrap(),
        0.35887291677671896,
        2e-15,
        8e-14,
        "published Y1",
    );
    assert_close(
        y0(10000.0).unwrap(),
        0.003647805558986606,
        2e-15,
        8e-14,
        "published Y0",
    );
    assert_close(
        i0(1.0).unwrap(),
        1.2660658777520084,
        0.0,
        8e-14,
        "published I0",
    );
    assert_close(
        i_n(3, 2.3).unwrap(),
        0.349223221159309,
        0.0,
        8e-14,
        "published I3",
    );
    assert_close(
        k0(1.0).unwrap(),
        0.42102443824070834,
        0.0,
        8e-14,
        "published K0",
    );
    assert_close(
        k0(100.0).unwrap(),
        4.656628229175902e-45,
        0.0,
        8e-14,
        "published large K0",
    );
}

/// Checks published intermediate and scaled anchors.
#[test]
fn published_intermediate_and_scaled_anchors() {
    // Published SciPy documentation examples, not library execution:
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ive.html
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.kve.html
    // Full scalar literals use 8e-14 relative; rounded array output uses
    // 5.1e-9 absolute, half the last printed decimal plus roundoff.
    assert_close(
        i_n_scaled(3, 1000.0).unwrap(),
        0.01256056218254712,
        0.0,
        8e-14,
        "published large scaled I",
    );
    assert_close(
        kn_scaled(3, 1000.0).unwrap(),
        0.03980696128440973,
        0.0,
        8e-14,
        "published large scaled K",
    );
    for (n, x, v) in [
        (0, 3.0, 0.24300035),
        (1, 3.0, 0.19682671),
        (2, 3.0, 0.11178255),
        (0, 2.0, 0.30850832),
        (1, 2.0, 0.21526929),
        (2, 2.0, 0.09323903),
    ] {
        assert_close(
            i_n_scaled(n, x).unwrap(),
            v,
            5.1e-9,
            0.0,
            "published intermediate I",
        );
    }
    for (n, x, v) in [
        (0, 3.0, 0.6977616),
        (1, 3.0, 0.80656348),
        (2, 3.0, 1.23547058),
        (0, 10.0, 0.39163193),
        (1, 10.0, 0.41076657),
        (2, 10.0, 0.47378525),
    ] {
        assert_close(
            kn_scaled(n, x).unwrap(),
            v,
            5.1e-9,
            0.0,
            "published intermediate K",
        );
    }
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yv.html
    for (n, x, v) in [
        (0, 3.0, 0.37685001),
        (1, 3.0, 0.32467442),
        (0, 8.0, 0.22352149),
        (1, 8.0, -0.15806046),
    ] {
        assert_close(
            yn(n, x).unwrap(),
            v,
            5.1e-9,
            0.0,
            "published intermediate Y",
        );
    }
}

/// Checks larger order anchors and transition.
#[test]
fn larger_order_anchors_and_transition() {
    // Independent rational partial sums at x=1 (DLMF 10.2.2/10.25.2).
    // Five correction terms leave <3e-15 relative; no implementation call
    // supplies a reference value. factorial(20)=2432902008176640000 exactly.
    let lead = 2.0_f64.powi(-20) / 2_432_902_008_176_640_000.0;
    let c1 = 1.0 / (4.0 * 21.0);
    let c2 = 1.0 / (32.0 * 21.0 * 22.0);
    let c3 = 1.0 / (384.0 * 21.0 * 22.0 * 23.0);
    let c4 = 1.0 / (6144.0 * 21.0 * 22.0 * 23.0 * 24.0);
    let c5 = c4 / (4.0 * 5.0 * 25.0);
    assert_close(
        jn(20, 1.0).unwrap(),
        lead * (1.0 - c1 + c2 - c3 + c4 - c5),
        0.0,
        1e-13,
        "J20 rational anchor",
    );
    assert_close(
        i_n(20, 1.0).unwrap(),
        lead * (1.0 + c1 + c2 + c3 + c4 + c5),
        0.0,
        1e-13,
        "I20 rational anchor",
    );
    // Dominant finite sums, DLMF 10.8.1/10.31.1. The logarithmic J/I
    // terms are <1e-45 relative; the next finite term is <5e-17 relative.
    let lead = 2.0_f64.powi(20) * 121_645_100_408_832_000.0;
    let c1 = 1.0 / (4.0 * 19.0);
    let c2 = 1.0 / (32.0 * 19.0 * 18.0);
    let c3 = 1.0 / (384.0 * 19.0 * 18.0 * 17.0);
    let c4 = 1.0 / (6144.0 * 19.0 * 18.0 * 17.0 * 16.0);
    let c5 = c4 / (4.0 * 5.0 * 15.0);
    let c6 = c5 / (4.0 * 6.0 * 14.0);
    assert_close(
        yn(20, 1.0).unwrap(),
        -lead / std::f64::consts::PI * (1.0 + c1 + c2 + c3 + c4 + c5 + c6),
        0.0,
        1e-13,
        "Y20 finite sum anchor",
    );
    assert_close(
        kn(20, 1.0).unwrap(),
        lead / 2.0 * (1.0 - c1 + c2 - c3 + c4 - c5 + c6),
        0.0,
        1e-13,
        "K20 finite sum anchor",
    );
    // DLMF 10.19.8 leading Airy transition term, Ai(0) from 9.2.3.
    // Its O(n^-4/3) relative correction permits a 1e-4 allowance at n=1000.
    let leading = 2.0_f64.cbrt() * 0.3550280538878172 / 10.0;
    assert_close(
        jn(1000, 1000.0).unwrap(),
        leading,
        0.0,
        1e-4,
        "large J transition",
    );
    assert_close(
        yn(1000, 1000.0).unwrap(),
        -3.0_f64.sqrt() * leading,
        0.0,
        1e-4,
        "large Y transition",
    );
    for n in [20, 1000] {
        for x in [f64::from(n) - 1.0, f64::from(n), f64::from(n) + 1.0] {
            let wronskian = jn(n, x).unwrap() * yn(n - 1, x).unwrap()
                - jn(n - 1, x).unwrap() * yn(n, x).unwrap();
            assert_close(
                wronskian,
                2.0 / (std::f64::consts::PI * x),
                0.0,
                2e-12,
                "order transition Wronskian",
            );
            assert!(i_n_scaled(n, x).unwrap() > 0.0);
            assert!(kn_scaled(n, x).unwrap() > 0.0);
        }
    }
}

/// Checks small values and zeros.
#[test]
fn small_values_and_zeros() {
    // Published numerical examples at x=1 (rounded to f64):
    // https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.special.jv.html
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.yv.html
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.i1.html
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.k1.html
    // J2 is published in mpmath's besseljzero examples (URL above).
    // Order-two Y/I/K references are arithmetic on these published base
    // values, using DLMF 10.6.1/10.29.1, not calls to this implementation.
    let y0_ref = 0.088256964215677;
    let y1_ref = -0.7812128213002889;
    let i0_ref = 1.2660658777520084;
    let i1_ref = 0.5651591039924851;
    let k0_ref = 0.42102443824070834;
    let k1_ref = 0.6019072301972346;
    for (v, y) in [
        (j0(1.0), 0.7651976865579666),
        (j1(1.0), 0.44005058574493355),
        (jn(2, 1.0), 0.11490348493190048),
        (y0(1.0), y0_ref),
        (y1(1.0), y1_ref),
        (yn(2, 1.0), 2.0 * y1_ref - y0_ref),
        (i1(1.0), i1_ref),
        (i_n(2, 1.0), i0_ref - 2.0 * i1_ref),
        (k1(1.0), k1_ref),
        (kn(2, 1.0), k0_ref + 2.0 * k1_ref),
    ] {
        assert_close(v.unwrap(), y, 0.0, 6e-14, "small Bessel anchor");
    }
    // Published zero examples, same mpmath URL, #besseljzero/#besselyzero.
    // Rounded roots have small argument error; use absolute bounds.
    assert_close(j0(2.404825557695773).unwrap(), 0.0, 3e-15, 0.0, "J0 zero");
    assert_close(j1(3.8317059702075125).unwrap(), 0.0, 3e-15, 0.0, "J1 zero");
    assert_close(y0(0.8935769662791675).unwrap(), 0.0, 3e-15, 0.0, "Y0 zero");
    for x in [f64::from_bits(2), 1e-300, 1e-20] {
        assert_close(j1(x).unwrap(), x / 2.0, 0.0, 2e-15, "tiny J1");
        assert_close(i1(x).unwrap(), x / 2.0, 0.0, 2e-15, "tiny I1");
        assert_close(
            y0(x).unwrap(),
            2.0 / std::f64::consts::PI * (x.ln() - std::f64::consts::LN_2 + 0.5772156649015329),
            0.0,
            4e-15,
            "tiny Y0",
        );
        assert_close(
            k0(x).unwrap(),
            -(x.ln() - std::f64::consts::LN_2 + 0.5772156649015329),
            0.0,
            4e-15,
            "tiny K0",
        );
    }
}

/// Checks recurrence parity and wronskians.
#[test]
fn recurrence_parity_and_wronskians() {
    for x in [0.1_f64, 1.0, 2.0, 7.5, 20.0, 31.0, 32.0, 100.0, 1000.0] {
        for n in [1, 2, 5, 20, 50] {
            let k = 2.0 * f64::from(n) / x;
            assert_close(
                jn(n - 1, x).unwrap() + jn(n + 1, x).unwrap(),
                k * jn(n, x).unwrap(),
                2e-14,
                1e-12,
                "J recurrence",
            );
            assert_close(
                yn(n - 1, x).unwrap() + yn(n + 1, x).unwrap(),
                k * yn(n, x).unwrap(),
                2e-14,
                1e-12,
                "Y recurrence",
            );
            assert_close(
                i_n_scaled(n - 1, x).unwrap() - i_n_scaled(n + 1, x).unwrap(),
                k * i_n_scaled(n, x).unwrap(),
                2e-14,
                1e-12,
                "I recurrence",
            );
            assert_close(
                kn_scaled(n + 1, x).unwrap() - kn_scaled(n - 1, x).unwrap(),
                k * kn_scaled(n, x).unwrap(),
                2e-14,
                1e-12,
                "K recurrence",
            );
            assert_close(
                jn(n, -x).unwrap(),
                if n % 2 == 0 {
                    jn(n, x).unwrap()
                } else {
                    -jn(n, x).unwrap()
                },
                0.0,
                0.0,
                "J parity",
            );
            assert_close(
                i_n_scaled(n, -x).unwrap(),
                if n % 2 == 0 {
                    i_n_scaled(n, x).unwrap()
                } else {
                    -i_n_scaled(n, x).unwrap()
                },
                0.0,
                0.0,
                "I parity",
            );
        }
        // DLMF 10.5.2 and 10.28.2, with scaled I/K exponentials cancelling.
        assert_close(
            j1(x).unwrap() * y0(x).unwrap() - j0(x).unwrap() * y1(x).unwrap(),
            2.0 / (std::f64::consts::PI * x),
            0.0,
            8e-14,
            "JY Wronskian",
        );
        assert_close(
            i0_scaled(x).unwrap() * k1_scaled(x).unwrap()
                + i1_scaled(x).unwrap() * k0_scaled(x).unwrap(),
            1.0 / x,
            0.0,
            8e-14,
            "IK Wronskian",
        );
    }
}

/// Checks scaled range and recovery.
#[test]
fn scaled_range_and_recovery() {
    for x in [0.1_f64, 1.0, 10.0, 32.0, 100.0, 700.0] {
        for n in [0, 1, 2, 20] {
            assert_close(
                i_n_scaled(n, x).unwrap(),
                (-x).exp() * i_n(n, x).unwrap(),
                0.0,
                2e-13,
                "I scaling",
            );
            assert_close(
                kn_scaled(n, x).unwrap(),
                x.exp() * kn(n, x).unwrap(),
                0.0,
                2e-13,
                "K scaling",
            );
        }
    }
    for n in [0, 1, 2, 20, 1000] {
        assert_eq!(i_n(n, 1e6), Err(Error::Overflow));
        assert_eq!(kn(n, 1e6), Ok(0.0));
        assert!(i_n_scaled(n, 1e6).unwrap() > 0.0);
        assert!(kn_scaled(n, 1e6).unwrap() > 0.0);
    }
    // Regression: a scaled I underflow must not discard a finite I result.
    assert!(i_n(1000, 400.0).unwrap() > 0.0);
    assert_eq!(i_n_scaled(1000, 400.0), Ok(0.0));
    // At x=1e6, fixed small-order leading approximations have <4e-7 relative
    // remainder (DLMF 10.40). A separate broad asymptotic check.
    assert_close(
        i0_scaled(1e6).unwrap(),
        1.0 / (2.0 * std::f64::consts::PI * 1e6).sqrt(),
        0.0,
        2e-7,
        "large I scale",
    );
    assert_close(
        k1_scaled(1e6).unwrap(),
        (std::f64::consts::PI / (2.0 * 1e6)).sqrt(),
        0.0,
        4e-7,
        "large K scale",
    );
}

/// Checks boundaries and checked errors.
#[test]
fn boundaries_and_checked_errors() {
    for boundary in [1.0_f64, 32.0] {
        for x in [
            f64::from_bits(boundary.to_bits() - 1),
            boundary,
            f64::from_bits(boundary.to_bits() + 1),
        ] {
            for (actual, expected) in [
                (j0(x), j0(boundary)),
                (j1(x), j1(boundary)),
                (y0(x), y0(boundary)),
                (y1(x), y1(boundary)),
                (i0_scaled(x), i0_scaled(boundary)),
                (i1_scaled(x), i1_scaled(boundary)),
                (k0_scaled(x), k0_scaled(boundary)),
                (k1_scaled(x), k1_scaled(boundary)),
            ] {
                assert_close(
                    actual.unwrap(),
                    expected.unwrap(),
                    3e-14,
                    1e-13,
                    "base branch",
                );
            }
        }
    }
    for n in [0, 1, 2, 20, 1000] {
        assert!(jn(n, 1e6).unwrap().is_finite());
        assert!(yn(n, 1e6).unwrap().is_finite());
        for f in [jn, yn, i_n, i_n_scaled, kn, kn_scaled] {
            assert_eq!(f(n, f64::NAN), Err(Error::Domain));
            assert_eq!(f(n, f64::INFINITY), Err(Error::UnsupportedRange));
            assert_eq!(
                f(n, f64::from_bits(1e6_f64.to_bits() + 1)),
                Err(Error::UnsupportedRange)
            );
            assert_eq!(f(1001, 1.0), Err(Error::UnsupportedRange));
            assert_eq!(f(u32::MAX, 1.0), Err(Error::UnsupportedRange));
        }
        for f in [yn, kn, kn_scaled] {
            for x in [0.0, -0.0, -1.0, f64::NEG_INFINITY] {
                assert_eq!(f(n, x), Err(Error::Domain));
            }
        }
        for f in [jn, i_n, i_n_scaled] {
            assert_eq!(f(n, f64::NEG_INFINITY), Err(Error::UnsupportedRange));
            let y = f(n, -0.0).unwrap();
            assert_eq!(y, if n == 0 { 1.0 } else { 0.0 });
            if n % 2 == 1 {
                assert!(y.is_sign_negative());
            }
        }
    }
    assert_eq!(yn(1, f64::from_bits(1)), Err(Error::Overflow));
    assert_eq!(kn(1, f64::from_bits(1)), Err(Error::Overflow));
    assert_eq!(yn(1000, 1.0), Err(Error::Overflow));
    assert_eq!(kn(1000, 1.0), Err(Error::Overflow));
    assert_eq!(jn(1000, 1.0), Ok(0.0));
    assert_eq!(i_n(1000, 1.0), Ok(0.0));
}

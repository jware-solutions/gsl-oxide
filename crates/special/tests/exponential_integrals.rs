//! Independent anchors, real continuations and range regressions.
use gsl_oxide_special::{exponential_integrals::*, Error};
mod support;
use support::assert_close;

#[test]
fn independent_values() {
    // mpmath's published examples (read only, no oracle executed):
    // https://mpmath.readthedocs.io/en/latest/functions/expintegrals.html
    // E1(4), E1(2.5), E1(6.25), Re E1(-2.5), Ei(pi).
    for (x, expected) in [
        (4.0, 0.0037793524098489065),
        (2.5, 0.024914917870269736),
        (6.25, 0.0002704758872637179),
        (-2.5, -7.073765894578601),
    ] {
        assert_close(e1(x).unwrap(), expected, 0.0, 6e-15, "E1 examples");
    }
    assert_close(
        ei(std::f64::consts::PI).unwrap(),
        10.92837438933141,
        0.0,
        3e-15,
        "Ei(pi)",
    );

    // A&S ch.5 / DLMF 6.6.1-2, 8.19.8 (decimal evaluations).
    // Additional spot values at 50.5 are published in Boost.Math's
    // test_expint.hpp, attributed there to functions.wolfram.com.
    // Only numerical facts are used; this test code is independent.
    for (x, expected, relative) in [
        (0.125, -1.3732085249429833, 3e-15),
        (0.5, 0.4542199048631736, 3e-15),
        (1.0, 1.8951178163559368, 3e-15),
        (-0.125, -1.6234256405841688, 3e-15),
        (-0.5, -0.5597735947761608, 3e-15),
        (-1.0, -0.21938393439552029, 3e-15),
        (50.5, 1.727631956029118e20, 3e-15),
        (-50.5, -2.2723713293221935e-24, 3e-15),
    ] {
        assert_close(ei(x).unwrap(), expected, 0.0, relative, "Ei anchors");
    }
    for (x, expected) in [
        (1.0, 0.14849550677592205),
        (1.5, 0.07310078653848085),
        (2.0, 0.03753426182049045),
        (-1.0, 0.8231640121031085),
    ] {
        assert_close(e2(x).unwrap(), expected, 0.0, 6e-15, "E2 anchors");
    }
    // Gompertz constant exp(1)*E1(1), DLMF 6.2.1 / 6.9.1.
    assert_close(
        e1_scaled(1.0).unwrap(),
        0.5963473623231941,
        0.0,
        3e-15,
        "scaled E1",
    );
    assert_close(
        e2_scaled(1.0).unwrap(),
        0.4036526376768059,
        0.0,
        3e-15,
        "scaled E2",
    );
    assert_close(
        ei_scaled(1.0).unwrap(),
        0.697174883235066,
        0.0,
        3e-15,
        "scaled Ei",
    );
}

#[test]
fn zero_nonfinite_and_tiny_arguments() {
    for x in [0.0, -0.0] {
        for f in [e1, e1_scaled, ei, ei_scaled] {
            assert_eq!(f(x), Err(Error::Domain));
        }
        assert_eq!(e2(x), Ok(1.0));
        assert_eq!(e2_scaled(x), Ok(1.0));
    }
    for f in [e1, e2, ei, e1_scaled, e2_scaled, ei_scaled] {
        assert_eq!(f(f64::NAN), Err(Error::Domain));
        assert_eq!(f(f64::INFINITY), Err(Error::UnsupportedRange));
        assert_eq!(f(f64::NEG_INFINITY), Err(Error::UnsupportedRange));
    }
    // DLMF 6.6.1-2: omitted correction is O(|x|), far below an ulp.
    for magnitude in [
        f64::from_bits(1),
        f64::from_bits(100),
        f64::MIN_POSITIVE,
        1e-100,
    ] {
        for x in [magnitude, -magnitude] {
            let leading = 0.5772156649015329 + magnitude.ln();
            for f in [e1, e1_scaled] {
                assert_close(f(x).unwrap(), -leading, 0.0, 2e-15, "tiny E1");
            }
            for f in [ei, ei_scaled] {
                assert_close(f(x).unwrap(), leading, 0.0, 2e-15, "tiny Ei");
            }
            assert_eq!(e2(x), Ok(1.0));
            assert_eq!(e2_scaled(x), Ok(1.0));
        }
    }
}

#[test]
fn large_scaled_and_unscaled_range() {
    // DLMF 6.12.1-2 / 8.20.1, first four asymptotic terms at 1e6;
    // the omitted term changes the result by <2e-22 relative.
    for (f, expected) in [
        (
            e1_scaled as fn(f64) -> _,
            (1.0 - 1e-6 + 2e-12 - 6e-18) / 1e6,
        ),
        (e2_scaled, (1.0 - 2e-6 + 6e-12 - 24e-18) / 1e6),
        (ei_scaled, (1.0 + 1e-6 + 2e-12 + 6e-18) / 1e6),
    ] {
        assert_close(f(1e6).unwrap(), expected, 0.0, 8e-16, "large scaled anchor");
    }
    // Real continuation requires a different E2 expansion on the negative axis.
    assert_close(
        e2_scaled(-1e6).unwrap(),
        -(1.0 + 2e-6 + 6e-12 + 24e-18) / 1e6,
        0.0,
        8e-16,
        "negative scaled E2 cancellation regression",
    );
    for x in [1e20, f64::MAX] {
        assert_close(
            e1_scaled(x).unwrap(),
            1.0 / x,
            f64::from_bits(1),
            3e-16,
            "extreme scaled E1",
        );
        assert_close(
            e2_scaled(x).unwrap(),
            1.0 / x,
            f64::from_bits(1),
            3e-16,
            "extreme scaled E2",
        );
        assert_close(
            ei_scaled(x).unwrap(),
            1.0 / x,
            f64::from_bits(1),
            3e-16,
            "extreme scaled Ei",
        );
        assert_close(
            e1_scaled(-x).unwrap(),
            -1.0 / x,
            f64::from_bits(1),
            3e-16,
            "negative extreme E1",
        );
        assert_close(
            e2_scaled(-x).unwrap(),
            -1.0 / x,
            f64::from_bits(1),
            3e-16,
            "negative extreme E2",
        );
        assert_close(
            ei_scaled(-x).unwrap(),
            -1.0 / x,
            f64::from_bits(1),
            3e-16,
            "negative extreme Ei",
        );
    }
    for x in [1000.0, f64::MAX] {
        assert_eq!(e1(x).unwrap().to_bits(), 0.0_f64.to_bits());
        assert_eq!(e2(x).unwrap().to_bits(), 0.0_f64.to_bits());
        assert_eq!(ei(-x).unwrap().to_bits(), (-0.0_f64).to_bits());
        assert_eq!(ei(x), Err(Error::Overflow));
        assert_eq!(e1(-x), Err(Error::Overflow));
        assert_eq!(e2(-x), Err(Error::Overflow));
    }
    // Regressions: exp(710) overflows, but these products remain finite.
    for x in [710.0, 716.0] {
        assert!(ei(x).unwrap().is_finite());
        assert!(e1(-x).unwrap().is_finite());
        assert!(e2(-x).unwrap().is_finite());
        assert!(e1(x).unwrap() > 0.0 && e1(x).unwrap() < f64::MIN_POSITIVE);
        assert!(e2(x).unwrap() > 0.0 && e2(x).unwrap() < f64::MIN_POSITIVE);
    }
    // Independent asymptotic plus elementary exponential, not the production
    // kernel: at 710, terms through x^-6 leave <6e-17 relative error.
    let x = 710.0_f64;
    let polynomial = 1.0
        + 1.0 / x
        + 2.0 / x.powi(2)
        + 6.0 / x.powi(3)
        + 24.0 / x.powi(4)
        + 120.0 / x.powi(5)
        + 720.0 / x.powi(6);
    let expected = (x / 2.0).exp() / x * polynomial * (x / 2.0).exp();
    assert_close(
        ei(x).unwrap(),
        expected,
        0.0,
        2e-15,
        "finite exponential overflow regression",
    );
    assert!(e1(737.0).unwrap() > 0.0);
    assert_eq!(e1(740.0), Ok(0.0));
}

#[test]
fn root_centering_and_branch_neighbors() {
    // DLMF 6.13.1 root 0.372507410781366634461991866580...
    // Its distance from the binary64 input is retained in a second part.
    // Ei(r+d)=exp(r)/r*d+O(d^2); the quadratic term here is <4e-32.
    let root = 0.3725074107813666_f64;
    for bits in [root.to_bits() - 1, root.to_bits(), root.to_bits() + 1] {
        let x = f64::from_bits(bits);
        let distance = (x - root) - 1.3140183414386029e-17;
        let expected = root.exp() / root * distance;
        assert_close(ei(x).unwrap(), expected, 5e-32, 3e-15, "Ei root neighbors");
        assert_close(
            e1(-x).unwrap(),
            -expected,
            5e-32,
            3e-15,
            "E1 continuation root",
        );
    }
    for boundary in [0.25_f64, 0.5, 1.0, 50.0, 700.0, 1400.0] {
        for sign in [-1.0, 1.0] {
            for f in [e1_scaled, e2_scaled, ei_scaled] {
                let center = f(sign * boundary).unwrap();
                for bits in [boundary.to_bits() - 1, boundary.to_bits() + 1] {
                    let value = f(sign * f64::from_bits(bits)).unwrap();
                    assert_close(value, center, 2e-15, 2e-14, "scaled branch neighbors");
                }
            }
        }
    }
    // The unscaled paths also switch exponential reconstruction at 700
    // and 1400. Adjacent inputs must retain subnormals and error classes.
    for boundary in [700.0_f64, 1400.0] {
        for bits in [
            boundary.to_bits() - 1,
            boundary.to_bits(),
            boundary.to_bits() + 1,
        ] {
            let x = f64::from_bits(bits);
            if boundary == 700.0 {
                for f in [e1, e2] {
                    assert_close(
                        f(x).unwrap(),
                        f(boundary).unwrap(),
                        f64::from_bits(2),
                        3e-13,
                        "unscale boundary",
                    );
                }
                assert!(ei(x).unwrap().is_finite());
                assert!(e2(-x).unwrap().is_finite());
            } else {
                assert_eq!(e1(x), Ok(0.0));
                assert_eq!(e2(x), Ok(0.0));
                assert_eq!(ei(x), Err(Error::Overflow));
                assert_eq!(e1(-x), Err(Error::Overflow));
                assert_eq!(e2(-x), Err(Error::Overflow));
            }
        }
    }
    // DLMF 8.19.8, evaluated as an independent E2 power series,
    // rather than the production continuation through Ei. The tail after
    // 40 terms is <1e-42 here. Cancellation near the E2 zero needs an
    // absolute allowance for rounding of gamma and log.
    let y = 1.347_f64;
    let mut expected = 1.0 + y * (1.0 - 0.5772156649015329 - y.ln());
    let mut power_over_factorial = y;
    for k in 2..=40 {
        power_over_factorial *= y / f64::from(k);
        expected -= power_over_factorial / f64::from(k - 1);
    }
    assert_close(
        e2(-y).unwrap(),
        expected,
        2e-15,
        0.0,
        "negative E2 near root",
    );
}

#[test]
fn deterministic_identities() {
    for x in [
        -50.5_f64, -10.0, -2.5, -0.4, -0.01, 0.01, 0.4, 1.0, 2.5, 10.0, 50.5,
    ] {
        assert_close(
            e1(x).unwrap(),
            -ei(-x).unwrap(),
            0.0,
            3e-15,
            "E1/Ei reflection",
        );
        // This grid is well conditioned for the E2 recurrence (unlike huge x).
        assert_close(
            e2(x).unwrap(),
            (-x).exp() - x * e1(x).unwrap(),
            2e-15,
            5e-13,
            "E2 recurrence",
        );
        for (plain, scaled, exponent) in [
            (e1 as fn(f64) -> _, e1_scaled as fn(f64) -> _, x),
            (e2, e2_scaled, x),
            (ei, ei_scaled, -x),
        ] {
            assert_close(
                scaled(x).unwrap(),
                exponent.exp() * plain(x).unwrap(),
                0.0,
                3e-15,
                "scaling",
            );
        }
    }
}

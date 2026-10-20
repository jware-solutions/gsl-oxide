//! Gamma and beta public contracts and numerical anchors.
use gsl_oxide_special::{
    beta::{beta, ln_beta},
    gamma::{gamma, ln_gamma},
    Error,
};
mod support;
use support::assert_close;

#[test]
fn published_gamma_and_log_values() {
    // Numerical examples (not generated here):
    // https://mpmath.readthedocs.io/en/latest/functions/gamma.html#gamma-function
    // and #loggamma. 5e-14 covers recurrence/log rounding, 15-digit anchors.
    for (x, expected) in [
        (0.5, std::f64::consts::PI.sqrt()),
        (1.5, 0.886226925452758),
        (3.5, 3.32335097044784),
        (3.0_f64.sqrt(), 0.9151022969730863),
    ] {
        assert_close(gamma(x).unwrap(), expected, 0.0, 5e-14, "gamma anchor");
    }
    for (x, l, s) in [
        (13.2, 20.494_004_194_566_04, 1),
        (-0.5, 1.2655121234846454, -1),
        (-1.5, 0.860047015376481, 1),
        (-2.5, -0.05624371649767405, -1),
        (1e30, 6.807755278982137e31, 1),
        (1e300, 6.897755278982137e302, 1),
    ] {
        let (v, sign) = ln_gamma(x).unwrap();
        assert_eq!(sign, s);
        assert_close(v, l, 3e-14, 5e-15, "log gamma anchor");
    }
    // Exact factorials and local derivative -Euler at 1 (DLMF 5.4.12).
    for (x, y) in [(1.0, 1.0), (2.0, 1.0), (3.0, 2.0), (6.0, 120.0)] {
        assert_close(gamma(x).unwrap(), y, 0.0, 5e-14, "factorial");
    }
    for base in [1.0_f64, 2.0] {
        for bits in [base.to_bits() - 1, base.to_bits(), base.to_bits() + 1] {
            let x = f64::from_bits(bits);
            let slope = if base == 1.0 {
                -0.5772156649015329
            } else {
                1.0 - 0.5772156649015329
            };
            assert_close(
                ln_gamma(x).unwrap().0,
                slope * (x - base),
                2e-31,
                4e-15,
                "near one/two",
            );
        }
    }
}

#[test]
fn gamma_reflection_poles_and_range() {
    for x in [-10.25, -3.75, -0.25, 0.25, 2.75, 20.5] {
        assert_close(
            gamma(x + 1.0).unwrap(),
            x * gamma(x).unwrap(),
            0.0,
            8e-14,
            "recurrence",
        );
    }
    for x in [-4.25, -2.5, -0.25, 0.25, 0.75] {
        assert_close(
            gamma(x).unwrap() * gamma(1.0 - x).unwrap(),
            std::f64::consts::PI / (std::f64::consts::PI * x).sin(),
            0.0,
            1e-13,
            "reflection",
        );
    }
    for k in [1_u64, 2, 3, 16, 100] {
        let pole = -(k as f64);
        for bits in [pole.to_bits() - 1, pole.to_bits() + 1] {
            let x = f64::from_bits(bits);
            // Residue (-1)^k/k!; relative correction O(distance*log(k)).
            let factorial = (1..=k).fold(1.0, |s, j| s * j as f64);
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            assert_close(
                gamma(x).unwrap() * (x - pole),
                sign / factorial,
                0.0,
                3e-13,
                "pole residue",
            );
        }
    }
    for x in [
        0.0,
        -0.0,
        -1.0,
        -100.0,
        -f64::MAX,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        assert_eq!(gamma(x), Err(Error::Domain));
        assert_eq!(ln_gamma(x), Err(Error::Domain));
    }
    assert!(gamma(171.5).unwrap().is_finite());
    assert!(gamma(171.624).unwrap().is_finite());
    assert_eq!(gamma(171.625), Err(Error::Overflow));
    assert_eq!(gamma(172.0), Err(Error::Overflow));
    assert!(ln_gamma(172.0).unwrap().0 > 710.0);
    assert_eq!(ln_gamma(f64::MAX), Err(Error::Overflow));
    assert_eq!(gamma(-180.5).unwrap().to_bits(), (-0.0_f64).to_bits());
    assert!(ln_gamma(-180.5).unwrap().0 < -745.0);
    // Half-integer recurrence from Gamma(-1/2)=-2sqrt(pi), evaluated by
    // successive division. The final reference is subnormal; 171 divisions
    // allow 4e-13 relative, with no absolute tolerance accepting zero.
    let mut half = -2.0 * std::f64::consts::PI.sqrt();
    for k in 1..=171 {
        half /= -(f64::from(k) + 0.5);
    }
    let subnormal = gamma(-171.5).unwrap();
    assert!(subnormal > 0.0 && subnormal < f64::MIN_POSITIVE);
    assert_close(subnormal, half, 0.0, 4e-13, "subnormal gamma half-integer");
    assert!(gamma(1e-300).unwrap() > 9.99e299);
    assert_close(
        ln_gamma(-f64::from_bits(1)).unwrap().0,
        -f64::from_bits(1).ln(),
        0.0,
        2e-15,
        "tiny negative log",
    );
}

#[test]
fn beta_independent_anchors_and_extremes() {
    // Published beta examples, same gamma.html#beta-function URL above.
    for (a, b, y) in [
        (5.0, 2.0, 1.0 / 30.0),
        (1.5, 2.0, 4.0 / 15.0),
        (2.5, 1.5, std::f64::consts::PI / 16.0),
        (
            std::f64::consts::PI,
            std::f64::consts::E,
            0.0378902987812122,
        ),
        (0.5, 0.5, std::f64::consts::PI),
    ] {
        assert_close(beta(a, b).unwrap(), y, 0.0, 6e-14, "beta anchor");
        assert_close(
            ln_beta(a, b).unwrap(),
            y.ln(),
            6e-14,
            1e-14,
            "log beta anchor",
        );
    }
    for b in [
        f64::from_bits(1),
        1e-300,
        0.01,
        1.0,
        1000.0,
        1e100,
        f64::MAX,
    ] {
        assert_close(
            ln_beta(1.0, b).unwrap(),
            -b.ln(),
            2e-13,
            2e-15,
            "B(1,b) log",
        );
        if (1.0 / b).is_finite() {
            assert_close(beta(1.0, b).unwrap(), 1.0 / b, 0.0, 1e-13, "B(1,b)");
        }
    }
    // B(2,b)=1/(b(b+1)); even gamma factors and a+b may overflow.
    for b in [200.0_f64, 1e10, 1e100, f64::MAX] {
        assert_close(
            ln_beta(2.0, b).unwrap(),
            -b.ln() - (b + 1.0).ln(),
            2e-12,
            3e-15,
            "unequal beta",
        );
    }
    for (a, b) in [
        (0.01, 0.02),
        (16.0, 17.0),
        (200.0, 300.0),
        (1e10, 1e12),
        (1e308, 1e308),
    ] {
        assert_eq!(ln_beta(a, b), ln_beta(b, a));
        assert!(ln_beta(a, b).unwrap().is_finite());
    }
    assert!(beta(200.0, 3.0).unwrap() > 0.0);
    assert_eq!(beta(1000.0, 1000.0).unwrap(), 0.0);
    assert!(ln_beta(1000.0, 1000.0).unwrap() < -1000.0);
    assert_eq!(beta(1e-310, 1e-310), Err(Error::Overflow));
    assert_eq!(ln_beta(f64::MAX, f64::MAX), Err(Error::Overflow));
    for x in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(beta(x, 1.0), Err(Error::Domain));
        assert_eq!(ln_beta(1.0, x), Err(Error::Domain));
    }
}

#[test]
fn gamma_and_beta_branch_neighbors() {
    for boundary in [0.25_f64, 0.75, 1.25, 1.75, 2.25, 16.0] {
        let center = gamma(boundary).unwrap();
        for bits in [boundary.to_bits() - 1, boundary.to_bits() + 1] {
            let x = f64::from_bits(bits);
            assert_close(
                gamma(x).unwrap(),
                center,
                0.0,
                6e-14,
                "gamma branch neighbors",
            );
            // Independent B(x,2)=1/(x(x+1)) at beta's recurrence boundaries.
            let expected = 1.0 / (x * (x + 1.0));
            assert_close(
                beta(x, 2.0).unwrap(),
                expected,
                0.0,
                8e-14,
                "beta branch neighbors",
            );
            assert_close(
                ln_beta(x, 2.0).unwrap(),
                expected.ln(),
                8e-14,
                0.0,
                "log beta branch neighbors",
            );
        }
    }
}

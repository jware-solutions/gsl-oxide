//! Complementary tails, numerical anchors, boundaries and domains.
use gsl_oxide_special::{
    error_functions::{erf, erfc},
    incomplete_beta::{complement, regularized},
    incomplete_gamma::{p, q},
    Error,
};
mod support;
use support::assert_close;

#[test]
fn erfc_anchors_and_tail_bounds() {
    // Published erfc(10) example and erf(0.5) example:
    // https://mpmath.readthedocs.io/en/latest/functions/expintegrals.html#erfc
    for (x, y) in [(0.5, 0.479_500_122_186_953_5), (10.0, 2.08848758376254e-45)] {
        assert_close(erfc(x), y, 0.0, 3e-14, "erfc anchor");
    }
    // Independent first-omitted-term enclosure from DLMF 7.12.1.
    // At 27 the interval still excludes zero: this tests a subnormal tail.
    for x in [6.0_f64, 10.0, 20.0, 26.0, 26.5, 27.0, 27.2] {
        let upper = (-x * x).exp() / (x * std::f64::consts::PI.sqrt());
        let lower = upper * (1.0 - 0.5 / (x * x));
        let y = erfc(x);
        let step = f64::from_bits(1);
        assert!(
            y > 0.0 && y >= lower - step && y <= upper + step,
            "tail x={x} y={y}"
        );
    }
    for x in [
        0.0,
        f64::from_bits(1),
        1e-20,
        0.5,
        1.0,
        3.0,
        6.0,
        27.0,
        28.0,
        f64::MAX,
    ] {
        assert_close(erfc(-x), 2.0 - erfc(x), 0.0, 2e-15, "erfc symmetry");
        assert_close(erf(x) + erfc(x), 1.0, 4e-15, 0.0, "erf complement");
    }
    for boundary in [1.0_f64, 28.0] {
        for x in [
            f64::from_bits(boundary.to_bits() - 1),
            boundary,
            f64::from_bits(boundary.to_bits() + 1),
        ] {
            if boundary == 1.0 {
                assert_close(erfc(x), 1.0 - erf(x), 3e-15, 0.0, "erfc boundary");
            } else {
                assert_eq!(erfc(x), 0.0);
            }
        }
    }
    for x in [0.0, -0.0] {
        assert_eq!(erfc(x), 1.0);
    }
    assert_eq!(erfc(f64::INFINITY).to_bits(), 0.0_f64.to_bits());
    assert_eq!(erfc(f64::NEG_INFINITY), 2.0);
    assert!(erfc(f64::NAN).is_nan());
    assert_eq!(erfc(27.5), 0.0);
}

#[test]
fn incomplete_gamma_anchors_and_closed_forms() {
    // Published unregularized values divided by known Gamma(10)=9!,
    // Gamma(4.5)=105sqrt(pi)/16, Gamma(0.5)=sqrt(pi).
    // https://mpmath.readthedocs.io/en/latest/functions/expintegrals.html#gammainc
    assert_close(
        q(10.0, 100.0).unwrap(),
        4.083_660_630_910_611e-26 / 362880.0,
        0.0,
        5e-14,
        "upper gamma anchor",
    );
    let g = 105.0 * std::f64::consts::PI.sqrt() / 16.0;
    assert_close(
        p(4.5, 2.0).unwrap(),
        1.0304254272321143 / g,
        0.0,
        5e-14,
        "lower gamma anchor",
    );
    assert_close(
        q(4.5, 2.0).unwrap(),
        10.601302969335335 / g,
        0.0,
        5e-14,
        "upper gamma anchor",
    );
    assert_close(
        p(0.5, 2.0).unwrap(),
        1.6918067329451983 / std::f64::consts::PI.sqrt(),
        0.0,
        5e-14,
        "half gamma anchor",
    );
    for x in [1e-300_f64, 1e-10, 0.1, 1.0, 5.0, 100.0, 700.0] {
        assert_close(
            p(1.0, x).unwrap(),
            -(-x).exp_m1(),
            0.0,
            2e-15,
            "P exponential",
        );
        assert_close(q(1.0, x).unwrap(), (-x).exp(), 0.0, 2e-15, "Q exponential");
        assert_close(
            p(0.5, x * x).unwrap(),
            erf(x),
            1e-14,
            4e-14,
            "half gamma erf",
        );
        assert_close(
            q(0.5, x * x).unwrap(),
            erfc(x),
            0.0,
            2e-13,
            "half gamma erfc",
        );
    }
    for n in [2, 5, 20, 100] {
        for x in [0.01_f64, 1.0, 10.0, 100.0] {
            let mut term = 1.0;
            let mut sum = term;
            for k in 1..n {
                term *= x / f64::from(k);
                sum += term;
            }
            assert_close(
                q(f64::from(n), x).unwrap(),
                (-x).exp() * sum,
                0.0,
                2e-13,
                "integer gamma",
            );
        }
    }
    // Tiny Q for tiny shape: Q(a,1)/a -> E1(1), with O(a) correction.
    // Published E1(1) https://mpmath.readthedocs.io/en/latest/functions/expintegrals.html#e1
    assert_close(
        q(1e-20, 1.0).unwrap(),
        2.193_839_343_955_203e-21,
        0.0,
        2e-14,
        "small-shape Q",
    );
    assert!(p(2.0, 1e-150).unwrap() > 0.0);
    assert_close(p(2.0, 1e-150).unwrap(), 5e-301, 0.0, 2e-13, "tiny P");
}

#[test]
fn incomplete_gamma_regions_and_errors() {
    for a in [0.1_f64, 0.25, 0.5, 0.99, 1.01, 10.0, 16.0, 1000.0, 1e6] {
        for x in [
            1e-10,
            a * 0.9,
            a,
            a + 1.0,
            f64::from_bits((a + 1.0).to_bits() - 1),
            f64::from_bits((a + 1.0).to_bits() + 1),
            a * 1.1,
            1e308,
        ] {
            let (v, w) = (p(a, x).unwrap(), q(a, x).unwrap());
            assert_close(v + w, 1.0, 3e-15, 0.0, "gamma complement");
            assert!((0.0..=1.0).contains(&v) && (0.0..=1.0).contains(&w));
        }
        assert_eq!(p(a, 0.0), Ok(0.0));
        assert_eq!(q(a, -0.0), Ok(1.0));
        assert_eq!(p(a, f64::INFINITY), Ok(1.0));
        assert_eq!(q(a, f64::INFINITY), Ok(0.0));
    }
    // Central-limit expansion at x=a: P=1/2+1/(3sqrt(2pi a))+O(a^-3/2),
    // DLMF 8.12. The first correction is 1/(540 sqrt(2pi) a^1.5).
    assert_close(
        p(1e6, 1e6).unwrap(),
        0.5001329807608725,
        3e-13,
        0.0,
        "large transition",
    );
    for a in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(p(a, 0.0), Err(Error::Domain));
        assert_eq!(q(a, 1.0), Err(Error::Domain));
    }
    for x in [-1.0, f64::NEG_INFINITY, f64::NAN] {
        assert_eq!(q(1.0, x), Err(Error::Domain));
    }
    assert_eq!(p(1e6 + 1.0, 1.0), Err(Error::UnsupportedRange));
}

#[test]
fn incomplete_beta_closed_forms_and_anchors() {
    // Published I_.25(1,2)=.4375:
    // https://mpmath.readthedocs.io/en/latest/functions/gamma.html#betainc
    assert_close(
        regularized(1.0, 2.0, 0.25).unwrap(),
        0.4375,
        0.0,
        2e-15,
        "published beta CDF",
    );
    // Polynomial integrals, independent exact rationals (DLMF 8.17.5).
    for (a, b, x, y) in [
        (2.0, 3.0, 0.5, 11.0 / 16.0),
        (2.0, 2.0, 0.25, 5.0 / 32.0),
        (3.0, 4.0, 0.5, 21.0 / 32.0),
        (0.5, 0.5, 0.25, 1.0 / 3.0),
    ] {
        assert_close(
            regularized(a, b, x).unwrap(),
            y,
            0.0,
            8e-14,
            "beta integral anchor",
        );
        assert_close(
            complement(a, b, x).unwrap(),
            1.0 - y,
            0.0,
            8e-14,
            "beta upper anchor",
        );
    }
    for x in [
        f64::from_bits(1),
        1e-300,
        1e-10,
        0.1,
        0.5,
        0.9,
        f64::from_bits(1.0_f64.to_bits() - 1),
    ] {
        assert_close(regularized(1.0, 1.0, x).unwrap(), x, 0.0, 3e-15, "uniform");
        for a in [0.1_f64, 2.0, 1000.0] {
            assert_close(
                regularized(a, 1.0, x).unwrap(),
                x.powf(a),
                0.0,
                4e-14,
                "power CDF",
            );
            assert_close(
                complement(1.0, a, x).unwrap(),
                (a * (-x).ln_1p()).exp(),
                0.0,
                2e-15,
                "power survival",
            );
        }
    }
    // B(2,3): I=6x²-8x³+3x⁴. Tiny nonzero reference cannot accept zero.
    assert_close(
        regularized(2.0, 3.0, 1e-150).unwrap(),
        6e-300,
        0.0,
        3e-13,
        "tiny beta lower",
    );
    let x = f64::from_bits(1.0_f64.to_bits() - 1);
    let d = 1.0 - x;
    assert_close(
        complement(3.0, 2.0, x).unwrap(),
        6.0 * d * d - 8.0 * d * d * d + 3.0 * d.powi(4),
        0.0,
        1e-13,
        "tiny beta upper",
    );
}

#[test]
fn published_beta_tails_and_large_shape_anchor() {
    // Published scalar examples, 5e-14 relative for log/CF rounding:
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betainc.html
    for (a, b, x, y) in [
        (1.4, 3.1, 0.5, 0.8148904036225296),
        (2.2, 3.1, 0.4, 0.49339638807619446),
    ] {
        assert_close(
            regularized(a, b, x).unwrap(),
            y,
            0.0,
            5e-14,
            "published beta lower",
        );
        assert_close(
            complement(a, b, x).unwrap(),
            1.0 - y,
            0.0,
            5e-14,
            "published beta upper",
        );
    }
    // Published tiny upper tails (9 significant digits, relative 3e-9):
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.betaincc.html
    for (x, y) in [
        (0.9, 2.05746321e-9),
        (0.99, 1.97259354e-17),
        (0.999, 1.96467954e-25),
    ] {
        assert_close(
            complement(0.5, 8.0, x).unwrap(),
            y,
            0.0,
            3e-9,
            "published tiny beta upper",
        );
    }
    // Exact binomial identity I_.5(n,n+1)=.5+.5*C(2n,n)/4^n
    // (DLMF 8.17.5). Central-binomial Stirling expansion, DLMF 5.11,
    // through 1/(128n²) has <1e-22 absolute truncation here. 1e-9 absolute
    // allows the documented large balanced log-factor cancellation.
    let n = 999_999.0;
    let expected = 0.5
        + (1.0 - 1.0 / (8.0 * n) + 1.0 / (128.0 * n * n))
            / (2.0 * (std::f64::consts::PI * n).sqrt());
    assert_close(
        regularized(n, n + 1.0, 0.5).unwrap(),
        expected,
        1e-9,
        0.0,
        "large beta binomial anchor",
    );
    assert_close(
        complement(n, n + 1.0, 0.5).unwrap(),
        1.0 - expected,
        1e-9,
        0.0,
        "large beta complement anchor",
    );
    for a in [1000.0_f64, 1e6] {
        let x = 1.0 - 1.0 / a;
        assert_close(
            regularized(a, 2.0, x).unwrap(),
            x.powf(a) * (1.0 + a * (1.0 - x)),
            0.0,
            2e-12,
            "unequal large beta closed form",
        );
    }
}

#[test]
fn incomplete_beta_regions_and_errors() {
    for (a, b) in [
        (0.1, 0.2),
        (2.0, 30.0),
        (1e5, 2.0),
        (1000.0, 1000.0),
        (1e6, 1e6),
    ] {
        let boundary: f64 = (a + 1.0) / (a + b + 2.0);
        let mut previous = 0.0;
        for x in [0.0, 1e-10, 0.1, 0.5, 0.9, 1.0] {
            let y = regularized(a, b, x).unwrap();
            assert!(y >= previous);
            previous = y;
            assert_close(
                y + complement(a, b, x).unwrap(),
                1.0,
                3e-15,
                0.0,
                "beta complement",
            );
            if x > 1e-9 && x < 1.0 {
                assert_close(
                    y,
                    complement(b, a, 1.0 - x).unwrap(),
                    2e-9,
                    0.0,
                    "beta symmetry",
                );
            }
        }
        let center = regularized(a, b, boundary).unwrap();
        for bits in [boundary.to_bits() - 1, boundary.to_bits() + 1] {
            assert_close(
                regularized(a, b, f64::from_bits(bits)).unwrap(),
                center,
                2e-9,
                0.0,
                "beta boundary",
            );
        }
    }
    for a in [0.1, 2.0, 1e6] {
        assert_eq!(regularized(a, a, 0.5), Ok(0.5));
    }
    for a in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(regularized(a, 1.0, 0.0), Err(Error::Domain));
        assert_eq!(complement(1.0, a, 1.0), Err(Error::Domain));
    }
    for x in [-0.1, 1.1, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(regularized(2.0, 3.0, x), Err(Error::Domain));
    }
    assert_eq!(
        regularized(1e6 + 1.0, 1.0, 0.0),
        Err(Error::UnsupportedRange)
    );
}

#[test]
fn tiny_shape_tail_regressions() {
    // Integrating t^(a-1)(1-t) gives I_x(a,2)=x^a[1+a(1-x)].
    let a = 1e-20;
    let x = 0.25_f64;
    let expected = -(a * x.ln()).exp_m1() - a * (1.0 - x) * (a * x.ln()).exp();
    assert_close(
        complement(a, 2.0, x).unwrap(),
        expected,
        0.0,
        2e-14,
        "tiny shape beta Q",
    );
    assert_eq!(regularized(a, 2.0, x), Ok(1.0));
    assert_close(
        regularized(2.0, a, 1.0 - x).unwrap(),
        expected,
        0.0,
        2e-14,
        "tiny shape beta symmetric P",
    );
    // x^a/Gamma(a+1) approximates gamma P at x=1e-300 to O(x).
    // Use a=1/8 and compare two lower tails after changing x by 2^8:
    // the ratio is 2, with O(x) error, and both must remain positive.
    let p1 = p(0.125, 1e-300).unwrap();
    assert!(p1 > 0.0);
    assert_close(
        p(0.125, 256e-300).unwrap(),
        2.0 * p1,
        0.0,
        4e-14,
        "small shape small gamma P",
    );
    let expected = 1e-300_f64.powf(0.125) * (1.0 + 0.125);
    assert_close(
        regularized(0.125, 2.0, 1e-300).unwrap(),
        expected,
        0.0,
        3e-14,
        "small shape tiny beta P",
    );
    // The old lower-series switch at x=1 must retain the same tiny Q
    // immediately above 1. Q(a,x)/a -> E1(x); |E1(1+h)-E1(1)| < h.
    let x = f64::from_bits(1.0_f64.to_bits() + 1);
    assert_close(
        q(1e-10, x).unwrap() / 1e-10,
        0.21938393439552027,
        5e-11,
        0.0,
        "tiny Q above one",
    );
}

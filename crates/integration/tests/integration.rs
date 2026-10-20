//! Independently written known-answer tests from elementary antiderivatives.
use gsl_oxide_integration::{integrate, integrate_with_breakpoints, Error};
use std::f64::consts::PI;
mod support;
use support::{assert_close, check};

#[test]
fn elementary_integrals_and_estimates() {
    // Integral x^n on [0,1] = 1/(n+1); sin and exp from their antiderivatives.
    for n in 0..=24 {
        check(
            integrate(|x| x.powi(n), 0.0, 1.0, 1e-12, 1e-12, 20_000),
            1.0 / f64::from(n + 1),
            1e-12,
            1e-12,
        );
    }
    for (f, b, expected) in [
        (f64::sin as fn(f64) -> f64, PI, 2.0),
        (f64::exp, 1.0, std::f64::consts::E - 1.0),
    ] {
        let (value, error) = check(
            integrate(f, 0.0, b, 1e-12, 1e-12, 20_000),
            expected,
            1e-12,
            1e-12,
        );
        assert!((value - expected).abs() <= error + 8.0 * f64::EPSILON * expected);
    }
    check(
        integrate(|x| x.exp(), 0.0, 1.0, 0.0, 1e-12, 10_000),
        std::f64::consts::E - 1.0,
        0.0,
        1e-12,
    );
    check(
        integrate(|_| 1.0, 0.0, 1.0, 0.0, 64.0 * f64::EPSILON, 15),
        1.0,
        0.0,
        64.0 * f64::EPSILON,
    );
}

#[test]
fn transformations_reversal_and_additivity() {
    let (area, error) = check(
        integrate(|x| (x - 7.0).powi(4), 3.0, 9.0, 1e-11, 1e-12, 10_000),
        (4.0_f64.powi(5) + 2.0_f64.powi(5)) / 5.0,
        1e-11,
        1e-12,
    );
    let reverse = integrate(|x| (x - 7.0).powi(4), 9.0, 3.0, 1e-11, 1e-12, 10_000).unwrap();
    assert_eq!(reverse, (-area, error));
    let left = integrate(f64::exp, -1.0, 0.3, 1e-12, 1e-12, 10_000)
        .unwrap()
        .0;
    let right = integrate(f64::exp, 0.3, 2.0, 1e-12, 1e-12, 10_000)
        .unwrap()
        .0;
    let all = integrate(f64::exp, -1.0, 2.0, 1e-12, 1e-12, 10_000)
        .unwrap()
        .0;
    assert_close(left + right, all, 2e-12, 1e-12);
    for (a, b) in [(-8.0, 3.0), (2.0, 2.1), (1000.0, 1004.0)] {
        check(
            integrate(|x| ((x - a) / (b - a)).powi(6), a, b, 1e-12, 1e-12, 10_000),
            (b - a) / 7.0,
            1e-12,
            1e-12,
        );
    }
}

#[test]
fn cancellation_and_global_breakpoint_accuracy() {
    check(
        integrate(|x| x.powi(9), -1.0, 1.0, 1e-13, 0.0, 1000),
        0.0,
        1e-13,
        0.0,
    );
    check(
        integrate(|x| 1e6 * x + 0.125, -1.0, 1.0, 1e-8, 0.0, 1000),
        0.25,
        1e-8,
        0.0,
    );
    let f = |x: f64| {
        if x < 0.0 {
            1e6 + x.exp()
        } else {
            -1e6 + x.exp()
        }
    };
    check(
        integrate_with_breakpoints(f, -1.0, 1.0, &[0.0], 2e-8, 0.0, 1000),
        std::f64::consts::E - (-1.0_f64).exp(),
        2e-8,
        0.0,
    );
    // Locally relative-accurate large segments cannot certify a tiny residual.
    assert_eq!(
        integrate_with_breakpoints(
            |x| if x < 0.0 { 1e12 } else { -1e12 + 1.0 },
            -1.0,
            1.0,
            &[0.0],
            0.0,
            1e-12,
            1000
        ),
        Err(Error::Roundoff)
    );
    check(
        integrate_with_breakpoints(
            |x| if x < 0.0 { 1e12 } else { -1e12 + 1.0 },
            -1.0,
            1.0,
            &[0.0],
            0.1,
            0.0,
            1000,
        ),
        1.0,
        0.1,
        0.0,
    );
}

#[test]
fn discontinuity_and_unequal_segment_difficulty() {
    for (a, b, sign) in [(0.0, 1.0, 1.0), (1.0, 0.0, -1.0)] {
        check(
            integrate_with_breakpoints(
                |x| if x < 0.3 { 1.0 } else { 4.0 },
                a,
                b,
                &[0.3],
                1e-12,
                0.0,
                30,
            ),
            sign * 3.1,
            1e-12,
            0.0,
        );
    }
    // Most work is on the exponential segment. The constant segment shares
    // the total target and budget rather than getting its own allowance.
    let expected = 1.0 + (20.0_f64.exp() - 1.0) / 20.0;
    check(
        integrate_with_breakpoints(
            |x| if x < 0.0 { 1.0 } else { (20.0 * x).exp() },
            -1.0,
            1.0,
            &[0.0],
            1e-8,
            1e-12,
            20_000,
        ),
        expected,
        1e-8,
        1e-12,
    );
}

#[test]
fn peaks_with_analytic_arctangent_answers() {
    for (width, points) in [(0.03, &[][..]), (1e-6, &[0.29999, 0.3, 0.30001][..])] {
        let expected = ((1.0_f64 - 0.3) / width).atan() - ((0.0_f64 - 0.3) / width).atan();
        check(
            integrate_with_breakpoints(
                |x| width / ((x - 0.3).powi(2) + width * width),
                0.0,
                1.0,
                points,
                1e-10,
                1e-11,
                100_000,
            ),
            expected,
            1e-10,
            1e-11,
        );
    }
}

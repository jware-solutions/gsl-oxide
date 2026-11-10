//! Endpoint singularity accuracy, conditioning limits, and failure regressions.
use gsl_oxide_integration::{integrate_singular, Error};
use std::f64::consts::PI;
mod support;
use support::check;

/// Checks endpoint roots and logarithms.
#[test]
fn endpoint_roots_and_logarithms() {
    // Elementary antiderivatives: 2sqrt(x), x*log(x)-x; the beta(1/2,1/2)
    // integral equals pi by x=sin(theta)^2. No numerical oracle is used.
    for (a, b, sign) in [(0.0, 1.0, 1.0), (1.0, 0.0, -1.0)] {
        check(
            integrate_singular(|x| 1.0 / x.sqrt(), a, b, 1e-11, 1e-11, 50_000),
            2.0 * sign,
            1e-11,
            1e-11,
        );
        check(
            integrate_singular(f64::ln, a, b, 1e-11, 1e-11, 50_000),
            -sign,
            1e-11,
            1e-11,
        );
        // A one-argument f64 callback loses endpoint distance near x=1.
        // Tail and mapping estimates retain this roughly sqrt(epsilon) limit.
        check(
            integrate_singular(|x| 1.0 / (1.0 - x).sqrt(), a, b, 2e-7, 0.0, 50_000),
            2.0 * sign,
            2e-7,
            0.0,
        );
        check(
            integrate_singular(|x| 1.0 / (x * (1.0 - x)).sqrt(), a, b, 3e-7, 0.0, 50_000),
            PI * sign,
            3e-7,
            0.0,
        );
        check(
            integrate_singular(|x| (1.0 - x).ln(), a, b, 1e-10, 0.0, 50_000),
            -sign,
            1e-10,
            0.0,
        );
    }
}

/// Checks nonsingular inputs and shifted endpoint.
#[test]
fn nonsingular_inputs_and_shifted_endpoint() {
    check(
        integrate_singular(|x| x * x, 0.0, 1.0, 1e-12, 1e-12, 30_000),
        1.0 / 3.0,
        1e-12,
        1e-12,
    );
    check(
        integrate_singular(f64::sin, 0.0, PI, 1e-11, 1e-11, 30_000),
        2.0,
        1e-11,
        1e-11,
    );
    check(
        integrate_singular(|x| 1.0 / (x - 2.0).sqrt(), 2.0, 6.0, 5e-7, 0.0, 50_000),
        4.0,
        5e-7,
        0.0,
    );
}

/// Checks unresolved endpoint contribution and divergence fail.
#[test]
fn unresolved_endpoint_contribution_and_divergence_fail() {
    assert!(matches!(
        integrate_singular(|x| 1.0 / (1.0 - x).sqrt(), 0.0, 1.0, 1e-12, 0.0, 100_000),
        Err(Error::Roundoff)
    ));
    assert!(integrate_singular(|x| 1.0 / x, 0.0, 1.0, 1e-10, 0.0, 2000).is_err());
}

/// Checks representative estimator quality and tiny interval resolution.
#[test]
fn representative_estimator_quality_and_tiny_interval_resolution() {
    for (f, expected, tolerance) in [
        ((|x: f64| 1.0 / x.sqrt()) as fn(f64) -> f64, 2.0, 1e-11),
        (f64::ln, -1.0, 1e-11),
        (|x| 1.0 / (x * (1.0 - x)).sqrt(), PI, 3e-7),
    ] {
        let (value, error) = check(
            integrate_singular(f, 0.0, 1.0, tolerance, 0.0, 50_000),
            expected,
            tolerance,
            0.0,
        );
        assert!((value - expected).abs() <= error + 32.0 * f64::EPSILON * expected.abs());
    }
    // Rounding a subnormal half width must not silently change the interval.
    let width = f64::from_bits(10_001);
    assert_eq!(
        integrate_singular(|_| 1e308, 0.0, width, 1e-25, 0.0, 10_000),
        Err(Error::Roundoff)
    );
}

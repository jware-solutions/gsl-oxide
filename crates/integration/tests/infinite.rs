//! Known improper integrals from elementary antiderivatives and the Gaussian integral.
use gsl_oxide_integration::{
    integrate_infinite, integrate_lower_infinite, integrate_upper_infinite,
};
use std::f64::consts::PI;
mod support;
use support::check;

/// Checks exponentials gaussians and heavy tails.
#[test]
fn exponentials_gaussians_and_heavy_tails() {
    check(
        integrate_upper_infinite(|x| (-x).exp(), 0.0, 1e-12, 1e-12, 30_000),
        1.0,
        1e-12,
        1e-12,
    );
    check(
        integrate_lower_infinite(f64::exp, 0.0, 1e-12, 1e-12, 30_000),
        1.0,
        1e-12,
        1e-12,
    );
    check(
        integrate_upper_infinite(|x| (-x * x).exp(), 0.0, 1e-12, 1e-12, 30_000),
        PI.sqrt() / 2.0,
        1e-12,
        1e-12,
    );
    check(
        integrate_infinite(|x| (-x * x).exp(), 1e-12, 1e-12, 30_000),
        PI.sqrt(),
        1e-12,
        1e-12,
    );
    check(
        integrate_upper_infinite(|x| 1.0 / (1.0 + x * x), 0.0, 1e-12, 1e-12, 30_000),
        PI / 2.0,
        1e-12,
        1e-12,
    );
    check(
        integrate_infinite(|x| 1.0 / (1.0 + x * x), 1e-12, 1e-12, 30_000),
        PI,
        1e-12,
        1e-12,
    );
}

/// Checks shifted boundaries and independently convergent asymmetric tails.
#[test]
fn shifted_boundaries_and_independently_convergent_asymmetric_tails() {
    check(
        integrate_upper_infinite(|x| (3.0 - x).exp(), 3.0, 1e-12, 1e-12, 30_000),
        1.0,
        1e-12,
        1e-12,
    );
    check(
        integrate_lower_infinite(|x| (x - 2.0).exp(), 2.0, 1e-12, 1e-12, 30_000),
        1.0,
        1e-12,
        1e-12,
    );
    check(
        integrate_upper_infinite(
            |x| 1.0 / (1.0 + (x - 1.0).powi(2)),
            2.0,
            1e-12,
            1e-12,
            30_000,
        ),
        PI / 4.0,
        1e-12,
        1e-12,
    );
    check(
        integrate_lower_infinite(
            |x| 1.0 / (1.0 + (x - 1.0).powi(2)),
            2.0,
            1e-12,
            1e-12,
            30_000,
        ),
        3.0 * PI / 4.0,
        1e-12,
        1e-12,
    );
    check(
        integrate_infinite(
            |x| {
                if x < 0.0 {
                    (2.0 * x).exp()
                } else {
                    3.0 * (-x).exp()
                }
            },
            1e-11,
            1e-12,
            30_000,
        ),
        3.5,
        1e-11,
        1e-12,
    );
    check(
        integrate_infinite(|x| (-(x - 2.0).powi(2)).exp(), 1e-11, 1e-12, 30_000),
        PI.sqrt(),
        1e-11,
        1e-12,
    );
    check(
        integrate_infinite(|x| 1.0 / (1.0 + (x - 3.0).powi(2)), 1e-11, 1e-12, 30_000),
        PI,
        1e-11,
        1e-12,
    );
}

/// Checks opposite divergent tails are not a principal value.
#[test]
fn opposite_divergent_tails_are_not_a_principal_value() {
    for f in [f64::signum as fn(f64) -> f64, |x| x / (1.0 + x * x)] {
        assert!(integrate_infinite(f, 1e-10, 1e-10, 3000).is_err());
    }
    assert!(integrate_upper_infinite(|x| 1.0 / x, 1.0, 1e-12, 1e-12, 1000).is_err());
}

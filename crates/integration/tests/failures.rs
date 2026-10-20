//! Validation order, exact original-callback budgets, and explicit numerical failures.
use gsl_oxide_integration::{
    integrate, integrate_infinite, integrate_lower_infinite, integrate_singular,
    integrate_upper_infinite, integrate_with_breakpoints, Error,
};

type Runner = fn(&mut dyn FnMut(f64) -> f64, f64, f64, usize) -> Result<(f64, f64), Error>;

fn runners() -> [Runner; 6] {
    [
        |f, a, r, n| integrate(f, 0.0, 1.0, a, r, n),
        |f, a, r, n| integrate_with_breakpoints(f, 0.0, 1.0, &[0.25, 0.75], a, r, n),
        |f, a, r, n| integrate_upper_infinite(f, 0.0, a, r, n),
        |f, a, r, n| integrate_lower_infinite(f, 0.0, a, r, n),
        |f, a, r, n| integrate_infinite(f, a, r, n),
        |f, a, r, n| integrate_singular(f, 0.0, 1.0, a, r, n),
    ]
}

#[test]
fn every_method_checks_tolerances_before_callbacks() {
    for run in runners() {
        for (abs, rel) in [
            (0.0, 0.0),
            (0.0, 1e-16),
            (-1.0, 1.0),
            (1.0, -1.0),
            (f64::NAN, 1.0),
            (1.0, f64::NAN),
            (f64::INFINITY, 1.0),
            (1.0, f64::INFINITY),
            (f64::NEG_INFINITY, 1.0),
        ] {
            assert_eq!(
                run(&mut |_| panic!("invalid arguments sampled"), abs, rel, 1000),
                Err(Error::InvalidInput)
            );
        }
        assert_eq!(
            run(&mut |_| panic!("zero budget sampled"), 1e-10, 1e-10, 0),
            Err(Error::EvaluationLimit)
        );
    }
}

#[test]
fn finite_bounds_and_breakpoint_validation_is_complete() {
    for nonfinite in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for (a, b) in [(nonfinite, 1.0), (0.0, nonfinite)] {
            assert_eq!(
                integrate(|_| panic!("endpoint"), a, b, 1e-10, 0.0, 100),
                Err(Error::InvalidInput)
            );
            assert_eq!(
                integrate_singular(|_| panic!("endpoint"), a, b, 1e-10, 0.0, 100),
                Err(Error::InvalidInput)
            );
        }
        assert_eq!(
            integrate_upper_infinite(|_| panic!("endpoint"), nonfinite, 1e-10, 0.0, 100),
            Err(Error::InvalidInput)
        );
        assert_eq!(
            integrate_lower_infinite(|_| panic!("endpoint"), nonfinite, 1e-10, 0.0, 100),
            Err(Error::InvalidInput)
        );
    }
    for points in [
        &[0.0][..],
        &[1.0],
        &[-0.1],
        &[1.1],
        &[0.2, 0.2],
        &[0.8, 0.2],
        &[0.2, f64::NAN],
        &[0.2, f64::INFINITY],
        &[f64::NEG_INFINITY],
    ] {
        for (a, b) in [(0.0, 1.0), (1.0, 0.0)] {
            assert_eq!(
                integrate_with_breakpoints(
                    |_| panic!("breakpoints"),
                    a,
                    b,
                    points,
                    1e-10,
                    0.0,
                    1000
                ),
                Err(Error::InvalidInput)
            );
        }
    }
    assert_eq!(
        integrate_with_breakpoints(|_| panic!("zero interval"), 1.0, 1.0, &[0.5], 1e-10, 0.0, 0),
        Err(Error::InvalidInput)
    );
}

#[test]
fn zero_length_and_empty_breakpoint_contracts() {
    for x in [-f64::MAX, -1.0, -0.0, 0.0, 1.0, f64::MAX] {
        assert_eq!(
            integrate(|_| panic!("zero length"), x, x, 1e-10, 0.0, 0),
            Ok((0.0, 0.0))
        );
        assert_eq!(
            integrate_singular(|_| panic!("zero length"), x, x, 1e-10, 0.0, 0),
            Ok((0.0, 0.0))
        );
        assert_eq!(
            integrate_with_breakpoints(|_| panic!("zero length"), x, x, &[], 1e-10, 0.0, 0),
            Ok((0.0, 0.0))
        );
    }
    assert_eq!(
        integrate(f64::exp, 0.0, 1.0, 1e-12, 1e-12, 1000),
        integrate_with_breakpoints(f64::exp, 0.0, 1.0, &[], 1e-12, 1e-12, 1000)
    );
    // Even a zero-length interval must have valid tolerances.
    assert_eq!(
        integrate(|_| panic!("validation order"), 0.0, 0.0, 0.0, 0.0, 0),
        Err(Error::InvalidInput)
    );
}

#[test]
fn exact_budget_boundaries_for_every_public_method() {
    for (index, run) in runners().into_iter().enumerate() {
        let function = |x: f64| (-x * x).exp();
        let mut calls = 0;
        let reference = run(
            &mut |x| {
                calls += 1;
                function(x)
            },
            1e-10,
            1e-10,
            50_000,
        )
        .unwrap();
        assert!(calls > 1);
        let required = calls;
        for limit in [0, 1, required - 1, required] {
            calls = 0;
            let result = run(
                &mut |x| {
                    calls += 1;
                    assert!(x.is_finite());
                    function(x)
                },
                1e-10,
                1e-10,
                limit,
            );
            assert!(calls <= limit, "method {index} cap {limit}: {calls}");
            if limit == required {
                assert_eq!(result, Ok(reference));
                assert_eq!(calls, required);
            } else {
                assert_eq!(
                    result,
                    Err(Error::EvaluationLimit),
                    "method {index} cap {limit}"
                );
                assert_eq!(calls, limit);
            }
        }
    }
    // One initial finite rule, three breakpoint segments, and both infinite
    // initial tails all count their full 15-sample rules.
    let mut calls = 0;
    assert!(integrate(
        |_| {
            calls += 1;
            1.0
        },
        0.0,
        1.0,
        1e-12,
        0.0,
        15
    )
    .is_ok());
    assert_eq!(calls, 15);
    calls = 0;
    assert!(integrate_with_breakpoints(
        |_| {
            calls += 1;
            1.0
        },
        0.0,
        1.0,
        &[0.1, 0.9],
        1e-12,
        0.0,
        45
    )
    .is_ok());
    assert_eq!(calls, 45);
    calls = 0;
    assert_eq!(
        integrate_infinite(
            |x| {
                calls += 1;
                (-x * x).exp()
            },
            1e-12,
            0.0,
            29
        ),
        Err(Error::EvaluationLimit)
    );
    assert_eq!(calls, 29);
}

#[test]
fn nan_and_both_infinities_are_never_discarded() {
    for run in runners() {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for fail_at in [1, 7, 17] {
                let mut calls = 0;
                let result = run(
                    &mut |x| {
                        calls += 1;
                        if calls == fail_at {
                            invalid
                        } else {
                            (-x * x).exp()
                        }
                    },
                    1e-12,
                    1e-12,
                    20_000,
                );
                assert_eq!(result, Err(Error::NonFiniteIntegrand));
                assert_eq!(calls, fail_at);
            }
        }
    }
}

#[test]
fn overflow_roundoff_subdivision_and_exhaustion() {
    // Regression: every weighted sample may underflow even though the integral
    // is representable. A relative epsilon floor alone would certify zero.
    for run in [integrate as fn(_, _, _, _, _, _) -> _, integrate_singular] {
        assert!(run(
            |_| f64::from_bits(2),
            0.0,
            1.0,
            f64::from_bits(1),
            0.0,
            10_000
        )
        .is_err());
    }
    assert_eq!(
        integrate(|_| f64::MAX, -1.0, 1.0, 1e-10, 0.0, 1000),
        Err(Error::NumericalBreakdown)
    );
    assert_eq!(
        integrate(|_| 1.0, 0.0, 1.0, 1e-30, 0.0, 1000),
        Err(Error::Roundoff)
    );
    assert_eq!(
        integrate(|x| x, -1.0, 1.0, 0.0, 1e-12, 1000),
        Err(Error::Roundoff)
    );
    // Success is checked before asking an already adequate interval to split.
    let small_end = f64::from_bits(1.0_f64.to_bits() + 256);
    assert!(integrate(|_| 1.0, 1.0, small_end, 1e-12, 0.0, 15).is_ok());
    assert_eq!(
        integrate(
            |_| panic!("collapsed nodes"),
            1.0,
            f64::from_bits(1.0_f64.to_bits() + 1),
            1e-10,
            0.0,
            1000
        ),
        Err(Error::Subdivision)
    );
    // The initial rule fits, but required child nodes eventually collapse.
    let end = f64::from_bits(1.0_f64.to_bits() + 512);
    assert_eq!(
        integrate(
            |x| if x.to_bits() % 2 == 0 { 1.0 } else { 2.0 },
            1.0,
            end,
            1e-24,
            0.0,
            10_000
        ),
        Err(Error::Subdivision)
    );
    assert_eq!(
        integrate(|x| (100.0 * x).sin(), 0.0, 1.0, 1e-12, 0.0, 44),
        Err(Error::EvaluationLimit)
    );
    assert_eq!(
        integrate_upper_infinite(|_| panic!("transform collapsed"), f64::MAX, 1e-10, 0.0, 100),
        Err(Error::Subdivision)
    );
    for error in [
        Error::InvalidInput,
        Error::NonFiniteIntegrand,
        Error::EvaluationLimit,
        Error::Roundoff,
        Error::Subdivision,
        Error::NumericalBreakdown,
    ] {
        assert!(!error.to_string().is_empty());
        let as_error: &dyn std::error::Error = &error;
        assert!(as_error.source().is_none());
    }
}

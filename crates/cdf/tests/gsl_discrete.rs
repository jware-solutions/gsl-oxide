//! The discrete PDF-sum, complementary-tail, and range checks from GSL 2.8 cdf/test.c.
mod support;
use gsl_oxide_cdf::*;
use gsl_oxide_randist as randist;
use support::{assert_relative, TEST_TOL2};

fn check_discrete(
    maximum: u32,
    pdf: impl Fn(u32) -> f64,
    cdf: impl Fn(u32) -> f64,
    survival: impl Fn(u32) -> f64,
) {
    let mut sum = 0.0;
    for k in 0..=maximum {
        sum += pdf(k);
        assert_relative(cdf(k), sum, TEST_TOL2, &format!("CDF sum at k={k}"));
    }
    let mut sum = survival(maximum);
    for k in (0..=maximum).rev() {
        assert_relative(
            survival(k),
            sum,
            TEST_TOL2,
            &format!("survival sum at k={k}"),
        );
        sum += pdf(k);
    }
    check_complements(0, maximum, cdf, survival);
}

fn check_complements(
    minimum: u32,
    maximum: u32,
    cdf: impl Fn(u32) -> f64,
    survival: impl Fn(u32) -> f64,
) {
    for k in minimum..=maximum {
        let p = cdf(k);
        let q = survival(k);
        assert_relative(p + q, 1.0, f64::EPSILON, &format!("P+Q at k={k}"));
        assert!((0.0..=1.0).contains(&p), "P at k={k}: {p}");
        assert!((0.0..=1.0).contains(&q), "Q at k={k}: {q}");
    }
}

#[test]
fn binomial() {
    check_discrete(
        5,
        |k| randist::binomial_pdf(k, 0.3, 5).unwrap(),
        |k| binomial_cdf(k, 0.3, 5).unwrap(),
        |k| binomial_survival(k, 0.3, 5).unwrap(),
    );
}
#[test]
fn poisson() {
    check_discrete(
        25,
        |k| randist::poisson_pdf(k, 2.3).unwrap(),
        |k| poisson_cdf(k, 2.3).unwrap(),
        |k| poisson_survival(k, 2.3).unwrap(),
    );
}
#[test]
fn geometric() {
    check_discrete(
        25,
        |k| randist::geometric_pdf(k, 0.3).unwrap(),
        |k| geometric_cdf(k, 0.3).unwrap(),
        |k| geometric_survival(k, 0.3).unwrap(),
    );
}
#[test]
fn negative_binomial() {
    check_discrete(
        15,
        |k| randist::negative_binomial_pdf(k, 0.3, 5.3).unwrap(),
        |k| negative_binomial_cdf(k, 0.3, 5.3).unwrap(),
        |k| negative_binomial_survival(k, 0.3, 5.3).unwrap(),
    );
}
#[test]
fn pascal() {
    check_discrete(
        15,
        |k| randist::pascal_pdf(k, 0.3, 5).unwrap(),
        |k| pascal_cdf(k, 0.3, 5).unwrap(),
        |k| pascal_survival(k, 0.3, 5).unwrap(),
    );
}
#[test]
fn hypergeometric() {
    check_discrete(
        26,
        |k| randist::hypergeometric_pdf(k, 7, 19, 13).unwrap(),
        |k| hypergeometric_cdf(k, 7, 19, 13).unwrap(),
        |k| hypergeometric_survival(k, 7, 19, 13).unwrap(),
    );
}
#[test]
fn hypergeometric_large_populations() {
    // These two restricted ranges are enabled in GSL; its optional full
    // 13,250,475-count scan is disabled by default and is not duplicated here.
    for (minimum, maximum) in [(3500, 3600), (13_247_474, 13_250_474)] {
        check_complements(
            minimum,
            maximum,
            |k| hypergeometric_cdf(k, 76200, 13174274, 678090).unwrap(),
            |k| hypergeometric_survival(k, 76200, 13174274, 678090).unwrap(),
        );
    }
}

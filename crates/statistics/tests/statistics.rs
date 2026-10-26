//! Independent elementary anchors from centered-sum definitions in the GSL 2.8
//! statistics manual, not from executing GSL or another numerical library.
mod support;

use gsl_oxide_statistics::*;
use support::{close, positive_zero};

type Query = fn(&[f64]) -> Result<f64, Error>;
const QUERIES: [Query; 5] = [
    mean,
    population_variance,
    sample_variance,
    population_standard_deviation,
    sample_standard_deviation,
];

#[test]
fn elementary_anchors() {
    // Sum of centered squares: 10, 32, and 5, respectively.
    let b = 2.0_f64.powi(40);
    for (data, m, pop, sample) in [
        (&[1., 2., 3., 4., 5.][..], 3., 2., 2.5),
        (&[2., 4., 4., 4., 5., 5., 7., 9.][..], 5., 4., 32. / 7.),
        (&[b, b + 1., b + 2., b + 3.][..], b + 1.5, 1.25, 5. / 3.),
        (&[-3., -1., 1., 3.][..], 0., 5., 20. / 3.),
    ] {
        close(mean(data).unwrap(), m, 0., 2e-16);
        close(population_variance(data).unwrap(), pop, 0., 4e-16);
        close(sample_variance(data).unwrap(), sample, 0., 4e-16);
        close(
            population_standard_deviation(data).unwrap(),
            pop.sqrt(),
            0.,
            4e-16,
        );
        close(
            sample_standard_deviation(data).unwrap(),
            sample.sqrt(),
            0.,
            4e-16,
        );
    }
}

#[test]
fn empty_singleton_constants_and_signed_zero() {
    for query in QUERIES {
        assert_eq!(query(&[]), Err(Error::EmptyInput));
    }
    for x in [0., -0., -3., f64::MAX, -f64::MAX, f64::from_bits(1)] {
        assert_eq!(mean(&[x]).unwrap().to_bits(), x.to_bits());
        positive_zero(population_variance(&[x]).unwrap());
        positive_zero(population_standard_deviation(&[x]).unwrap());
        assert_eq!(sample_variance(&[x]), Err(Error::InsufficientData));
        assert_eq!(
            sample_standard_deviation(&[x]),
            Err(Error::InsufficientData)
        );
        for query in &QUERIES[1..] {
            positive_zero(query(&[x; 7]).unwrap());
        }
    }
    for query in &QUERIES[1..] {
        positive_zero(query(&[-0., 0., -0., 0.]).unwrap());
    }
}

#[test]
fn nonfinite_and_overflow() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for query in QUERIES {
            assert_eq!(query(&[1., bad, 2.]), Err(Error::NonFiniteData));
        }
    }
    // Raw sums/squares overflow, but the means and constant variances do not.
    assert_eq!(mean(&[f64::MAX; 3]), Ok(f64::MAX));
    close(mean(&[1e308, 1e308, 0., 0.]).unwrap(), 5e307, 0., 2e-16);
    assert_eq!(mean(&[-f64::MAX, f64::MAX]), Ok(0.));
    for query in &QUERIES[1..] {
        assert_eq!(query(&[-f64::MAX, f64::MAX]), Err(Error::NumericalRange));
    }
    // Batch divides scaled M2 before restoring units: M2=2e308 overflows,
    // population variance=1e308 fits; sample variance and thus deviation fail.
    close(
        population_variance(&[-1e154, 1e154]).unwrap(),
        1e308,
        0.,
        4e-16,
    );
    assert_eq!(
        sample_variance(&[-1e154, 1e154]),
        Err(Error::NumericalRange)
    );
    assert_eq!(
        sample_standard_deviation(&[-1e154, 1e154]),
        Err(Error::NumericalRange)
    );
}

#[test]
fn cancellation_extreme_and_tiny_scales() {
    close(mean(&[1e16, 1., -1e16]).unwrap(), 1. / 3., 0., 2e-16);
    close(mean(&[1e308, 1., -1e308]).unwrap(), 1. / 3., 0., 2e-15);
    close(
        mean(&[1e308, 1e-300, -1e308]).unwrap(),
        1e-300 / 3.,
        0.,
        2e-16,
    );
    close(
        mean(&[f64::MAX, f64::MAX, -f64::MAX, -f64::MAX, 1e-300]).unwrap(),
        1e-300 / 5.,
        0.,
        3e-16,
    );
    let tiny = f64::from_bits(1);
    assert_eq!(mean(&[tiny, 3. * tiny]), Ok(2. * tiny));
    positive_zero(population_variance(&[-tiny, tiny]).unwrap());
    // Squared subnormal *results* are retained; no near-constant threshold.
    let s = 2.0_f64.powi(-530);
    close(
        population_variance(&[-s, s]).unwrap(),
        2.0_f64.powi(-1060),
        0.,
        0.,
    );
    let data = [1e150, -1e150, 1., -1., 1e-150, -1e-150];
    positive_zero(mean(&data).unwrap());
    close(population_variance(&data).unwrap(), 1e300 / 3., 0., 5e-16);
}

#[test]
fn translation_and_scaling() {
    let base = [-4., -1., 0., 2., 8.];
    let v = sample_variance(&base).unwrap();
    for scale in [-16., 0.125, 1e100, 1e-100] {
        let data = base.map(|x| x * scale);
        close(
            sample_variance(&data).unwrap(),
            v * scale * scale,
            0.,
            7e-16,
        );
    }
    let shifted = base.map(|x| x + 2.0_f64.powi(40));
    close(sample_variance(&shifted).unwrap(), v, 0., 4e-16);
}

#[test]
fn paired_anchors_and_contracts() {
    let x = [1., 2., 3.];
    for (y, cov, r) in [
        ([2., 4., 6.], 2., 1.),
        ([6., 4., 2.], -2., -1.),
        ([1., 0., 1.], 0., 0.),
    ] {
        close(sample_covariance(&x, &y).unwrap(), cov, 0., 3e-16);
        close(pearson_correlation(&x, &y).unwrap(), r, 0., 3e-16);
        assert_eq!(sample_covariance(&x, &y), sample_covariance(&y, &x));
    }
    assert_eq!(sample_covariance(&x, &[5.; 3]), Ok(0.));
    assert_eq!(sample_covariance(&[f64::MAX; 3], &x), Ok(0.));
    assert_eq!(
        pearson_correlation(&[f64::MAX; 3], &x),
        Err(Error::ZeroVariance)
    );
    assert_eq!(pearson_correlation(&x, &[5.; 3]), Err(Error::ZeroVariance));
    assert_eq!(pearson_correlation(&[5.; 3], &x), Err(Error::ZeroVariance));
    for query in [sample_covariance, pearson_correlation] {
        assert_eq!(query(&x, &[1., 2.]), Err(Error::LengthMismatch));
        assert_eq!(query(&[], &[]), Err(Error::InsufficientData));
        assert_eq!(query(&[1.], &[2.]), Err(Error::InsufficientData));
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(query(&x, &[1., bad, 3.]), Err(Error::NonFiniteData));
            assert_eq!(query(&[1., bad, 3.], &x), Err(Error::NonFiniteData));
        }
    }
}

#[test]
fn paired_scaling_and_normalization_range() {
    // Perfect linear relationships, with moments/product of deviations outside
    // f64's range. Normalized coordinates still give the Pearson coefficient.
    for (a, b) in [
        (1e200, 1e200),
        (1e-200, 1e-200),
        (1e200, 1e-200),
        (1e-160, 1e-160),
    ] {
        let x = [-a, 0., a];
        let y = [-b, 0., b];
        close(pearson_correlation(&x, &y).unwrap(), 1., 0., 4e-16);
        close(
            pearson_correlation(&x, &y.map(|v| -v)).unwrap(),
            -1.,
            0.,
            4e-16,
        );
    }
    close(
        pearson_correlation(&[-f64::MAX, f64::MAX], &[-1., 1.]).unwrap(),
        1.,
        0.,
        4e-16,
    );
    // Opposite-extreme subtraction from a rounded nonzero mean overflows.
    close(
        pearson_correlation(&[-f64::MAX, f64::MAX, f64::MAX], &[-1., 1., 1.]).unwrap(),
        1.,
        0.,
        5e-16,
    );
    let x = [-1e200, 0., 1e200];
    let y = [-1e-200, 0., 1e-200];
    close(sample_covariance(&x, &y).unwrap(), 1., 0., 4e-16);
    assert_eq!(sample_covariance(&x, &x), Err(Error::NumericalRange));
    let base = [-4., -1., 0., 2., 8.];
    let paired = [7., 2., 9., -3., 4.];
    let cov = sample_covariance(&base, &paired).unwrap();
    let r = pearson_correlation(&base, &paired).unwrap();
    let shifted = base.map(|v| 2.0_f64.powi(40) - 2. * v);
    let scaled = paired.map(|v| 4. * v + 7.);
    close(
        sample_covariance(&shifted, &scaled).unwrap(),
        -8. * cov,
        0.,
        8e-16,
    );
    close(
        pearson_correlation(&shifted, &scaled).unwrap(),
        -r,
        0.,
        8e-16,
    );
    let t = f64::from_bits(1);
    close(
        pearson_correlation(&[-t, t], &[-t, t]).unwrap(),
        1.,
        0.,
        4e-16,
    );
}

#[test]
fn finite_comparison_helper_rejects_invalid_arguments() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.] {
        assert!(std::panic::catch_unwind(|| close(0., 0., bad, 0.)).is_err());
        assert!(std::panic::catch_unwind(|| close(0., 0., 0., bad)).is_err());
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(std::panic::catch_unwind(|| close(bad, 0., 0., 0.)).is_err());
        assert!(std::panic::catch_unwind(|| close(0., bad, 0., 0.)).is_err());
    }
    close(1.01, 1., 0.02, 0.);
    close(1.01, 1., 0., 0.02);
    assert!(std::panic::catch_unwind(|| close(1.1, 1., 0.01, 0.01)).is_err());
}

#[test]
fn errors_have_display_and_standard_error() {
    for error in [
        Error::EmptyInput,
        Error::InsufficientData,
        Error::LengthMismatch,
        Error::NonFiniteData,
        Error::InvalidWeight,
        Error::ZeroTotalWeight,
        Error::WeightRange,
        Error::ZeroVariance,
        Error::NumericalRange,
        Error::CountOverflow,
    ] {
        let standard: &dyn std::error::Error = &error;
        assert!(!standard.to_string().is_empty());
    }
}

//! Estimators follow GSL 2.8 statistics, "Weighted Samples". All anchors below
//! are independently calculated from elementary centered sums.
mod support;

use gsl_oxide_statistics::*;
use support::{close, positive_zero};

type Query = fn(&[f64], &[f64]) -> Result<f64, Error>;
const QUERIES: [Query; 3] = [
    weighted_mean,
    weighted_population_variance,
    weighted_reliability_variance,
];

/// Checks independent weighted anchors.
#[test]
fn independent_weighted_anchors() {
    // W=4, W2=10, mean=13/4, M2_w=27/4. Reliability denominator=3/2,
    // whereas the deliberately unsupported frequency denominator would be 3.
    for (query, expected) in QUERIES.into_iter().zip([13. / 4., 27. / 16., 9. / 2.]) {
        close(query(&[1., 4.], &[1., 3.]).unwrap(), expected, 0., 4e-16);
    }
    // W=6, mean=1, M2=30, W2=14: variance=5, reliability=90/11.
    close(
        weighted_mean(&[-3., 0., 3.], &[1., 2., 3.]).unwrap(),
        1.,
        0.,
        4e-16,
    );
    close(
        weighted_population_variance(&[-3., 0., 3.], &[1., 2., 3.]).unwrap(),
        5.,
        0.,
        4e-16,
    );
    close(
        weighted_reliability_variance(&[-3., 0., 3.], &[1., 2., 3.]).unwrap(),
        90. / 11.,
        0.,
        4e-16,
    );
}

/// Checks degenerate weights and removal before value validation.
#[test]
fn degenerate_weights_and_removal_before_value_validation() {
    for query in QUERIES {
        assert_eq!(query(&[], &[]), Err(Error::ZeroTotalWeight));
        assert_eq!(
            query(&[f64::NAN, f64::INFINITY], &[0., -0.]),
            Err(Error::ZeroTotalWeight)
        );
        assert_eq!(query(&[1., 2.], &[0., -0.]), Err(Error::ZeroTotalWeight));
        assert_eq!(query(&[1., 2.], &[1.]), Err(Error::LengthMismatch));
        for bad in [-1., f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(query(&[1., 2.], &[1., bad]), Err(Error::InvalidWeight));
        }
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(query(&[1., bad], &[1., 1.]), Err(Error::NonFiniteData));
            let expected = query(&[1., 4.], &[1., 3.]);
            assert_eq!(query(&[1., bad, 4.], &[1., -0., 3.]), expected);
        }
    }
    assert_eq!(weighted_mean(&[f64::NAN, -2.], &[0., 3.]), Ok(-2.));
    positive_zero(weighted_population_variance(&[f64::INFINITY, -2.], &[-0., 3.]).unwrap());
    assert_eq!(
        weighted_reliability_variance(&[f64::NAN, -2.], &[0., 3.]),
        Err(Error::InsufficientData)
    );
    for x in [0., -0., f64::MAX, -f64::MAX, f64::from_bits(1)] {
        assert_eq!(weighted_mean(&[x], &[2.]).unwrap().to_bits(), x.to_bits());
        for query in &QUERIES[1..] {
            positive_zero(query(&[x; 3], &[1., 2., 4.]).unwrap());
        }
    }
}

/// Checks equal weights and common rescaling.
#[test]
fn equal_weights_and_common_rescaling() {
    let b = 2.0_f64.powi(40);
    for data in [
        [1., 2., 3., 4.],
        [b, b + 1., b + 2., b + 3.],
        [-7., -2., 3., 5.],
    ] {
        for w in [1., 1e300, 1e-300, f64::MAX] {
            for (query, expected) in QUERIES.into_iter().zip([
                mean(&data).unwrap(),
                population_variance(&data).unwrap(),
                sample_variance(&data).unwrap(),
            ]) {
                close(query(&data, &[w; 4]).unwrap(), expected, 0., 5e-16);
            }
        }
    }
    let data = [-7., -2., 3., 5.];
    let weights = [1., 2., 3., 7.];
    for factor in [1e-300, 1e300, 0.125, 16.] {
        for query in QUERIES {
            close(
                query(&data, &weights.map(|w| w * factor)).unwrap(),
                query(&data, &weights).unwrap(),
                0.,
                8e-16,
            );
        }
    }
}

/// Checks strong imbalance and weight range.
#[test]
fn strong_imbalance_and_weight_range() {
    // For two observations a,b and positive weights u,v, the reliability
    // variance is (a-b)^2/2 regardless of the weight ratio (algebraic anchor).
    // A literal W-W2/W rounds to zero for these dominant weights.
    for tiny in [1e-20, 1e-100, 1e-300, f64::MIN_POSITIVE] {
        for weights in [[1., tiny], [tiny, 1.]] {
            close(
                weighted_reliability_variance(&[1., 4.], &weights).unwrap(),
                4.5,
                0.,
                7e-16,
            );
        }
        close(
            weighted_mean(&[0., 1.], &[1., tiny]).unwrap(),
            tiny / (1. + tiny),
            0.,
            4e-16,
        );
        close(
            weighted_population_variance(&[0., 1.], &[1., tiny]).unwrap(),
            tiny / ((1. + tiny) * (1. + tiny)),
            0.,
            7e-16,
        );
    }
    for weights in [[f64::MAX, f64::from_bits(1)], [1., 1e-310]] {
        for query in QUERIES {
            assert_eq!(query(&[1., 4.], &weights), Err(Error::WeightRange));
        }
    }
    // Both raw W and W2 overflow; normalization retains ordinary estimators.
    for query in QUERIES {
        close(
            query(&[1., 4.], &[1e308, 1e308]).unwrap(),
            query(&[1., 4.], &[1., 1.]).unwrap(),
            0.,
            0.,
        );
    }
}

/// Checks weighted value range and identities.
#[test]
fn weighted_value_range_and_identities() {
    close(
        weighted_mean(&[1e308, 1e-300, -1e308], &[1.; 3]).unwrap(),
        1e-300 / 3.,
        0.,
        2e-16,
    );
    close(
        weighted_mean(
            &[f64::MAX, f64::MAX, -f64::MAX, -f64::MAX, 1e-300],
            &[1.; 5],
        )
        .unwrap(),
        1e-300 / 5.,
        0.,
        3e-16,
    );
    close(
        weighted_mean(&[1e308, 1e308, -1e308], &[1e308; 3]).unwrap(),
        1e308 / 3.,
        0.,
        4e-16,
    );
    let data = [-4., -1., 3., 8.];
    let weights = [2., 3., 7., 1.];
    // Exact 2^-1060 and 2^-1059 in units of the smallest subnormal, 2^-1074.
    // powi with these negative exponents can underflow to zero on some platforms.
    for (query, tiny_expected) in QUERIES[1..]
        .iter()
        .zip([f64::from_bits(1 << 14), f64::from_bits(1 << 15)])
    {
        let expected = query(&data, &weights).unwrap();
        close(
            query(&data.map(|x| x + 2.0_f64.powi(40)), &weights).unwrap(),
            expected,
            0.,
            6e-16,
        );
        close(
            query(&data.map(|x| -16. * x), &weights).unwrap(),
            expected * 256.,
            0.,
            6e-16,
        );
        assert_eq!(
            query(&[-f64::MAX, f64::MAX], &[1., 1.]),
            Err(Error::NumericalRange)
        );
        let tiny = 2.0_f64.powi(-530);
        close(
            query(&[-tiny, tiny], &[1., 1.]).unwrap(),
            tiny_expected,
            0.,
            0.,
        );
    }
}

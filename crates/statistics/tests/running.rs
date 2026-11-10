//! Independent centered-sum anchors from GSL statistics/rstat definitions.
//! Merge checks also exercise Chan/Pébay's count-weighted pairwise identity
//! (Pébay, SAND2008-6212, §2), without an external oracle.
mod support;

use gsl_oxide_statistics::*;
use support::{close, positive_zero};

/// Builds a running statistics accumulator by pushing all observations in `data`.
fn accumulate(data: &[f64]) -> RunningStatistics {
    let mut result = RunningStatistics::new();
    for &value in data {
        result.push(value).unwrap();
    }
    result
}

/// Checks accumulator `s` against expected count, mean, and variances using the provided tolerances.
fn anchor(
    s: &RunningStatistics,
    count: usize,
    m: f64,
    pop: f64,
    sample: f64,
    abs_mean: f64,
    rel: f64,
) {
    assert_eq!(s.count(), count);
    // Querying repeatedly must not change state.
    for _ in 0..3 {
        close(s.mean().unwrap(), m, abs_mean, 3e-16);
        close(s.population_variance().unwrap(), pop, 0., rel);
        close(s.sample_variance().unwrap(), sample, 0., rel);
        close(
            s.population_standard_deviation().unwrap(),
            pop.sqrt(),
            0.,
            rel,
        );
        close(
            s.sample_standard_deviation().unwrap(),
            sample.sqrt(),
            0.,
            rel,
        );
    }
}

/// Compares accumulator `s` against batch statistics of `data` with relative tolerance `rel`.
fn compare(s: &RunningStatistics, data: &[f64], rel: f64) {
    anchor(
        s,
        data.len(),
        mean(data).unwrap(),
        population_variance(data).unwrap(),
        sample_variance(data).unwrap(),
        1e-15,
        rel,
    );
}

/// Checks independent anchors and deterministic orders.
#[test]
fn independent_anchors_and_deterministic_orders() {
    let b = 2.0_f64.powi(40);
    let cases: [(&[f64], f64, f64, f64); 4] = [
        (&[1., 2., 3., 4., 5.], 3., 2., 2.5),
        (&[2., 4., 4., 4., 5., 5., 7., 9.], 5., 4., 32. / 7.),
        (&[b, b + 1., b + 2., b + 3.], b + 1.5, 1.25, 5. / 3.),
        (&[-3., -1., 1., 3.], 0., 5., 20. / 3.),
    ];
    for (data, m, pop, sample) in cases {
        let orders = [
            data.to_vec(),
            data.iter().rev().copied().collect(),
            data.iter()
                .step_by(2)
                .chain(data.iter().skip(1).step_by(2))
                .copied()
                .collect(),
        ];
        for order in orders {
            anchor(
                &accumulate(&order),
                data.len(),
                m,
                pop,
                sample,
                1e-15,
                1e-15,
            );
        }
    }
}

/// Checks every split including empty and unequal chunks.
#[test]
fn every_split_including_empty_and_unequal_chunks() {
    let b = 2.0_f64.powi(40);
    for data in [
        vec![1., 2., 3., 4., 5.],
        vec![2., 4., 4., 4., 5., 5., 7., 9.],
        vec![b, b + 1., b + 2., b + 3., b + 4.],
        vec![-7., 19., 2., -3., 0., 14., -11.],
    ] {
        for split in 0..=data.len() {
            let mut left = accumulate(&data[..split]);
            let right = accumulate(&data[split..]);
            let before = format!("{right:?}");
            left.merge(&right).unwrap();
            assert_eq!(format!("{right:?}"), before);
            compare(&left, &data, 2e-15);
            left.push(data[0]).unwrap();
            let mut extended = data.clone();
            extended.push(data[0]);
            compare(&left, &extended, 2e-15);
        }
    }
    // Two chunk means are 0 and 10, with counts 1 and 9: merged mean is 9,
    // M2=90, population variance=9, sample variance=10.
    let mut singleton = accumulate(&[0.]);
    singleton.merge(&accumulate(&[10.; 9])).unwrap();
    anchor(&singleton, 10, 9., 9., 10., 0., 4e-16);
}

/// Checks singleton chunks sequential and balanced trees.
#[test]
fn singleton_chunks_sequential_and_balanced_trees() {
    let data = [-8., -3., 0., 1., 2., 5., 6., 13.];
    let mut sequential = RunningStatistics::new();
    for value in data {
        sequential.merge(&RunningStatistics::new()).unwrap();
        sequential.merge(&accumulate(&[value])).unwrap();
    }
    compare(&sequential, &data, 2e-15);
    let mut chunks: Vec<_> = data.iter().map(|&v| accumulate(&[v])).collect();
    while chunks.len() > 1 {
        chunks = chunks
            .chunks(2)
            .map(|pair| {
                let mut left = pair[0].clone();
                if pair.len() == 2 {
                    left.merge(&pair[1]).unwrap();
                }
                left
            })
            .collect();
    }
    compare(&chunks[0], &data, 2e-15);
    let mut uneven = accumulate(&data[..1]);
    uneven.merge(&accumulate(&data[1..6])).unwrap();
    uneven.merge(&accumulate(&data[6..])).unwrap();
    compare(&uneven, &data, 2e-15);
}

/// Checks difficult inputs and sub ulp mean residual.
#[test]
fn difficult_inputs_and_sub_ulp_mean_residual() {
    // Compensation must preserve delta's low part as well as mean increments.
    let data = [1e16, 1., -1e16];
    close(accumulate(&data).mean().unwrap(), 1. / 3., 0., 4e-16);
    compare(&accumulate(&data), &data, 6e-16);
    for split in 0..=data.len() {
        let mut left = accumulate(&data[..split]);
        left.merge(&accumulate(&data[split..])).unwrap();
        compare(&left, &data, 6e-16);
    }
    let b = 2.0_f64.powi(52);
    // True mean B+1/2 rounds to B, but M2 remains 1/2.
    let pair = [b, b + 1.];
    anchor(&accumulate(&pair), 2, b, 0.25, 0.5, 0., 0.);
    close(population_variance(&pair).unwrap(), 0.25, 0., 0.);
    let mut merged = accumulate(&pair[..1]);
    merged.merge(&accumulate(&pair[1..])).unwrap();
    anchor(&merged, 2, b, 0.25, 0.5, 0., 0.);
    for scale in [1e150, 1e-150, 2.0_f64.powi(-530)] {
        let data = [-scale, 0., scale];
        // Finite tiny M2 may suffer a few subnormal ulps in online updates.
        let s = accumulate(&data);
        close(s.mean().unwrap(), 0., scale.abs() * 2e-16, 0.);
        close(
            s.sample_variance().unwrap(),
            scale * scale,
            f64::from_bits(4),
            8e-16,
        );
        close(
            s.population_variance().unwrap(),
            (2. / 3.) * scale * scale,
            f64::from_bits(4),
            8e-16,
        );
    }
    let data = [1e150, -1e150, 1., -1., 1e-150, -1e-150];
    let s = accumulate(&data);
    close(s.mean().unwrap(), 0., 1e-15, 0.);
    close(s.population_variance().unwrap(), 1e300 / 3., 0., 7e-16);
}

/// Checks degenerate inputs reset and reuse.
#[test]
fn degenerate_inputs_reset_and_reuse() {
    let mut s = RunningStatistics::default();
    assert_eq!(s.count(), 0);
    assert_eq!(s.mean(), Err(Error::EmptyInput));
    assert_eq!(s.population_variance(), Err(Error::EmptyInput));
    assert_eq!(s.sample_variance(), Err(Error::EmptyInput));
    assert_eq!(s.population_standard_deviation(), Err(Error::EmptyInput));
    assert_eq!(s.sample_standard_deviation(), Err(Error::EmptyInput));
    s.merge(&RunningStatistics::new()).unwrap();
    assert_eq!(s.count(), 0);
    for value in [0., -0., -2., f64::MAX, -f64::MAX, f64::from_bits(1)] {
        s.push(value).unwrap();
        assert_eq!(s.mean().unwrap().to_bits(), value.to_bits());
        positive_zero(s.population_variance().unwrap());
        positive_zero(s.population_standard_deviation().unwrap());
        assert_eq!(s.sample_variance(), Err(Error::InsufficientData));
        assert_eq!(s.sample_standard_deviation(), Err(Error::InsufficientData));
        s.push(value).unwrap();
        positive_zero(s.sample_variance().unwrap());
        positive_zero(s.sample_standard_deviation().unwrap());
        s.merge(&accumulate(&[value; 3])).unwrap();
        positive_zero(s.sample_variance().unwrap());
        assert_eq!(s.count(), 5);
        s.reset();
        assert_eq!(s.count(), 0);
    }
    for x in [1., 2., 3., 4., 5.] {
        s.push(x).unwrap();
    }
    anchor(&s, 5, 3., 2., 2.5, 0., 4e-16);
}

/// Checks rejected pushes and merges preserve state and allow continuation.
#[test]
fn rejected_pushes_and_merges_preserve_state_and_allow_continuation() {
    for initial in [vec![], vec![1., 2., 3.]] {
        let mut s = accumulate(&initial);
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let before = format!("{s:?}");
            assert_eq!(s.push(bad), Err(Error::NonFiniteData));
            assert_eq!(format!("{s:?}"), before);
        }
        s.push(4.).unwrap();
        assert_eq!(s.count(), initial.len() + 1);
    }
    for (initial, value) in [(0., 1e308), (f64::MAX, -f64::MAX), (-1e154, 1e154)] {
        let mut s = accumulate(&[initial]);
        let before = format!("{s:?}");
        assert_eq!(s.push(value), Err(Error::NumericalRange));
        assert_eq!(format!("{s:?}"), before);
        let other = accumulate(&[value]);
        let other_before = format!("{other:?}");
        assert_eq!(s.merge(&other), Err(Error::NumericalRange));
        assert_eq!(format!("{s:?}"), before);
        assert_eq!(format!("{other:?}"), other_before);
        s.push(initial).unwrap();
        assert_eq!(s.count(), 2);
        positive_zero(s.sample_variance().unwrap());
    }
}

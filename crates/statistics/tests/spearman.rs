//! GSL 2.8 Spearman regression cases and independent rank-correlation contracts.
// The two GSL regression cases below are adapted from GSL 2.8
// statistics/test_float_source.c and statistics/test_int_source.c.
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jim Davies, Brian Gough
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
// Rust adaptation: Copyright (c) 2026 JWare Solutions.
// 2026-10-25: Ported Spearman inputs, expected values, and relative tolerances
// to contiguous f64 slices and Rust assertions; added independent contract tests.
mod support;

use gsl_oxide_statistics::{spearman_correlation, Error};
use support::{close, positive_zero};

/// Checks gsl floating point regression.
#[test]
fn gsl_floating_point_regression() {
    let x = [
        0.0421, 0.0941, 0.1064, 0.0242, 0.1331, 0.0773, 0.0243, 0.0815, 0.1186, 0.0356, 0.0728,
        0.0999, 0.0614, 0.0479,
    ];
    let y = [
        0.1081, 0.0986, 0.1566, 0.1961, 0.1125, 0.1942, 0.1079, 0.1021, 0.1583, 0.1673, 0.1675,
        0.1856, 0.1688, 0.1512,
    ];
    close(
        spearman_correlation(&x, &y).unwrap(),
        -0.160_439_560_439_560_45,
        0.,
        1e-10,
    );
}

/// Checks gsl integer regression with ties.
#[test]
fn gsl_integer_regression_with_ties() {
    let x = [
        17., 18., 16., 18., 12., 20., 18., 20., 20., 22., 20., 10., 8., 12., 16., 16., 18., 20.,
        18., 21.,
    ];
    let y = [
        19., 20., 22., 24., 10., 25., 20., 22., 21., 23., 20., 10., 12., 14., 12., 20., 22., 24.,
        23., 17.,
    ];
    close(
        spearman_correlation(&x, &y).unwrap(),
        0.644239670474577,
        0.,
        1e-10,
    );
}

/// Checks average ranks for ties.
#[test]
fn average_ranks_for_ties() {
    // Hand-computed centered rank sums. The first case has M2_x=1.5,
    // M2_y=2, C=1.5; the next three have M2_x=M2_y=4.5, C=3.75.
    for (x, y, expected) in [
        (&[1., 1., 2.][..], &[1., 2., 3.][..], 3.0_f64.sqrt() / 2.),
        (&[1., 1., 2., 3.][..], &[1., 2., 2., 3.][..], 5. / 6.),
        (&[1., 2., 3., 3.][..], &[1., 2., 2., 3.][..], 5. / 6.),
        (&[-0., 0., 1., 2.][..], &[1., 2., 2., 3.][..], 5. / 6.),
    ] {
        close(spearman_correlation(x, y).unwrap(), expected, 5e-16, 5e-16);
        close(spearman_correlation(y, x).unwrap(), expected, 5e-16, 5e-16);
    }
    positive_zero(spearman_correlation(&[1., 1., 2., 2.], &[1., 2., 1., 2.]).unwrap());
}

/// Checks monotonic association and extreme values.
#[test]
fn monotonic_association_and_extreme_values() {
    let x = [4.0_f64, -2., 8., 0., 1.];
    let cubic = x.map(|v| v.powi(3));
    close(spearman_correlation(&x, &cubic).unwrap(), 1., 0., 5e-16);
    close(
        spearman_correlation(&x, &cubic.map(|v| -v)).unwrap(),
        -1.,
        0.,
        5e-16,
    );
    for y in [[10., 20.], [20., 10.]] {
        let expected = if y[0] < y[1] { 1. } else { -1. };
        close(
            spearman_correlation(&[1., 2.], &y).unwrap(),
            expected,
            0.,
            5e-16,
        );
    }
    let tiny = f64::from_bits(1);
    let extreme = [-f64::MAX, -1., -tiny, -0., 0., tiny, 1., f64::MAX];
    let labels = [-4., -3., -2., 0., -0., 2., 3., 4.];
    close(
        spearman_correlation(&extreme, &labels).unwrap(),
        1.,
        0.,
        5e-16,
    );
    // Adjacent representable values must not be treated as ties.
    let next = f64::from_bits(1.0_f64.to_bits() + 1);
    close(
        spearman_correlation(&[1., next, 2.], &[1., 2., 3.]).unwrap(),
        1.,
        0.,
        5e-16,
    );
}

/// Checks pairing order symmetry and input preservation.
#[test]
fn pairing_order_symmetry_and_input_preservation() {
    let x = [8., -0., 3., 8., 5., 0., 1.];
    let y = [6., 1., 7., 2., 3., 1., 5.];
    let before_x = x.map(f64::to_bits);
    let before_y = y.map(f64::to_bits);
    let expected = spearman_correlation(&x, &y).unwrap();
    close(
        spearman_correlation(&y, &x).unwrap(),
        expected,
        5e-16,
        5e-16,
    );
    for offset in 0..x.len() {
        let mut px = x;
        let mut py = y;
        px.rotate_left(offset);
        py.rotate_left(offset);
        close(
            spearman_correlation(&px, &py).unwrap(),
            expected,
            5e-16,
            5e-16,
        );
    }
    assert_eq!(x.map(f64::to_bits), before_x);
    assert_eq!(y.map(f64::to_bits), before_y);
}

/// Checks untied permutations match rank difference formula.
#[test]
fn untied_permutations_match_rank_difference_formula() {
    // With no ties rho = 1 - 6*sum((rank_x-rank_y)^2)/(n*(n^2-1)).
    // Check every permutation of five ranks against this independent identity.
    /// Recursively checks Spearman correlation for permutations of `y` from index `start`.
    fn check(y: &mut [f64; 5], start: usize) {
        if start == y.len() {
            let x = [1., 2., 3., 4., 5.];
            let squared_differences: f64 =
                x.iter().zip(y.iter()).map(|(a, b)| (a - b).powi(2)).sum();
            let expected = 1. - squared_differences / 20.;
            close(spearman_correlation(&x, y).unwrap(), expected, 5e-16, 5e-16);
        } else {
            for i in start..y.len() {
                y.swap(start, i);
                check(y, start + 1);
                y.swap(start, i);
            }
        }
    }
    check(&mut [1., 2., 3., 4., 5.], 0);
}

/// Checks invalid input contracts.
#[test]
fn invalid_input_contracts() {
    assert_eq!(spearman_correlation(&[], &[1.]), Err(Error::LengthMismatch));
    assert_eq!(
        spearman_correlation(&[1., 2.], &[1.]),
        Err(Error::LengthMismatch)
    );
    assert_eq!(spearman_correlation(&[], &[]), Err(Error::InsufficientData));
    assert_eq!(
        spearman_correlation(&[1.], &[2.]),
        Err(Error::InsufficientData)
    );
    for constant in [[1., 1., 1.], [-0., 0., -0.], [f64::MAX; 3]] {
        assert_eq!(
            spearman_correlation(&constant, &[1., 2., 3.]),
            Err(Error::ZeroVariance)
        );
        assert_eq!(
            spearman_correlation(&[1., 2., 3.], &constant),
            Err(Error::ZeroVariance)
        );
        assert_eq!(
            spearman_correlation(&constant, &constant),
            Err(Error::ZeroVariance)
        );
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            spearman_correlation(&[1., bad, 3.], &[1., 2., 3.]),
            Err(Error::NonFiniteData)
        );
        assert_eq!(
            spearman_correlation(&[1., 2., 3.], &[1., bad, 3.]),
            Err(Error::NonFiniteData)
        );
        assert_eq!(
            spearman_correlation(&[1.; 3], &[1., bad, 3.]),
            Err(Error::NonFiniteData)
        );
    }
}

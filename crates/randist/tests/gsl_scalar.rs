//! Rust ports of the GSL 2.8 random-distribution tests.
// Ported from GSL 2.8 randist/test.c (GPL-3.0-or-later).
// Copyright (C) James Theiler, Brian Gough and GSL contributors.
#![warn(clippy::missing_docs_in_private_items)]
/// Rust-specific parameter, boundary, and numerical-range checks.
mod contracts;
/// Fixed random source and statistical comparison helpers.
mod support;
use gsl_oxide_randist::*;
use support::*;

/// Checks GSL's Dirichlet marginal histogram for concentrations `[2.5, 5.0]`.
#[test]
fn gsl_dirichlet() {
    continuous(
        |r| dirichlet(r, &[2.5, 5.0]).unwrap()[0],
        |x| {
            if x <= 0.0 || x >= 1.0 {
                0.0
            } else {
                dirichlet_density(&[2.5, 5.0], &[x, 1.0 - x]).unwrap()
            }
        },
    );
}

/// Checks the small-concentration Dirichlet histogram against analytic beta bin probabilities.
#[test]
fn gsl_dirichlet_small() {
    // GSL's endpoint singularity has power -0.9975, beyond the integration
    // crate's supported powers. Integrate its beta marginal analytically.
    let mut p = [0.0; 100];
    for (i, p) in p.iter_mut().enumerate() {
        let a = (-5.0 + i as f64 * 0.1).clamp(0.0, 1.0);
        let b = (-5.0 + (i + 1) as f64 * 0.1).clamp(0.0, 1.0);
        *p = gsl_oxide_special::incomplete_beta::regularized(0.0025, 0.005, b).unwrap()
            - gsl_oxide_special::incomplete_beta::regularized(0.0025, 0.005, a).unwrap();
    }
    continuous_masses(|r| dirichlet(r, &[0.0025, 0.005]).unwrap()[0], p);
    for x in [0.01, 0.3, 0.9] {
        relative(
            dirichlet_density(&[0.0025, 0.005], &[x, 1.0 - x]).unwrap(),
            beta_density(x, 0.0025, 0.005).unwrap(),
            1e-12,
        );
    }
}

/// Checks GSL's centered and scaled gamma histogram at shape `2^32`.
#[test]
fn gsl_gamma_vlarge() {
    continuous(
        |r| (gamma(r, 4294967296.0, 6.32899304917e-10).unwrap() - 2.71828181565) * 1e4,
        |x| gamma_density(x / 1e4 + 2.71828181565, 4294967296.0, 6.32899304917e-10).unwrap() / 1e4,
    );
}
/// Checks the centered and scaled Knuth gamma histogram at shape `2^32`.
#[test]
fn gsl_gamma_knuth_vlarge() {
    continuous(
        |r| (gamma_knuth(r, 4294967296.0, 6.32899304917e-10).unwrap() - 2.71828181565) * 1e4,
        |x| gamma_density(x / 1e4 + 2.71828181565, 4294967296.0, 6.32899304917e-10).unwrap() / 1e4,
    );
}

/// Checks the first bivariate-normal marginal against its normal density.
#[test]
fn gsl_bivariate_gaussian1() {
    continuous(
        |r| bivariate_normal(r, 3.0, 2.0, 0.3).unwrap().0,
        |x| normal_density(x, 3.0).unwrap(),
    );
}
/// Checks the second bivariate-normal marginal using GSL's numerical marginalization.
#[test]
fn gsl_bivariate_gaussian2() {
    continuous(
        |r| bivariate_normal(r, 3.0, 2.0, 0.3).unwrap().1,
        |y| {
            // Preserve GSL's ten-point numerical marginalization regression.
            (0..10)
                .map(|i| {
                    bivariate_normal_density(-10.0 + i as f64 * 2.0, y, 3.0, 2.0, 0.3).unwrap()
                        * 2.0
                })
                .sum()
        },
    );
}
/// Checks the sum of positively correlated normal components against its normal density.
#[test]
fn gsl_bivariate_gaussian3() {
    continuous(
        |r| {
            let (x, y) = bivariate_normal(r, 3.0, 2.0, 0.3).unwrap();
            x + y
        },
        |x| normal_density(x, 16.6_f64.sqrt()).unwrap(),
    );
}
/// Checks the sum of negatively correlated normal components against its normal density.
#[test]
fn gsl_bivariate_gaussian4() {
    continuous(
        |r| {
            let (x, y) = bivariate_normal(r, 3.0, 2.0, -0.5).unwrap();
            x + y
        },
        |x| normal_density(x, 7.0_f64.sqrt()).unwrap(),
    );
}

/// Returns the uniform angle density at `x` in radians.
///
/// The reference value is `1/(2*pi)` on `[-pi,pi]` and zero outside it.
fn angle_density(x: f64) -> f64 {
    if x.abs() <= std::f64::consts::PI {
        0.5 / std::f64::consts::PI
    } else {
        0.0
    }
}
/// Checks that disk-based circle directions have a uniform azimuth.
#[test]
fn gsl_dir2d() {
    continuous(
        |r| {
            let [x, y] = direction_2d(r);
            y.atan2(x)
        },
        angle_density,
    );
}
/// Checks that trigonometric circle directions have a uniform azimuth.
#[test]
fn gsl_dir2d_trig_method() {
    continuous(
        |r| {
            let [x, y] = direction_2d_trigonometric(r);
            y.atan2(x)
        },
        angle_density,
    );
}
/// Checks uniform azimuth of spherical directions in the xy plane.
#[test]
fn gsl_dir3dxy() {
    continuous(
        |r| {
            let [x, y, _] = direction_3d(r);
            y.atan2(x)
        },
        angle_density,
    );
}
/// Checks uniform azimuth of spherical directions in the yz plane.
#[test]
fn gsl_dir3dyz() {
    continuous(
        |r| {
            let [_, y, z] = direction_3d(r);
            z.atan2(y)
        },
        angle_density,
    );
}
/// Checks uniform azimuth of spherical directions in the zx plane.
#[test]
fn gsl_dir3dzx() {
    continuous(
        |r| {
            let [x, _, z] = direction_3d(r);
            x.atan2(z)
        },
        angle_density,
    );
}

/// Checks alias-table frequencies for GSL weights `[0.59, 0.4, 0.01]`.
#[test]
fn gsl_discrete1() {
    let d = DiscreteDistribution::new(&[0.59, 0.4, 0.01]).unwrap();
    discrete(|r| d.sample(r) as u32, |k| d.pdf(k as usize));
}
/// Checks alias-table frequencies for GSL's unnormalized weights, including a zero weight.
#[test]
fn gsl_discrete2() {
    let d = DiscreteDistribution::new(&[1., 9., 3., 4., 5., 8., 6., 7., 2., 0.]).unwrap();
    discrete(|r| d.sample(r) as u32, |k| d.pdf(k as usize));
}
/// Checks alias-table frequencies for twenty equally weighted categories.
#[test]
fn gsl_discrete3() {
    let d = DiscreteDistribution::new(&[0.05; 20]).unwrap();
    discrete(|r| d.sample(r) as u32, |k| d.pdf(k as usize));
}
/// Returns the reference probability of category count `k` in 100 multinomial trials.
///
/// The category probability is 0.2, represented by weights `[0.4,1.6]`.
/// The caller supplies `k <= 100` so the complementary count is nonnegative.
fn multinomial_marginal(k: u32) -> f64 {
    multinomial_pdf(&[0.4, 1.6], &[k, 100 - k]).unwrap()
}
/// Checks the first multinomial marginal for unnormalized weights `[2,7,1]`.
#[test]
fn gsl_multinomial() {
    discrete(
        |r| multinomial(r, 100, &[2., 7., 1.]).unwrap()[0],
        multinomial_marginal,
    );
}
/// Checks the first multinomial marginal for GSL's ten-category weight vector.
#[test]
fn gsl_multinomial_large() {
    discrete(
        |r| {
            multinomial(
                r,
                100,
                &[0.2, 0.2, 0.17, 0.14, 0.12, 0.07, 0.05, 0.04, 0.01, 0.0],
            )
            .unwrap()[0]
        },
        multinomial_marginal,
    );
}

/// Checks GSL's normal, exponential, Cauchy, and categorical interval probabilities.
#[test]
fn gsl_interval_probabilities() {
    interval(standard_normal, 0.0, 100.0, 0.5);
    interval(standard_normal, -1.0, 1.0, 0.6826895);
    interval(standard_normal, 3.0, 3.5, 0.0011172689);
    interval(
        |r| standard_normal_tail(r, 3.0).unwrap(),
        3.0,
        3.5,
        0.0011172689 / 0.0013498981,
    );
    interval(
        |r| exponential(r, 2.0).unwrap(),
        0.0,
        1.0,
        1.0 - (-0.5_f64).exp(),
    );
    interval(|r| cauchy(r, 2.0).unwrap(), 0.0, 10000.0, 0.5);
    let d = DiscreteDistribution::new(&[0.59, 0.4, 0.01]).unwrap();
    for (a, b, p) in [(-0.5, 0.5, 0.59), (0.5, 1.5, 0.4), (1.5, 3.5, 0.01)] {
        interval(|r| d.sample(r) as f64, a, b, p);
    }
    let d = DiscreteDistribution::new(&[1., 9., 3., 4., 5., 8., 6., 7., 2., 0.]).unwrap();
    interval(|r| d.sample(r) as f64, -0.5, 0.5, 1.0 / 45.0);
    interval(|r| d.sample(r) as f64, 8.5, 9.5, 0.0);
    let d = DiscreteDistribution::new(&[0.05; 20]).unwrap();
    for (a, b, p) in [(-0.5, 0.5, 0.05), (0.5, 1.5, 0.05), (-0.5, 9.5, 0.5)] {
        interval(|r| d.sample(r) as f64, a, b, p);
    }
}

/// Checks the ten-coordinate Dirichlet means against their theoretical sampling errors.
#[test]
fn gsl_dirichlet_moments() {
    let mut r = Generator::new();
    let alpha: Vec<_> = (0..10).map(|_| exponential(&mut r, 0.1).unwrap()).collect();
    let total: f64 = alpha.iter().sum();
    let mut sums = [0.0; 10];
    for _ in 0..N {
        let theta = dirichlet(&mut r, &alpha).unwrap();
        assert!((theta.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        for (sum, x) in sums.iter_mut().zip(theta) {
            *sum += x;
        }
    }
    for k in 0..10 {
        let mean = alpha[k] / total;
        let sd = (alpha[k] * (1.0 - mean) / (total * (total + 1.0))).sqrt();
        assert!((N as f64).sqrt() * (sums[k] / N as f64 - mean).abs() / sd <= 3.0);
    }
}

/// Checks GSL's ten-category multinomial means and conservation of trial counts.
#[test]
fn gsl_multinomial_moments() {
    let mut r = Generator::new();
    let p = [0.2, 0.2, 0.17, 0.14, 0.12, 0.07, 0.05, 0.02, 0.02, 0.01];
    let mut sums = [0.0; 10];
    for _ in 0..N {
        let x = multinomial(&mut r, 100, &p).unwrap();
        assert_eq!(x.iter().sum::<u32>(), 100);
        for (sum, x) in sums.iter_mut().zip(x) {
            *sum += f64::from(x);
        }
    }
    for k in 0..10 {
        let mean = p[k] * 100.0;
        let sd = p[k] * (1.0 - p[k]) * 100.0;
        assert!((N as f64).sqrt() * (sums[k] / N as f64 - mean).abs() / sd <= 3.0);
    }
}

/// Checks the frequencies of all ten values in all ten shuffled positions.
#[test]
fn gsl_shuffle() {
    let mut r = Generator::new();
    let mut count = [[0.0_f64; 10]; 10];
    for _ in 0..N {
        let mut x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        shuffle(&mut r, &mut x);
        for (j, x) in x.into_iter().enumerate() {
            count[x][j] += 1.0;
        }
    }
    for row in count {
        for x in row {
            let expected = N as f64 / 10.0;
            assert!((x - expected).abs() <= 5.0 * expected.sqrt());
        }
    }
}
/// Checks selection frequencies and source-order preservation for three of ten items.
#[test]
fn gsl_choose() {
    let mut r = Generator::new();
    let mut count = [0.0_f64; 10];
    for _ in 0..N {
        let y = choose(&mut r, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 3).unwrap();
        assert!(y.windows(2).all(|w| w[0] < w[1]));
        for x in y {
            count[x] += 1.0;
        }
    }
    for x in count {
        let expected = 3.0 * N as f64 / 10.0;
        assert!((x - expected).abs() <= 5.0 * expected.sqrt());
    }
}

/// Checks the multivariate-normal log density against GSL's R reference value.
#[test]
fn gsl_multivariate_gaussian_log_pdf() {
    relative(
        multivariate_normal_log_density(&[0., 0.], &[1., 2.], &[2., 0., 1., 2.0_f64.sqrt()])
            .unwrap(),
        -3.565097837249263,
        1e-10,
    );
}
/// Checks the multivariate-normal density against GSL's R reference value.
#[test]
fn gsl_multivariate_gaussian_pdf() {
    relative(
        multivariate_normal_density(&[0., 0.], &[1., 2.], &[2., 0., 1., 2.0_f64.sqrt()]).unwrap(),
        0.028294217120391,
        1e-10,
    );
}
/// Checks multivariate-normal sampling with GSL's Hotelling statistic and covariance targets.
#[test]
fn gsl_multivariate_gaussian() {
    let mut r = Generator::new();
    let mut data = Vec::with_capacity(2 * N);
    for _ in 0..N {
        data.extend(multivariate_normal(&mut r, &[1., 2.], &[2., 0., 1., 2.0_f64.sqrt()]).unwrap());
    }
    let mean = multivariate_mean(&data, 2).unwrap();
    let covariance = multivariate_sample_covariance(&data, 2).unwrap();
    let x = mean[0] - 1.0;
    let y = mean[1] - 2.0;
    let t2 = N as f64
        * (covariance[3] * x * x - 2.0 * covariance[1] * x * y + covariance[0] * y * y)
        / (covariance[0] * covariance[3] - covariance[1] * covariance[2]);
    // For F(2,v), Q(f)=(v/(v+2f))^(v/2); invert it directly.
    let v = (N - 2) as f64;
    let f95 = 0.5 * v * ((-2.0 / v * 0.05_f64.ln()).exp() - 1.0);
    assert!(t2 < 2.0 * (N - 1) as f64 / v * f95);
    assert!((v / (v + 2.0 * t2)).powf(v / 2.0) > 0.05);
    for (&x, y) in covariance.iter().zip([4., 2., 2., 3.]) {
        relative(x, y, 0.02);
    }
}

/// Returns the row-major lower Cholesky factor of GSL's two-dimensional Wishart fixture.
///
/// The factored observation matrix has diagonal `[2.213322, 3.285779]` and
/// both off-diagonal entries equal to `1.453357`.
fn wishart_reference_factor() -> [f64; 4] {
    let a = 2.213322_f64.sqrt();
    let b = 1.453357 / a;
    [a, 0., b, (3.285779 - b * b).sqrt()]
}
/// Checks the Wishart log density against GSL's R reference value.
#[test]
fn gsl_wishart_log_pdf() {
    relative(
        wishart_log_density(
            &wishart_reference_factor(),
            3.,
            &[1., 0., 0.3, 0.91_f64.sqrt()],
            2,
        )
        .unwrap(),
        -4.931913612377813,
        1e-10,
    );
}
/// Checks the Wishart density against GSL's R reference value.
#[test]
fn gsl_wishart_pdf() {
    relative(
        wishart_density(
            &wishart_reference_factor(),
            3.,
            &[1., 0., 0.3, 0.91_f64.sqrt()],
            2,
        )
        .unwrap(),
        0.007212687778224,
        1e-10,
    );
}
/// Checks that one-dimensional Wishart samples have the mean and variance of chi-squared samples.
#[test]
fn gsl_wishart() {
    let mut r = Generator::new();
    for df in [1., 3., 5., 7., 9.] {
        let samples: Vec<_> = (0..N)
            .map(|_| wishart(&mut r, df, &[1.], 1).unwrap()[0])
            .collect();
        relative(gsl_oxide_statistics::mean(&samples).unwrap(), df, 0.1);
        relative(
            gsl_oxide_statistics::sample_variance(&samples).unwrap(),
            2.0 * df,
            0.1,
        );
    }
}

/// Checks GSL's sampling histogram for `beta(r, 2.0, 3.0)`.
#[test]
fn gsl_beta() {
    support::continuous(
        |r| beta(r, 2.0, 3.0).unwrap(),
        |x| beta_density(x, 2.0, 3.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `cauchy(r, 2.0)`.
#[test]
fn gsl_cauchy() {
    support::continuous(
        |r| cauchy(r, 2.0).unwrap(),
        |x| cauchy_density(x, 2.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `chi_squared(r, 13.0)`.
#[test]
fn gsl_chisq() {
    support::continuous(
        |r| chi_squared(r, 13.0).unwrap(),
        |x| chi_squared_density(x, 13.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `chi_squared(r, 2.0)`.
#[test]
fn gsl_chisqnu2() {
    support::continuous(
        |r| chi_squared(r, 2.0).unwrap(),
        |x| chi_squared_density(x, 2.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `erlang(r, 3.0, 4.0)`.
#[test]
fn gsl_erlang() {
    support::continuous(
        |r| erlang(r, 3.0, 4.0).unwrap(),
        |x| erlang_density(x, 3.0, 4.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `exponential(r, 2.0)`.
#[test]
fn gsl_exponential() {
    support::continuous(
        |r| exponential(r, 2.0).unwrap(),
        |x| exponential_density(x, 2.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `exponential_power(r, 3.7, 0.3)`.
#[test]
fn gsl_exppow0() {
    support::continuous(
        |r| exponential_power(r, 3.7, 0.3).unwrap(),
        |x| exponential_power_density(x, 3.7, 0.3).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `exponential_power(r, 3.7, 1.0)`.
#[test]
fn gsl_exppow1() {
    support::continuous(
        |r| exponential_power(r, 3.7, 1.0).unwrap(),
        |x| exponential_power_density(x, 3.7, 1.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `exponential_power(r, 3.7, 1.9)`.
#[test]
fn gsl_exppow1a() {
    support::continuous(
        |r| exponential_power(r, 3.7, 1.9).unwrap(),
        |x| exponential_power_density(x, 3.7, 1.9).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `exponential_power(r, 3.7, 2.0)`.
#[test]
fn gsl_exppow2() {
    support::continuous(
        |r| exponential_power(r, 3.7, 2.0).unwrap(),
        |x| exponential_power_density(x, 3.7, 2.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `exponential_power(r, 3.7, 3.5)`.
#[test]
fn gsl_exppow2a() {
    support::continuous(
        |r| exponential_power(r, 3.7, 3.5).unwrap(),
        |x| exponential_power_density(x, 3.7, 3.5).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `exponential_power(r, 3.7, 7.5)`.
#[test]
fn gsl_exppow2b() {
    support::continuous(
        |r| exponential_power(r, 3.7, 7.5).unwrap(),
        |x| exponential_power_density(x, 3.7, 7.5).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `fisher_f(r, 3.0, 4.0)`.
#[test]
fn gsl_fdist() {
    support::continuous(
        |r| fisher_f(r, 3.0, 4.0).unwrap(),
        |x| fisher_f_density(x, 3.0, 4.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `fisher_f(r, 8.0, 249.0)`.
#[test]
fn gsl_fdist_large() {
    support::continuous(
        |r| fisher_f(r, 8.0, 249.0).unwrap(),
        |x| fisher_f_density(x, 8.0, 249.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `uniform_interval(r, 3.0, 4.0)`.
#[test]
fn gsl_flat() {
    support::continuous(
        |r| uniform_interval(r, 3.0, 4.0).unwrap(),
        |x| uniform_interval_density(x, 3.0, 4.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 2.5, 2.17)`.
#[test]
fn gsl_gamma() {
    support::continuous(
        |r| gamma(r, 2.5, 2.17).unwrap(),
        |x| gamma_density(x, 2.5, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 1.0, 2.17)`.
#[test]
fn gsl_gamma1() {
    support::continuous(
        |r| gamma(r, 1.0, 2.17).unwrap(),
        |x| gamma_density(x, 1.0, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 10.0, 2.17)`.
#[test]
fn gsl_gamma_int() {
    support::continuous(
        |r| gamma(r, 10.0, 2.17).unwrap(),
        |x| gamma_density(x, 10.0, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 20.0, 2.17)`.
#[test]
fn gsl_gamma_large() {
    support::continuous(
        |r| gamma(r, 20.0, 2.17).unwrap(),
        |x| gamma_density(x, 20.0, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 0.92, 2.17)`.
#[test]
fn gsl_gamma_small() {
    support::continuous(
        |r| gamma(r, 0.92, 2.17).unwrap(),
        |x| gamma_density(x, 0.92, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 2.5, 2.17)`.
#[test]
fn gsl_gamma_mt() {
    support::continuous(
        |r| gamma(r, 2.5, 2.17).unwrap(),
        |x| gamma_density(x, 2.5, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 1.0, 2.17)`.
#[test]
fn gsl_gamma_mt1() {
    support::continuous(
        |r| gamma(r, 1.0, 2.17).unwrap(),
        |x| gamma_density(x, 1.0, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 10.0, 2.17)`.
#[test]
fn gsl_gamma_mt_int() {
    support::continuous(
        |r| gamma(r, 10.0, 2.17).unwrap(),
        |x| gamma_density(x, 10.0, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 20.0, 2.17)`.
#[test]
fn gsl_gamma_mt_large() {
    support::continuous(
        |r| gamma(r, 20.0, 2.17).unwrap(),
        |x| gamma_density(x, 20.0, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gamma(r, 0.92, 2.17)`.
#[test]
fn gsl_gamma_mt_small() {
    support::continuous(
        |r| gamma(r, 0.92, 2.17).unwrap(),
        |x| gamma_density(x, 0.92, 2.17).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `normal(r, 3.0)`.
#[test]
fn gsl_gaussian() {
    support::continuous(
        |r| normal(r, 3.0).unwrap(),
        |x| normal_density(x, 3.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `normal_ratio(r, 3.0)`.
#[test]
fn gsl_gaussian_ratio_method() {
    support::continuous(
        |r| normal_ratio(r, 3.0).unwrap(),
        |x| normal_density(x, 3.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `normal_ziggurat(r, 3.12)`.
#[test]
fn gsl_gaussian_ziggurat() {
    support::continuous(
        |r| normal_ziggurat(r, 3.12).unwrap(),
        |x| normal_density(x, 3.12).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `standard_normal(r)`.
#[test]
fn gsl_ugaussian() {
    support::continuous(standard_normal, |x| standard_normal_density(x).unwrap());
}

/// Checks GSL's sampling histogram for `standard_normal_ratio(r)`.
#[test]
fn gsl_ugaussian_ratio_method() {
    support::continuous(standard_normal_ratio, |x| {
        standard_normal_density(x).unwrap()
    });
}

/// Checks GSL's sampling histogram for `normal_tail(r, 1.7, 0.25)`.
#[test]
fn gsl_gaussian_tail() {
    support::continuous(
        |r| normal_tail(r, 1.7, 0.25).unwrap(),
        |x| normal_tail_density(x, 1.7, 0.25).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `normal_tail(r, -1.7, 5.0)`.
#[test]
fn gsl_gaussian_tail1() {
    support::continuous(
        |r| normal_tail(r, -1.7, 5.0).unwrap(),
        |x| normal_tail_density(x, -1.7, 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `normal_tail(r, 0.1, 2.0)`.
#[test]
fn gsl_gaussian_tail2() {
    support::continuous(
        |r| normal_tail(r, 0.1, 2.0).unwrap(),
        |x| normal_tail_density(x, 0.1, 2.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `standard_normal_tail(r, 3.0)`.
#[test]
fn gsl_ugaussian_tail() {
    support::continuous(
        |r| standard_normal_tail(r, 3.0).unwrap(),
        |x| standard_normal_tail_density(x, 3.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gumbel_type_1(r, 3.12, 4.56)`.
#[test]
fn gsl_gumbel1() {
    support::continuous(
        |r| gumbel_type_1(r, 3.12, 4.56).unwrap(),
        |x| gumbel_type_1_density(x, 3.12, 4.56).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `gumbel_type_2(r, 3.12, 4.56)`.
#[test]
fn gsl_gumbel2() {
    support::continuous(
        |r| gumbel_type_2(r, 3.12, 4.56).unwrap(),
        |x| gumbel_type_2_density(x, 3.12, 4.56).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `landau(r)`.
#[test]
fn gsl_landau() {
    support::continuous(landau, |x| landau_density(x).unwrap());
}

/// Checks GSL's sampling histogram for `stable_symmetric(r, 5.0, 1.0)`.
#[test]
fn gsl_levy1() {
    support::continuous(
        |r| stable_symmetric(r, 5.0, 1.0).unwrap(),
        |x| cauchy_density(x, 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `stable_symmetric(r, 5.0, 2.0)`.
#[test]
fn gsl_levy2() {
    support::continuous(
        |r| stable_symmetric(r, 5.0, 2.0).unwrap(),
        |x| normal_density(x, 2.0_f64.sqrt() * 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `stable_symmetric(r, 5.0, 1.01)`.
#[test]
fn gsl_levy1a() {
    support::continuous(
        |r| stable_symmetric(r, 5.0, 1.01).unwrap(),
        |x| cauchy_density(x, 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `stable_symmetric(r, 5.0, 1.99)`.
#[test]
fn gsl_levy2a() {
    support::continuous(
        |r| stable_symmetric(r, 5.0, 1.99).unwrap(),
        |x| normal_density(x, 2.0_f64.sqrt() * 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `stable_skewed(r, 5.0, 1.0, 0.0)`.
#[test]
fn gsl_levy_skew1() {
    support::continuous(
        |r| stable_skewed(r, 5.0, 1.0, 0.0).unwrap(),
        |x| cauchy_density(x, 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `stable_skewed(r, 5.0, 2.0, 0.0)`.
#[test]
fn gsl_levy_skew2() {
    support::continuous(
        |r| stable_skewed(r, 5.0, 2.0, 0.0).unwrap(),
        |x| normal_density(x, 2.0_f64.sqrt() * 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `stable_skewed(r, 5.0, 1.01, 0.0)`.
#[test]
fn gsl_levy_skew1a() {
    support::continuous(
        |r| stable_skewed(r, 5.0, 1.01, 0.0).unwrap(),
        |x| cauchy_density(x, 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `stable_skewed(r, 5.0, 1.99, 0.0)`.
#[test]
fn gsl_levy_skew2a() {
    support::continuous(
        |r| stable_skewed(r, 5.0, 1.99, 0.0).unwrap(),
        |x| normal_density(x, 2.0_f64.sqrt() * 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `stable_skewed(r, 5.0, 1.01, 0.001)`.
#[test]
fn gsl_levy_skew1b() {
    support::continuous(
        |r| stable_skewed(r, 5.0, 1.01, 0.001).unwrap(),
        |x| cauchy_density(x, 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `stable_skewed(r, 5.0, 1.99, 0.001)`.
#[test]
fn gsl_levy_skew2b() {
    support::continuous(
        |r| stable_skewed(r, 5.0, 1.99, 0.001).unwrap(),
        |x| normal_density(x, 2.0_f64.sqrt() * 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `logistic(r, 3.1)`.
#[test]
fn gsl_logistic() {
    support::continuous(
        |r| logistic(r, 3.1).unwrap(),
        |x| logistic_density(x, 3.1).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `lognormal(r, 2.7, 1.3)`.
#[test]
fn gsl_lognormal() {
    support::continuous(
        |r| lognormal(r, 2.7, 1.3).unwrap(),
        |x| lognormal_density(x, 2.7, 1.3).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `pareto(r, 1.9, 2.75)`.
#[test]
fn gsl_pareto() {
    support::continuous(
        |r| pareto(r, 1.9, 2.75).unwrap(),
        |x| pareto_density(x, 1.9, 2.75).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `rayleigh(r, 1.9)`.
#[test]
fn gsl_rayleigh() {
    support::continuous(
        |r| rayleigh(r, 1.9).unwrap(),
        |x| rayleigh_density(x, 1.9).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `rayleigh_tail(r, 2.7, 1.9)`.
#[test]
fn gsl_rayleigh_tail() {
    support::continuous(
        |r| rayleigh_tail(r, 2.7, 1.9).unwrap(),
        |x| rayleigh_tail_density(x, 2.7, 1.9).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `student_t(r, 1.75)`.
#[test]
fn gsl_tdist1() {
    support::continuous(
        |r| student_t(r, 1.75).unwrap(),
        |x| student_t_density(x, 1.75).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `student_t(r, 12.75)`.
#[test]
fn gsl_tdist2() {
    support::continuous(
        |r| student_t(r, 12.75).unwrap(),
        |x| student_t_density(x, 12.75).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `laplace(r, 2.75)`.
#[test]
fn gsl_laplace() {
    support::continuous(
        |r| laplace(r, 2.75).unwrap(),
        |x| laplace_density(x, 2.75).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `weibull(r, 3.14, 2.75)`.
#[test]
#[allow(clippy::approx_constant)] // GSL tests this scale, not pi.
fn gsl_weibull() {
    support::continuous(
        |r| weibull(r, 3.14, 2.75).unwrap(),
        |x| weibull_density(x, 3.14, 2.75).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `weibull(r, 2.97, 1.0)`.
#[test]
fn gsl_weibull1() {
    support::continuous(
        |r| weibull(r, 2.97, 1.0).unwrap(),
        |x| weibull_density(x, 2.97, 1.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `poisson(r, 5.0)`.
#[test]
fn gsl_poisson() {
    support::discrete(
        |r| poisson(r, 5.0).unwrap(),
        |x| poisson_pdf(x, 5.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `poisson(r, 30.0)`.
#[test]
fn gsl_poisson_large() {
    support::discrete(
        |r| poisson(r, 30.0).unwrap(),
        |x| poisson_pdf(x, 30.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `bernoulli(r, 0.3)`.
#[test]
fn gsl_bernoulli() {
    support::discrete(
        |r| bernoulli(r, 0.3).unwrap(),
        |x| bernoulli_pdf(x, 0.3).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `binomial(r, 0.3, 5)`.
#[test]
fn gsl_binomial() {
    support::discrete(
        |r| binomial(r, 0.3, 5).unwrap(),
        |x| binomial_pdf(x, 0.3, 5).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `binomial(r, 0.0, 8)`.
#[test]
fn gsl_binomial0() {
    support::discrete(
        |r| binomial(r, 0.0, 8).unwrap(),
        |x| binomial_pdf(x, 0.0, 8).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `binomial(r, 1.0, 8)`.
#[test]
fn gsl_binomial1() {
    support::discrete(
        |r| binomial(r, 1.0, 8).unwrap(),
        |x| binomial_pdf(x, 1.0, 8).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `binomial_knuth(r, 0.3, 5)`.
#[test]
fn gsl_binomial_knuth() {
    support::discrete(
        |r| binomial_knuth(r, 0.3, 5).unwrap(),
        |x| binomial_pdf(x, 0.3, 5).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `binomial(r, 0.3, 55)`.
#[test]
fn gsl_binomial_large() {
    support::discrete(
        |r| binomial(r, 0.3, 55).unwrap(),
        |x| binomial_pdf(x, 0.3, 55).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `binomial_knuth(r, 0.3, 55)`.
#[test]
fn gsl_binomial_large_knuth() {
    support::discrete(
        |r| binomial_knuth(r, 0.3, 55).unwrap(),
        |x| binomial_pdf(x, 0.3, 55).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `binomial(r, 0.3, 5500)`.
#[test]
fn gsl_binomial_huge() {
    support::discrete(
        |r| binomial(r, 0.3, 5500).unwrap(),
        |x| binomial_pdf(x, 0.3, 5500).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `binomial_knuth(r, 0.3, 5500)`.
#[test]
fn gsl_binomial_huge_knuth() {
    support::discrete(
        |r| binomial_knuth(r, 0.3, 5500).unwrap(),
        |x| binomial_pdf(x, 0.3, 5500).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `binomial(r, 1e-8, 1_u32 << 31)`.
#[test]
fn gsl_binomial_max() {
    support::discrete(
        |r| binomial(r, 1e-8, 1_u32 << 31).unwrap(),
        |x| binomial_pdf(x, 1e-8, 1_u32 << 31).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `geometric(r, 0.5)`.
#[test]
fn gsl_geometric() {
    support::discrete(
        |r| geometric(r, 0.5).unwrap(),
        |x| geometric_pdf(x, 0.5).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `geometric(r, 1.0)`.
#[test]
fn gsl_geometric1() {
    support::discrete(
        |r| geometric(r, 1.0).unwrap(),
        |x| geometric_pdf(x, 1.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `hypergeometric(r, 5, 7, 4)`.
#[test]
fn gsl_hypergeometric1() {
    support::discrete(
        |r| hypergeometric(r, 5, 7, 4).unwrap(),
        |x| hypergeometric_pdf(x, 5, 7, 4).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `hypergeometric(r, 5, 7, 11)`.
#[test]
fn gsl_hypergeometric2() {
    support::discrete(
        |r| hypergeometric(r, 5, 7, 11).unwrap(),
        |x| hypergeometric_pdf(x, 5, 7, 11).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `hypergeometric(r, 5, 7, 1)`.
#[test]
fn gsl_hypergeometric3() {
    support::discrete(
        |r| hypergeometric(r, 5, 7, 1).unwrap(),
        |x| hypergeometric_pdf(x, 5, 7, 1).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `hypergeometric(r, 5, 7, 20)`.
#[test]
fn gsl_hypergeometric4() {
    support::discrete(
        |r| hypergeometric(r, 5, 7, 20).unwrap(),
        |x| hypergeometric_pdf(x, 5, 7, 20).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `hypergeometric(r, 2, 7, 5)`.
#[test]
fn gsl_hypergeometric5() {
    support::discrete(
        |r| hypergeometric(r, 2, 7, 5).unwrap(),
        |x| hypergeometric_pdf(x, 2, 7, 5).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `hypergeometric(r, 2, 10, 3)`.
#[test]
fn gsl_hypergeometric6() {
    support::discrete(
        |r| hypergeometric(r, 2, 10, 3).unwrap(),
        |x| hypergeometric_pdf(x, 2, 10, 3).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `logarithmic_series(r, 0.4)`.
#[test]
fn gsl_logarithmic() {
    support::discrete(
        |r| logarithmic_series(r, 0.4).unwrap(),
        |x| logarithmic_series_pdf(x, 0.4).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `negative_binomial(r, 0.3, 20.0)`.
#[test]
fn gsl_negative_binomial() {
    support::discrete(
        |r| negative_binomial(r, 0.3, 20.0).unwrap(),
        |x| negative_binomial_pdf(x, 0.3, 20.0).unwrap(),
    );
}

/// Checks GSL's sampling histogram for `pascal(r, 0.8, 3)`.
#[test]
fn gsl_pascal() {
    support::discrete(
        |r| pascal(r, 0.8, 3).unwrap(),
        |x| pascal_pdf(x, 0.8, 3).unwrap(),
    );
}

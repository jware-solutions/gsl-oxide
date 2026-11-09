use super::support::{relative, Generator};
use gsl_oxide_randist::*;

/// Checks that invalid scalar parameters, weights, and matrix dimensions return errors.
#[test]
fn invalid_parameters_and_dimensions_are_checked() {
    let mut r = Generator::new();
    assert_eq!(normal(&mut r, 0.0), Err(Error::Domain));
    assert_eq!(normal_ratio(&mut r, -1.0), Err(Error::Domain));
    assert_eq!(normal_ziggurat(&mut r, f64::NAN), Err(Error::Domain));
    assert_eq!(normal_tail(&mut r, f64::INFINITY, 1.0), Err(Error::Domain));
    assert_eq!(beta(&mut r, 0.0, 1.0), Err(Error::Domain));
    assert_eq!(gamma(&mut r, -1.0, 1.0), Err(Error::Domain));
    assert_eq!(gamma_knuth(&mut r, 1.0, f64::INFINITY), Err(Error::Domain));
    assert_eq!(bernoulli(&mut r, 1.1), Err(Error::Domain));
    assert_eq!(binomial(&mut r, f64::NAN, 1), Err(Error::Domain));
    assert_eq!(geometric(&mut r, 0.0), Err(Error::Domain));
    assert_eq!(logarithmic_series(&mut r, 1.0), Err(Error::Domain));
    assert_eq!(negative_binomial(&mut r, 0.0, 1.0), Err(Error::Domain));
    assert_eq!(pascal(&mut r, 0.5, 0), Err(Error::Domain));
    assert_eq!(poisson(&mut r, -1.0), Err(Error::Domain));
    assert_eq!(fill_poisson(&mut r, &mut [], f64::NAN), Err(Error::Domain));
    assert_eq!(uniform_interval(&mut r, 1.0, 1.0), Err(Error::Domain));
    assert_eq!(bivariate_normal(&mut r, 1.0, 1.0, 1.1), Err(Error::Domain));
    assert_eq!(
        bivariate_normal_density(0.0, 0.0, 1.0, 1.0, -1.0),
        Err(Error::Domain)
    );
    assert_eq!(stable_skewed(&mut r, 1.0, 2.1, 0.0), Err(Error::Domain));
    assert_eq!(stable_skewed(&mut r, 1.0, 1.5, 1.1), Err(Error::Domain));
    assert_eq!(dirichlet(&mut r, &[]), Err(Error::EmptyInput));
    assert_eq!(dirichlet(&mut r, &[1.0, -1.0]), Err(Error::Domain));
    assert_eq!(dirichlet_density(&[1.0], &[f64::NAN]), Err(Error::Domain));
    assert_eq!(
        dirichlet_density(&[1.0], &[0.5, 0.5]),
        Err(Error::DimensionMismatch)
    );
    assert_eq!(dirichlet_density(&[1.0, 1.0], &[0.4, 0.4]), Ok(0.0));
    assert!(matches!(
        DiscreteDistribution::new(&[]),
        Err(Error::EmptyInput)
    ));
    for weights in [&[0.0, 0.0][..], &[-1.0, 1.0], &[1.0, f64::NAN]] {
        assert!(matches!(
            DiscreteDistribution::new(weights),
            Err(Error::Domain)
        ));
    }
    assert_eq!(multinomial(&mut r, 5, &[0.0]), Err(Error::Domain));
    assert_eq!(
        multinomial_pdf(&[1.0, 1.0], &[1]),
        Err(Error::DimensionMismatch)
    );
    assert_eq!(
        multivariate_normal(&mut r, &[0.0, 0.0], &[1.0]),
        Err(Error::DimensionMismatch)
    );
    assert_eq!(
        multivariate_normal(&mut r, &[0.0], &[-1.0]),
        Err(Error::Domain)
    );
    assert_eq!(
        multivariate_normal_log_density(&[], &[], &[]),
        Err(Error::EmptyInput)
    );
    assert_eq!(
        multivariate_normal_density(&[0.0], &[0.0, 0.0], &[1.0]),
        Err(Error::DimensionMismatch)
    );
    assert_eq!(
        multivariate_mean(&[0.0, 1.0, 2.0], 2),
        Err(Error::DimensionMismatch)
    );
    assert_eq!(
        multivariate_sample_covariance(&[0.0], 1),
        Err(Error::DimensionMismatch)
    );
    assert_eq!(
        wishart(&mut r, 1.0, &[1.0, 0.0, 0.0, 1.0], 2),
        Err(Error::Domain)
    );
    assert_eq!(
        wishart_log_density(&[1.0], 1.0, &[0.0], 1),
        Err(Error::Domain)
    );
    assert_eq!(wishart_density(&[], 1.0, &[], 0), Err(Error::EmptyInput));
    assert_eq!(choose(&mut r, &[1, 2], 3), Err(Error::SampleTooLarge));
    assert_eq!(
        sample_with_replacement::<_, i32>(&mut r, &[], 1),
        Err(Error::EmptyInput)
    );
    assert_eq!(direction_nd(&mut r, 0), Err(Error::EmptyInput));
}

/// Checks support boundaries, nonfinite observations, and degenerate discrete PDFs.
#[test]
fn scalar_support_and_nonfinite_observations() {
    let densities: [fn(f64) -> Result<f64, Error>; 23] = [
        |x| normal_density(x, 1.0),
        standard_normal_density,
        |x| normal_tail_density(x, 0.0, 1.0),
        |x| standard_normal_tail_density(x, 0.0),
        |x| gamma_density(x, 2.0, 1.0),
        |x| beta_density(x, 2.0, 2.0),
        |x| exponential_density(x, 1.0),
        |x| exponential_power_density(x, 1.0, 0.3),
        |x| cauchy_density(x, 1.0),
        |x| chi_squared_density(x, 3.0),
        |x| erlang_density(x, 1.0, 2.0),
        |x| fisher_f_density(x, 2.0, 3.0),
        |x| uniform_interval_density(x, -1.0, 1.0),
        |x| gumbel_type_1_density(x, 1.0, 1.0),
        |x| gumbel_type_2_density(x, 1.0, 1.0),
        |x| logistic_density(x, 1.0),
        |x| lognormal_density(x, 0.0, 1.0),
        |x| pareto_density(x, 2.0, 1.0),
        |x| rayleigh_density(x, 1.0),
        |x| rayleigh_tail_density(x, 1.0, 1.0),
        |x| student_t_density(x, 3.0),
        |x| laplace_density(x, 1.0),
        |x| weibull_density(x, 1.0, 1.0),
    ];
    for pdf in densities {
        assert_eq!(pdf(f64::NAN), Err(Error::Domain));
        assert_eq!(pdf(f64::INFINITY), Ok(0.0));
        assert_eq!(pdf(f64::NEG_INFINITY), Ok(0.0));
    }
    assert_eq!(landau_density(f64::NAN), Err(Error::Domain));
    assert_eq!(landau_density(f64::NEG_INFINITY), Ok(0.0));
    assert_eq!(landau_density(-1000.0), Ok(0.0));
    assert_eq!(beta_density(0.0, 0.5, 2.0), Ok(f64::INFINITY));
    relative(beta_density(0.0, 1.0, 2.0).unwrap(), 2.0, 1e-14);
    assert_eq!(beta_density(1.0, 2.0, 0.5), Ok(f64::INFINITY));
    assert_eq!(gamma_density(0.0, 0.5, 1.0), Ok(f64::INFINITY));
    assert_eq!(gamma_density(0.0, 1.0, 2.0), Ok(0.5));
    assert_eq!(weibull_density(0.0, 1.0, 0.5), Ok(f64::INFINITY));
    assert_eq!(fisher_f_density(0.0, 1.0, 5.0), Ok(f64::INFINITY));
    assert_eq!(uniform_interval_density(1.0, 0.0, 1.0), Ok(0.0));
    assert_eq!(geometric_pdf(0, 0.5), Ok(0.0));
    assert_eq!(logarithmic_series_pdf(0, 0.5), Ok(0.0));
    assert_eq!(bernoulli_pdf(2, 0.5), Ok(0.0));
    assert_eq!(binomial_pdf(11, 0.5, 10), Ok(0.0));
    assert_eq!(hypergeometric_pdf(1, 0, 0, 5), Ok(0.0));
    assert_eq!(poisson_pdf(0, 0.0), Ok(1.0));
    assert_eq!(poisson_pdf(1, 0.0), Ok(0.0));
    assert_eq!(negative_binomial_pdf(0, 1.0, 2.5), Ok(1.0));
    assert_eq!(negative_binomial_pdf(1, 1.0, 2.5), Ok(0.0));
    assert_eq!(multinomial_pdf(&[0.0, 1.0], &[0, 5]), Ok(1.0));
    assert_eq!(multinomial_pdf(&[0.0, 1.0], &[1, 4]), Ok(0.0));
}

/// Checks that finite sums of the tested discrete PDFs equal one within tolerance.
#[test]
fn normalized_probability_masses() {
    for p in [0.0, 0.01, 0.3, 0.5, 0.99, 1.0] {
        for n in [0, 1, 10, 100, 1000] {
            let sum: f64 = (0..=n).map(|k| binomial_pdf(k, p, n).unwrap()).sum();
            relative(sum, 1.0, 2e-11);
        }
    }
    for (a, b, t) in [(0, 0, 0), (5, 7, 4), (5, 7, 20), (20, 2, 7)] {
        let sum: f64 = (0..=a)
            .map(|k| hypergeometric_pdf(k, a, b, t).unwrap())
            .sum();
        relative(sum, 1.0, 1e-13);
    }
    for mean in [0.0, 0.1, 5.0, 30.0, 1000.0] {
        let sum: f64 = (0..2000).map(|k| poisson_pdf(k, mean).unwrap()).sum();
        relative(sum, 1.0, 1e-12);
    }
    relative(
        (1..1000)
            .map(|k| logarithmic_series_pdf(k, 0.8).unwrap())
            .sum(),
        1.0,
        1e-13,
    );
    relative(
        (0..1000)
            .map(|k| negative_binomial_pdf(k, 0.3, 2.7).unwrap())
            .sum(),
        1.0,
        1e-12,
    );
}

/// Checks extreme ranges, narrow intervals, and logarithmic density preservation.
#[test]
fn extreme_parameters_and_logarithmic_densities() {
    let mut r = Generator::new();
    assert_eq!(poisson(&mut r, f64::MAX), Err(Error::NumericalRange));
    assert_eq!(
        geometric(&mut r, f64::MIN_POSITIVE),
        Err(Error::NumericalRange)
    );
    let d = DiscreteDistribution::new(&[f64::MAX, 0.0, f64::MAX]).unwrap();
    assert_eq!(d.pdf(0), 0.5);
    assert_eq!(d.pdf(1), 0.0);
    assert_eq!(d.pdf(3), 0.0);
    for _ in 0..1000 {
        assert_ne!(d.sample(&mut r), 1);
    }
    relative(
        normal_tail_density(40.0, 40.0, 1.0).unwrap(),
        40.024_968_847_207_26,
        2e-14,
    );
    assert!(normal_tail(&mut r, 40.0, 1.0).unwrap() >= 40.0);
    relative(
        gamma_density(4294967296.0, 4294967296.0, 1.0).unwrap(),
        1.0 / (65536.0 * (2.0 * std::f64::consts::PI).sqrt()),
        1e-10,
    );
    assert!(fisher_f_density(1.0, 8.0, 249.0).unwrap().is_finite());
    let sample = dirichlet(&mut r, &[1e-10; 10]).unwrap();
    assert_eq!(sample.iter().sum::<f64>(), 1.0);
    for _ in 0..100 {
        let x = beta(&mut r, 1e-5, 1e-5).unwrap();
        assert!((0.0..=1.0).contains(&x));
        let x = uniform_interval(&mut r, -f64::MAX, f64::MAX).unwrap();
        assert!(x.is_finite());
        let a = 1.0;
        let b = f64::from_bits(1.0_f64.to_bits() + 1);
        assert_eq!(uniform_interval(&mut r, a, b), Ok(a));
    }
    let log = multivariate_normal_log_density(&[100.0], &[0.0], &[1.0]).unwrap();
    assert!(log.is_finite() && log < -5000.0);
    assert_eq!(
        multivariate_normal_density(&[100.0], &[0.0], &[1.0]),
        Ok(0.0)
    );
    let factor = [2.0, f64::NAN, 1.0, 1.0];
    assert!(multivariate_normal(&mut r, &[0.0, 0.0], &factor).is_ok());
}

/// Checks Poisson slice filling, integer-gamma moments, and selection frequencies.
#[test]
fn poisson_arrays_integer_gamma_and_replacement_sampling() {
    let mut r = Generator::new();
    let mut values = [99; 1024];
    fill_poisson(&mut r, &mut values, 0.0).unwrap();
    assert_eq!(values, [0; 1024]);
    fill_poisson(&mut r, &mut values, 5.0).unwrap();
    relative(
        values.iter().map(|&x| f64::from(x)).sum::<f64>() / 1024.0,
        5.0,
        0.1,
    );
    assert_eq!(gamma_integer(&mut r, 0), Ok(0.0));
    for shape in [1, 5, 12, 100] {
        let samples: Vec<_> = (0..20_000)
            .map(|_| gamma_integer(&mut r, shape).unwrap())
            .collect();
        relative(
            gsl_oxide_statistics::mean(&samples).unwrap(),
            f64::from(shape),
            0.03,
        );
        relative(
            gsl_oxide_statistics::sample_variance(&samples).unwrap(),
            f64::from(shape),
            0.05,
        );
    }
    let sample = sample_with_replacement(&mut r, &[0, 1, 2, 3], 100_000).unwrap();
    let mut counts = [0; 4];
    for x in sample {
        counts[x] += 1;
    }
    for x in counts {
        assert!((f64::from(x) - 25_000.0).abs() < 5.0 * 25000.0_f64.sqrt());
    }
    assert_eq!(choose(&mut r, &[1, 2], 2), Ok(vec![1, 2]));
    assert_eq!(choose::<_, i32>(&mut r, &[], 0), Ok(vec![]));
    assert_eq!(
        sample_with_replacement::<_, i32>(&mut r, &[], 0),
        Ok(vec![])
    );
    shuffle::<_, i32>(&mut r, &mut []);
}

/// Checks unit lengths and isotropic coordinate moments of direction samples.
#[test]
fn directions_have_unit_norm_and_isotropic_moments() {
    let mut r = Generator::new();
    for d in [1, 2, 3, 10] {
        let mut sums = vec![0.0; d];
        let mut squares = vec![0.0; d];
        for _ in 0..20_000 {
            let x = direction_nd(&mut r, d).unwrap();
            relative(x.iter().map(|x| x * x).sum(), 1.0, 1e-14);
            for i in 0..d {
                sums[i] += x[i];
                squares[i] += x[i] * x[i];
            }
        }
        for i in 0..d {
            assert!((sums[i] / 20_000.0).abs() < 5.0 / (20_000.0 * d as f64).sqrt());
            relative(squares[i] / 20_000.0, 1.0 / d as f64, 0.05);
        }
    }
    for _ in 0..100 {
        for v in [direction_2d(&mut r), direction_2d_trigonometric(&mut r)] {
            relative(v.iter().map(|x| x * x).sum(), 1.0, 1e-14);
        }
        relative(direction_3d(&mut r).iter().map(|x| x * x).sum(), 1.0, 1e-14);
    }
}

/// Checks three-dimensional Wishart moments and multivariate-estimator reference values.
#[test]
fn higher_dimensional_wishart_moments_and_covariance_estimates() {
    let mut r = Generator::new();
    let l = [1., 0., 0., 0.3, 2., 0., -0.2, 0.5, 1.];
    let expected = [4., 1.2, -0.8, 1.2, 16.36, 3.76, -0.8, 3.76, 5.16];
    let mut sums = [0.0; 9];
    for _ in 0..30_000 {
        let x = wishart(&mut r, 4.0, &l, 3).unwrap();
        assert!(x[0] > 0.0 && x[4] > 0.0 && x[8] > 0.0);
        assert_eq!(x[1], x[3]);
        assert_eq!(x[2], x[6]);
        assert_eq!(x[5], x[7]);
        assert!(x[0] * x[4] - x[1] * x[1] > 0.0);
        for i in 0..9 {
            sums[i] += x[i];
        }
    }
    for i in 0..9 {
        relative(sums[i] / 30_000.0, expected[i], 0.06);
    }
    let x = [4.348817, 4.711934, 2.995049, 1.190864, -3.793431, -1.357363];
    let mean = multivariate_mean(&x, 2).unwrap();
    relative(mean[0], 1.1834783333333332, 1e-12);
    relative(mean[1], 1.515145, 1e-12);
    let covariance = multivariate_sample_covariance(&x, 2).unwrap();
    for (actual, expected) in covariance
        .iter()
        .zip([19.03539, 11.91384, 11.91384, 9.28796])
    {
        relative(*actual, expected, 1e-6);
    }
}

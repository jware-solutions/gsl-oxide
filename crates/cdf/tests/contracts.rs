//! Rust domain contracts, support endpoints, and inverse checks for every CDF family.
use gsl_oxide_cdf::*;
mod support;
use support::assert_relative;

type Function = fn(f64) -> Result<f64, Error>;
struct Continuous {
    name: &'static str,
    cdf: Function,
    survival: Function,
    quantile: Function,
    inverse_survival: Function,
    minimum: f64,
    maximum: f64,
}
macro_rules! case {
    ($cdf:ident, $sf:ident, $inv:ident, $isf:ident, ($($p:expr),*), $min:expr, $max:expr) => {
        Continuous { name: stringify!($cdf), cdf: |x| $cdf(x, $($p),*), survival: |x| $sf(x, $($p),*),
            quantile: |p| $inv(p, $($p),*), inverse_survival: |p| $isf(p, $($p),*), minimum: $min, maximum: $max }
    };
}

fn continuous_cases() -> Vec<Continuous> {
    let inf = f64::INFINITY;
    vec![
        Continuous {
            name: "standard_normal",
            cdf: standard_normal_cdf,
            survival: standard_normal_survival,
            quantile: standard_normal_quantile,
            inverse_survival: standard_normal_inverse_survival,
            minimum: -inf,
            maximum: inf,
        },
        case!(
            normal_cdf,
            normal_survival,
            normal_quantile,
            normal_inverse_survival,
            (2.3),
            -inf,
            inf
        ),
        case!(
            gamma_cdf,
            gamma_survival,
            gamma_quantile,
            gamma_inverse_survival,
            (0.7, 2.3),
            0.0,
            inf
        ),
        case!(
            gamma_cdf,
            gamma_survival,
            gamma_quantile,
            gamma_inverse_survival,
            (27.0, 0.7),
            0.0,
            inf
        ),
        case!(
            beta_cdf,
            beta_survival,
            beta_quantile,
            beta_inverse_survival,
            (0.7, 1.3),
            0.0,
            1.0
        ),
        case!(
            chi_squared_cdf,
            chi_squared_survival,
            chi_squared_quantile,
            chi_squared_inverse_survival,
            (13.0),
            0.0,
            inf
        ),
        case!(
            student_t_cdf,
            student_t_survival,
            student_t_quantile,
            student_t_inverse_survival,
            (0.7),
            -inf,
            inf
        ),
        case!(
            student_t_cdf,
            student_t_survival,
            student_t_quantile,
            student_t_inverse_survival,
            (300.0),
            -inf,
            inf
        ),
        case!(
            fisher_f_cdf,
            fisher_f_survival,
            fisher_f_quantile,
            fisher_f_inverse_survival,
            (0.7, 1.3),
            0.0,
            inf
        ),
        case!(
            exponential_cdf,
            exponential_survival,
            exponential_quantile,
            exponential_inverse_survival,
            (2.3),
            0.0,
            inf
        ),
        case!(
            cauchy_cdf,
            cauchy_survival,
            cauchy_quantile,
            cauchy_inverse_survival,
            (2.3),
            -inf,
            inf
        ),
        case!(
            laplace_cdf,
            laplace_survival,
            laplace_quantile,
            laplace_inverse_survival,
            (2.3),
            -inf,
            inf
        ),
        case!(
            rayleigh_cdf,
            rayleigh_survival,
            rayleigh_quantile,
            rayleigh_inverse_survival,
            (2.3),
            0.0,
            inf
        ),
        case!(
            uniform_interval_cdf,
            uniform_interval_survival,
            uniform_interval_quantile,
            uniform_interval_inverse_survival,
            (-2.3, 4.7),
            -2.3,
            4.7
        ),
        case!(
            lognormal_cdf,
            lognormal_survival,
            lognormal_quantile,
            lognormal_inverse_survival,
            (1.3, 0.7),
            0.0,
            inf
        ),
        case!(
            gumbel_type_1_cdf,
            gumbel_type_1_survival,
            gumbel_type_1_quantile,
            gumbel_type_1_inverse_survival,
            (1.3, 2.7),
            -inf,
            inf
        ),
        case!(
            gumbel_type_2_cdf,
            gumbel_type_2_survival,
            gumbel_type_2_quantile,
            gumbel_type_2_inverse_survival,
            (1.3, 2.7),
            0.0,
            inf
        ),
        case!(
            weibull_cdf,
            weibull_survival,
            weibull_quantile,
            weibull_inverse_survival,
            (2.3, 0.7),
            0.0,
            inf
        ),
        case!(
            pareto_cdf,
            pareto_survival,
            pareto_quantile,
            pareto_inverse_survival,
            (1.3, 2.7),
            2.7,
            inf
        ),
        case!(
            logistic_cdf,
            logistic_survival,
            logistic_quantile,
            logistic_inverse_survival,
            (2.3),
            -inf,
            inf
        ),
    ]
}

#[test]
fn endpoints_and_inverse_round_trips() {
    for c in continuous_cases() {
        assert_eq!((c.cdf)(f64::NEG_INFINITY), Ok(0.0), "{}", c.name);
        assert_eq!((c.survival)(f64::NEG_INFINITY), Ok(1.0), "{}", c.name);
        assert_eq!((c.cdf)(f64::INFINITY), Ok(1.0), "{}", c.name);
        assert_eq!((c.survival)(f64::INFINITY), Ok(0.0), "{}", c.name);
        assert_eq!((c.quantile)(0.0), Ok(c.minimum), "{}", c.name);
        assert_eq!((c.quantile)(1.0), Ok(c.maximum), "{}", c.name);
        assert_eq!((c.inverse_survival)(0.0), Ok(c.maximum), "{}", c.name);
        assert_eq!((c.inverse_survival)(1.0), Ok(c.minimum), "{}", c.name);
        let mut previous = c.minimum;
        for p in [1e-6, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0 - 1e-6] {
            let x = (c.quantile)(p).unwrap();
            assert!(x >= previous, "{}: non-monotone quantile", c.name);
            previous = x;
            assert_relative((c.cdf)(x).unwrap(), p, 2e-9, c.name);
            let x = (c.inverse_survival)(p).unwrap();
            assert_relative((c.survival)(x).unwrap(), p, 2e-9, c.name);
        }
    }
    for function in [exponential_power_cdf, exponential_power_survival] {
        assert_eq!(function(f64::NAN, 1.0, 2.0), Err(Error::Domain));
    }
    assert_eq!(exponential_power_cdf(f64::NEG_INFINITY, 1.0, 2.0), Ok(0.0));
    assert_eq!(exponential_power_cdf(f64::INFINITY, 1.0, 2.0), Ok(1.0));
    assert_eq!(
        exponential_power_survival(f64::NEG_INFINITY, 1.0, 2.0),
        Ok(1.0)
    );
    assert_eq!(exponential_power_survival(f64::INFINITY, 1.0, 2.0), Ok(0.0));
}

#[test]
fn invalid_observations_and_probabilities_are_errors() {
    for c in continuous_cases() {
        assert_eq!((c.cdf)(f64::NAN), Err(Error::Domain), "{}", c.name);
        assert_eq!((c.survival)(f64::NAN), Err(Error::Domain), "{}", c.name);
        for p in [f64::NAN, f64::NEG_INFINITY, -0.01, 1.01, f64::INFINITY] {
            assert_eq!((c.quantile)(p), Err(Error::Domain), "{}", c.name);
            assert_eq!((c.inverse_survival)(p), Err(Error::Domain), "{}", c.name);
        }
    }
}

#[test]
fn invalid_parameters_are_checked_before_support_limits() {
    for bad in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
        for function in [
            normal_cdf,
            normal_survival,
            normal_quantile,
            normal_inverse_survival,
            exponential_cdf,
            exponential_survival,
            exponential_quantile,
            exponential_inverse_survival,
            rayleigh_cdf,
            rayleigh_survival,
            rayleigh_quantile,
            rayleigh_inverse_survival,
            cauchy_cdf,
            cauchy_survival,
            cauchy_quantile,
            cauchy_inverse_survival,
            laplace_cdf,
            laplace_survival,
            laplace_quantile,
            laplace_inverse_survival,
            logistic_cdf,
            logistic_survival,
            logistic_quantile,
            logistic_inverse_survival,
            chi_squared_cdf,
            chi_squared_survival,
            chi_squared_quantile,
            chi_squared_inverse_survival,
            student_t_cdf,
            student_t_survival,
            student_t_quantile,
            student_t_inverse_survival,
        ] {
            assert_eq!(function(0.0, bad), Err(Error::Domain));
        }
        for function in [
            gamma_cdf,
            gamma_survival,
            gamma_quantile,
            gamma_inverse_survival,
            beta_cdf,
            beta_survival,
            beta_quantile,
            beta_inverse_survival,
            fisher_f_cdf,
            fisher_f_survival,
            fisher_f_quantile,
            fisher_f_inverse_survival,
            gumbel_type_1_cdf,
            gumbel_type_1_survival,
            gumbel_type_1_quantile,
            gumbel_type_1_inverse_survival,
            gumbel_type_2_cdf,
            gumbel_type_2_survival,
            gumbel_type_2_quantile,
            gumbel_type_2_inverse_survival,
            weibull_cdf,
            weibull_survival,
            weibull_quantile,
            weibull_inverse_survival,
            pareto_cdf,
            pareto_survival,
            pareto_quantile,
            pareto_inverse_survival,
            exponential_power_cdf,
            exponential_power_survival,
        ] {
            assert_eq!(function(0.0, bad, 1.0), Err(Error::Domain));
            assert_eq!(function(0.0, 1.0, bad), Err(Error::Domain));
        }
        assert_eq!(lognormal_cdf(0.0, 1.0, bad), Err(Error::Domain));
        assert_eq!(lognormal_quantile(0.0, 1.0, bad), Err(Error::Domain));
    }
    for function in [
        uniform_interval_cdf,
        uniform_interval_survival,
        uniform_interval_quantile,
        uniform_interval_inverse_survival,
    ] {
        assert_eq!(function(0.0, 1.0, 1.0), Err(Error::Domain));
        assert_eq!(function(0.0, 2.0, 1.0), Err(Error::Domain));
        assert_eq!(
            function(0.0, -f64::MAX, f64::MAX),
            Err(Error::UnsupportedRange)
        );
    }
    assert_eq!(gamma_cdf(0.0, 1e8 + 1.0, 1.0), Err(Error::UnsupportedRange));
    assert_eq!(beta_cdf(0.0, 1e8 + 1.0, 1.0), Err(Error::UnsupportedRange));
}

#[test]
fn discrete_edges_and_widened_population_counts() {
    for bad in [f64::NAN, -0.1, 1.1, f64::INFINITY] {
        assert_eq!(binomial_cdf(0, bad, 0), Err(Error::Domain));
        assert_eq!(binomial_survival(0, bad, 0), Err(Error::Domain));
        assert_eq!(geometric_cdf(0, bad), Err(Error::Domain));
        assert_eq!(geometric_survival(0, bad), Err(Error::Domain));
        assert_eq!(negative_binomial_cdf(0, bad, 1.0), Err(Error::Domain));
        assert_eq!(negative_binomial_survival(0, bad, 1.0), Err(Error::Domain));
        assert_eq!(pascal_cdf(0, bad, 1), Err(Error::Domain));
        assert_eq!(pascal_survival(0, bad, 1), Err(Error::Domain));
    }
    assert_eq!(binomial_cdf(u32::MAX, 0.5, u32::MAX), Ok(1.0));
    assert_eq!(binomial_survival(u32::MAX, 0.5, u32::MAX), Ok(0.0));
    assert_eq!(binomial_cdf(0, 0.0, 5), Ok(1.0));
    assert_eq!(binomial_cdf(4, 1.0, 5), Ok(0.0));
    assert_eq!(geometric_cdf(0, 1.0), Ok(0.0));
    assert_eq!(geometric_cdf(1, 1.0), Ok(1.0));
    assert_eq!(geometric_survival(u32::MAX, 0.0), Ok(1.0));
    assert_eq!(negative_binomial_cdf(0, 1.0, 1.0), Ok(1.0));
    assert_eq!(negative_binomial_cdf(0, 0.0, 1.0), Ok(0.0));
    assert_eq!(negative_binomial_cdf(0, 0.5, 0.0), Err(Error::Domain));
    assert_eq!(pascal_cdf(0, 0.5, 0), Err(Error::Domain));
    assert_eq!(poisson_cdf(u32::MAX, 0.0), Ok(1.0));
    assert_eq!(poisson_survival(0, 0.0), Ok(0.0));
    assert_eq!(poisson_cdf(0, -1.0), Err(Error::Domain));
    assert_eq!(poisson_survival(0, f64::INFINITY), Err(Error::Domain));
    assert_eq!(hypergeometric_cdf(0, 0, 0, 0), Ok(1.0));
    assert_eq!(hypergeometric_survival(0, 0, 0, 0), Ok(0.0));
    assert_eq!(hypergeometric_cdf(0, 3, 4, 8), Err(Error::Domain));
    assert_eq!(hypergeometric_survival(0, 3, 4, 8), Err(Error::Domain));
    assert_eq!(hypergeometric_cdf(2, 7, 3, 6), Ok(0.0));
    assert_eq!(hypergeometric_survival(2, 7, 3, 6), Ok(1.0));
    assert_relative(
        hypergeometric_cdf(0, u32::MAX, u32::MAX, 1).unwrap(),
        0.5,
        2e-14,
        "widened population",
    );
}

#[test]
fn very_small_tails_and_overflow() {
    for p in [1e-300, 1e-100, 1e-20] {
        let x = standard_normal_inverse_survival(p).unwrap();
        assert_relative(
            standard_normal_survival(x).unwrap(),
            p,
            2e-12,
            "normal extreme tail",
        );
        let x = gamma_inverse_survival(p, 0.7, 2.3).unwrap();
        assert_relative(
            gamma_survival(x, 0.7, 2.3).unwrap(),
            p,
            2e-12,
            "gamma extreme tail",
        );
        let x = beta_quantile(p, 1.2, 1.3).unwrap();
        assert_relative(
            beta_cdf(x, 1.2, 1.3).unwrap(),
            p,
            2e-12,
            "beta extreme tail",
        );
    }
    assert!(standard_normal_survival(38.0).unwrap() > 0.0);
    assert!(student_t_survival(1e300, 1.0).unwrap() > 0.0);
    assert_eq!(gamma_quantile(0.1, 0.001, 1.0), Ok(0.0));
    assert_eq!(beta_quantile(0.1, 0.001, 1.0), Ok(0.0));
    assert!(gamma_quantile(0.4, 0.001, 1e300).unwrap() > 0.0);
    let x = gamma_quantile(0.4, 0.001, 1e300).unwrap();
    assert_relative(
        gamma_cdf(x, 0.001, 1e300).unwrap(),
        0.4,
        2e-13,
        "scaled gamma underflow",
    );
    assert_relative(
        weibull_cdf(1e300, 1e-300, 0.001).unwrap(),
        weibull_cdf(10.0_f64.powf(0.6), 1.0, 1.0).unwrap(),
        2e-14,
        "Weibull ratio overflow",
    );
    assert_relative(
        weibull_cdf(1e-300, 1e300, 0.001).unwrap(),
        weibull_cdf(10.0_f64.powf(-0.6), 1.0, 1.0).unwrap(),
        2e-14,
        "Weibull ratio underflow",
    );
    assert!(weibull_quantile(0.99, 1e-300, 0.002).unwrap().is_finite());
    assert!(pareto_quantile(0.9, 0.003, 1e-300).unwrap().is_finite());
    assert!(student_t_survival(f64::MAX, 1e-100).unwrap() > 0.0);
    assert!(beta_survival(1e-10, 1e-100, 1e8).unwrap() > 0.0);
    assert_eq!(lognormal_quantile(0.5, 1000.0, 1.0), Err(Error::Overflow));
    assert_eq!(cauchy_quantile(1e-300, 1e300), Err(Error::Overflow));
}

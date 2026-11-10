use crate::{finite_quantile, observation, positive, probability, Error};
use std::f64::consts::PI;

// Generate the four public operations together so validation and parameter
// conventions are identical for both tails and their inverses.
macro_rules! distribution {
    ($cdf:ident, $sf:ident, $inv:ident, $isf:ident,
     ($($param:ident),+), $description:literal, $validate:block,
     |$x:ident| $tails:block, |$p:ident, $upper:ident| $inverse:block) => {
        #[doc = concat!($description, " CDF `Pr(X <= x)`. See the crate's domain and error contract.")]
        pub fn $cdf(x: f64, $($param: f64),+) -> Result<f64, Error> {
            $validate
            observation(x)?;
            let $x = x;
            let tails: (f64, f64) = $tails;
            Ok(tails.0)
        }
        #[doc = concat!($description, " survival `Pr(X > x)`, evaluating small upper tails directly.")]
        pub fn $sf(x: f64, $($param: f64),+) -> Result<f64, Error> {
            $validate
            observation(x)?;
            let $x = x;
            let tails: (f64, f64) = $tails;
            Ok(tails.1)
        }
        #[doc = concat!($description, " quantile for probability `p` in `[0,1]`, including support endpoints.")]
        pub fn $inv(p: f64, $($param: f64),+) -> Result<f64, Error> {
            $validate
            probability(p)?;
            let $p = p;
            let $upper = false;
            finite_quantile($inverse, p)
        }
        #[doc = concat!($description, " survival inverse for probability `p` in `[0,1]`, retaining small tails.")]
        pub fn $isf(p: f64, $($param: f64),+) -> Result<f64, Error> {
            $validate
            probability(p)?;
            let $p = p;
            let $upper = true;
            finite_quantile($inverse, p)
        }
    };
}

/// Returns lower and upper probabilities from `log_survival`, evaluating the lower tail with `expm1`.
fn exponential_tails(log_survival: f64) -> (f64, f64) {
    (-log_survival.exp_m1(), log_survival.exp())
}

/// Returns the log lower-tail probability corresponding to `p`; `upper` identifies a supplied survival probability.
fn lower_log(p: f64, upper: bool) -> f64 {
    if upper {
        (-p).ln_1p()
    } else {
        p.ln()
    }
}

/// Returns the log survival probability corresponding to `p`; `upper` identifies a supplied survival probability.
fn upper_log(p: f64, upper: bool) -> f64 {
    if upper {
        p.ln()
    } else {
        (-p).ln_1p()
    }
}

/// Returns `scale * base^exponent`, switching to logarithms when the power overflows or underflows.
fn scaled_power(scale: f64, base: f64, exponent: f64) -> f64 {
    let power = base.powf(exponent);
    if power == 0.0 || power.is_infinite() {
        (scale.ln() + exponent * base.ln()).exp()
    } else {
        scale * power
    }
}

distribution!(
    exponential_cdf,
    exponential_survival,
    exponential_quantile,
    exponential_inverse_survival,
    (mean),
    "Exponential distribution with positive mean `mean`.",
    {
        positive(mean)?;
    },
    |x| { exponential_tails(-x.max(0.0) / mean) },
    |p, upper| { -mean * upper_log(p, upper) }
);

distribution!(
    rayleigh_cdf,
    rayleigh_survival,
    rayleigh_quantile,
    rayleigh_inverse_survival,
    (sigma),
    "Rayleigh distribution with positive scale `sigma`.",
    {
        positive(sigma)?;
    },
    |x| {
        let u = x.max(0.0) / sigma;
        exponential_tails(-0.5 * u * u)
    },
    |p, upper| { sigma * (-2.0 * upper_log(p, upper)).sqrt() }
);

distribution!(
    weibull_cdf,
    weibull_survival,
    weibull_quantile,
    weibull_inverse_survival,
    (scale, shape),
    "Weibull distribution with positive `scale` and `shape` (GSL a,b).",
    {
        positive(scale)?;
        positive(shape)?;
    },
    |x| {
        let ratio = x.max(0.0) / scale;
        let power = if x > 0.0 && (ratio == 0.0 || ratio.is_infinite()) {
            (shape * (x.ln() - scale.ln())).exp()
        } else {
            ratio.powf(shape)
        };
        exponential_tails(-power)
    },
    |p, upper| { scaled_power(scale, -upper_log(p, upper), 1.0 / shape) }
);

distribution!(
    pareto_cdf,
    pareto_survival,
    pareto_quantile,
    pareto_inverse_survival,
    (shape, minimum),
    "Pareto distribution with positive `shape` and support `[minimum,infinity)`.",
    {
        positive(shape)?;
        positive(minimum)?;
    },
    |x| {
        if x <= minimum {
            (0.0, 1.0)
        } else {
            let ratio = x / minimum;
            let log = if ratio.is_finite() {
                -shape * ratio.ln()
            } else {
                shape * (minimum.ln() - x.ln())
            };
            exponential_tails(log)
        }
    },
    |p, upper| {
        let exponent = -upper_log(p, upper) / shape;
        let power = exponent.exp();
        if power.is_infinite() {
            (minimum.ln() + exponent).exp()
        } else {
            minimum * power
        }
    }
);

distribution!(
    cauchy_cdf,
    cauchy_survival,
    cauchy_quantile,
    cauchy_inverse_survival,
    (scale),
    "Zero-centered Cauchy distribution with positive `scale`.",
    {
        positive(scale)?;
    },
    |x| {
        let u = x / scale;
        let small = if u.abs() > 1.0 {
            (scale / x.abs()).atan() / PI
        } else {
            0.5 - u.abs().atan() / PI
        };
        if x < 0.0 {
            (small, 1.0 - small)
        } else {
            (1.0 - small, small)
        }
    },
    |p, upper| {
        let value = if p == 0.5 {
            0.0
        } else if p < 0.5 {
            -scale / (PI * p).tan()
        } else {
            scale / (PI * (1.0 - p)).tan()
        };
        if upper {
            -value
        } else {
            value
        }
    }
);

distribution!(
    laplace_cdf,
    laplace_survival,
    laplace_quantile,
    laplace_inverse_survival,
    (scale),
    "Zero-centered Laplace distribution with positive `scale`.",
    {
        positive(scale)?;
    },
    |x| {
        let small = 0.5 * (-x.abs() / scale).exp();
        if x < 0.0 {
            (small, 1.0 - small)
        } else {
            (1.0 - small, small)
        }
    },
    |p, upper| {
        let value = if p <= 0.5 {
            scale * (2.0 * p).ln()
        } else {
            -scale * (2.0 * (1.0 - p)).ln()
        };
        if upper {
            -value
        } else {
            value
        }
    }
);

distribution!(
    logistic_cdf,
    logistic_survival,
    logistic_quantile,
    logistic_inverse_survival,
    (scale),
    "Zero-centered logistic distribution with positive `scale`.",
    {
        positive(scale)?;
    },
    |x| {
        let e = (-x.abs() / scale).exp();
        let small = e / (1.0 + e);
        if x < 0.0 {
            (small, 1.0 - small)
        } else {
            (1.0 - small, small)
        }
    },
    |p, upper| {
        let value = scale * (p.ln() - (-p).ln_1p());
        if upper {
            -value
        } else {
            value
        }
    }
);

distribution!(
    gumbel_type_1_cdf,
    gumbel_type_1_survival,
    gumbel_type_1_quantile,
    gumbel_type_1_inverse_survival,
    (rate, weight),
    "Gumbel type I, `P(x)=exp(-weight*exp(-rate*x))`, with positive parameters.",
    {
        positive(rate)?;
        positive(weight)?;
    },
    |x| {
        let log = -(weight.ln() - rate * x).exp();
        (log.exp(), -log.exp_m1())
    },
    |p, upper| { (weight.ln() - (-lower_log(p, upper)).ln()) / rate }
);

distribution!(
    gumbel_type_2_cdf,
    gumbel_type_2_survival,
    gumbel_type_2_quantile,
    gumbel_type_2_inverse_survival,
    (shape, weight),
    "Gumbel type II, `P(x)=exp(-weight*x^(-shape))` for `x>0`, with positive parameters.",
    {
        positive(shape)?;
        positive(weight)?;
    },
    |x| {
        if x <= 0.0 {
            (0.0, 1.0)
        } else {
            let log = -(weight.ln() - shape * x.ln()).exp();
            (log.exp(), -log.exp_m1())
        }
    },
    |p, upper| { ((weight.ln() - (-lower_log(p, upper)).ln()) / shape).exp() }
);

distribution!(
    uniform_interval_cdf,
    uniform_interval_survival,
    uniform_interval_quantile,
    uniform_interval_inverse_survival,
    (lower, upper_bound),
    "Uniform distribution on finite bounds `lower < upper_bound`; width must be finite.",
    {
        if !lower.is_finite() || !upper_bound.is_finite() || lower >= upper_bound {
            return Err(Error::Domain);
        }
        if !(upper_bound - lower).is_finite() {
            return Err(Error::UnsupportedRange);
        }
    },
    |x| {
        if x <= lower {
            (0.0, 1.0)
        } else if x >= upper_bound {
            (1.0, 0.0)
        } else {
            (
                (x - lower) / (upper_bound - lower),
                (upper_bound - x) / (upper_bound - lower),
            )
        }
    },
    |p, upper| {
        if upper {
            p * lower + (1.0 - p) * upper_bound
        } else {
            (1.0 - p) * lower + p * upper_bound
        }
    }
);

distribution!(
    lognormal_cdf,
    lognormal_survival,
    lognormal_quantile,
    lognormal_inverse_survival,
    (log_mean, log_sigma),
    "Lognormal distribution with finite `log_mean` and positive `log_sigma`.",
    {
        observation(log_mean)?;
        if !log_mean.is_finite() {
            return Err(Error::Domain);
        }
        positive(log_sigma)?;
    },
    |x| {
        if x <= 0.0 {
            (0.0, 1.0)
        } else {
            let z = (x.ln() - log_mean) / log_sigma;
            (
                crate::standard_normal_cdf(z)?,
                crate::standard_normal_survival(z)?,
            )
        }
    },
    |p, upper| {
        let z = if upper {
            crate::standard_normal_inverse_survival(p)?
        } else {
            crate::standard_normal_quantile(p)?
        };
        (log_mean + log_sigma * z).exp()
    }
);

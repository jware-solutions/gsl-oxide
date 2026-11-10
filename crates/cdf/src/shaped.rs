use crate::{
    finite_quantile, inverse::positive_quantile, observation, positive, probability, Error,
};
use gsl_oxide_special::{
    beta::ln_beta, gamma::ln_abs_gamma_with_sign, incomplete_beta, incomplete_gamma,
};

macro_rules! family {
    ($cdf:ident, $sf:ident, $inv:ident, $isf:ident, ($($arg:ident),+), $description:literal, $tails:ident, $inverse:ident) => {
        #[doc = concat!($description, " CDF `Pr(X <= x)`. See crate documentation for parameter limits and errors.")]
        pub fn $cdf(x: f64, $($arg: f64),+) -> Result<f64, Error> { Ok($tails(x, $($arg),+)?.0) }
        #[doc = concat!($description, " survival `Pr(X > x)`, retaining small upper tails.")]
        pub fn $sf(x: f64, $($arg: f64),+) -> Result<f64, Error> { Ok($tails(x, $($arg),+)?.1) }
        #[doc = concat!($description, " CDF inverse for `p` in `[0,1]`, including support endpoints.")]
        pub fn $inv(p: f64, $($arg: f64),+) -> Result<f64, Error> { $inverse(p, false, $($arg),+) }
        #[doc = concat!($description, " survival inverse for `p` in `[0,1]`, including support endpoints.")]
        pub fn $isf(p: f64, $($arg: f64),+) -> Result<f64, Error> { $inverse(p, true, $($arg),+) }
    };
}

/// Checks that gamma `shape` is positive, finite, and within the supported numerical limit.
fn gamma_shape(shape: f64) -> Result<(), Error> {
    positive(shape)?;
    if shape > 1e8 {
        Err(Error::UnsupportedRange)
    } else {
        Ok(())
    }
}

// Extend the special crate's gamma tails beyond its documented 1e6 limit.
// The deviance form of the Stirling prefactor retains precision near x=a;
// a positive series / modified Lentz fraction evaluates the smaller tail.
/// Returns the lower and upper gamma probabilities for positive `shape` and unit-scale observation `x`, or a range or convergence error.
pub(crate) fn gamma_unit_tails(shape: f64, x: f64) -> Result<(f64, f64), Error> {
    gamma_shape(shape)?;
    if shape <= 1e6 {
        return Ok((
            incomplete_gamma::lower_regularized(shape, x)?,
            incomplete_gamma::upper_regularized(shape, x)?,
        ));
    }
    if x == 0.0 {
        return Ok((0.0, 1.0));
    }
    if x.is_infinite() {
        return Ok((1.0, 0.0));
    }
    let u = (x - shape) / shape;
    let deviance = if u.abs() < 0.5 {
        let mut term = -u * u / 2.0;
        let mut sum = term;
        for k in 3..=64 {
            term *= -u * f64::from(k - 1) / f64::from(k);
            sum += term;
        }
        shape * sum
    } else {
        shape * (x.ln() - shape.ln()) + shape - x
    };
    let log_front = deviance + 0.5 * shape.ln() - 0.9189385332046727 - 1.0 / (12.0 * shape);
    let (log, lower) = if x < shape + 1.0 {
        let mut term = 1.0;
        let mut sum = 1.0;
        let mut converged = false;
        for n in 1..=262_144 {
            term *= x / (shape + f64::from(n));
            sum += term;
            let ratio = x / (shape + f64::from(n) + 1.0);
            if ratio < 1.0 && term * ratio / (1.0 - ratio) <= f64::EPSILON * sum {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(Error::NonConvergence);
        }
        (log_front - shape.ln() + sum.ln(), true)
    } else {
        let mut b = x + 1.0 - shape;
        let mut c = 1e300;
        let mut d = 1.0 / b;
        let mut h = d;
        let mut converged = false;
        for n in 1..=262_144 {
            let n = f64::from(n);
            let an = -n * (n - shape);
            b += 2.0;
            d = 1.0 / guard(b + an * d);
            c = guard(b + an / c);
            let delta = c * d;
            h *= delta;
            if (delta - 1.0).abs() <= 4.0 * f64::EPSILON {
                converged = true;
                break;
            }
        }
        if !converged || h <= 0.0 {
            return Err(Error::NonConvergence);
        }
        (log_front + h.ln(), false)
    };
    let small = log.exp();
    let large = -log.exp_m1();
    if !(0.0..=1.0).contains(&small) {
        return Err(Error::NonConvergence);
    }
    Ok(if lower {
        (small, large)
    } else {
        (large, small)
    })
}

/// Protects a continued-fraction denominator `x` from zero while preserving its sign.
fn guard(x: f64) -> f64 {
    if x.abs() < 1e-300 {
        1e-300_f64.copysign(x)
    } else {
        x
    }
}

/// Returns gamma lower and upper tails at `x` for positive `shape` and `scale`.
fn gamma_tails(x: f64, shape: f64, scale: f64) -> Result<(f64, f64), Error> {
    observation(x)?;
    gamma_shape(shape)?;
    positive(scale)?;
    let unit = x.max(0.0) / scale;
    if x > 0.0 && unit < 1e-300 {
        let log_p = shape * (x.ln() - scale.ln()) - log_gamma_one_plus(shape)?;
        return Ok((log_p.exp(), -log_p.exp_m1()));
    }
    gamma_unit_tails(shape, unit)
}

/// Returns `ln Gamma(1+shape)`, retaining small values when adding one would round away `shape`.
fn log_gamma_one_plus(shape: f64) -> Result<f64, Error> {
    if shape < 1e-8 {
        // DLMF 5.7.3, retaining the linear term when 1+shape rounds to 1.
        Ok(shape
            * (-0.5772156649015329
                + shape * (std::f64::consts::PI.powi(2) / 12.0 - shape * 0.4006856343865314)))
    } else {
        Ok(ln_abs_gamma_with_sign(shape + 1.0)?.0)
    }
}

/// Inverts the gamma tail at probability `p` for `shape, scale`; `upper` selects the survival inverse.
fn gamma_inverse(p: f64, upper: bool, shape: f64, scale: f64) -> Result<f64, Error> {
    probability(p)?;
    gamma_shape(shape)?;
    positive(scale)?;
    if p == 0.0 {
        return Ok(if upper { f64::INFINITY } else { 0.0 });
    }
    if p == 1.0 {
        return Ok(if upper { 0.0 } else { f64::INFINITY });
    }
    if shape == 1.0 {
        return finite_quantile(-scale * if upper { p.ln() } else { (-p).ln_1p() }, p);
    }
    let (target, upper) = if p > 0.5 {
        (1.0 - p, !upper)
    } else {
        (p, upper)
    };
    if !upper {
        // P(a,x) ~ x^a / Gamma(a+1). Below exp(-700), corrections are
        // negligible. Apply scale in log space before rounding: an underflowed
        // unit quantile can still give a representable scaled quantile.
        let log_x = (target.ln() + log_gamma_one_plus(shape)?) / shape;
        if log_x < -700.0 {
            return finite_quantile((log_x + scale.ln()).exp(), p);
        }
    }
    let unit = positive_quantile(target, upper, f64::MAX, |x| {
        let tails = gamma_unit_tails(shape, x)?;
        Ok(if upper { tails.1 } else { tails.0 })
    })?;
    finite_quantile(scale * unit, p)
}

family!(
    gamma_cdf,
    gamma_survival,
    gamma_quantile,
    gamma_inverse_survival,
    (shape, scale),
    "Gamma distribution with positive `shape` and `scale` (GSL a,b).",
    gamma_tails,
    gamma_inverse
);

/// Returns chi-squared lower and upper tails at `x` for positive degrees of freedom `degrees`.
fn chi_tails(x: f64, degrees: f64) -> Result<(f64, f64), Error> {
    positive(degrees)?;
    gamma_tails(x, degrees / 2.0, 2.0)
}
/// Inverts a chi-squared tail at `p` for `degrees`; `upper` selects the survival inverse.
fn chi_inverse(p: f64, upper: bool, degrees: f64) -> Result<f64, Error> {
    positive(degrees)?;
    gamma_inverse(p, upper, degrees / 2.0, 2.0)
}
family!(
    chi_squared_cdf,
    chi_squared_survival,
    chi_squared_quantile,
    chi_squared_inverse_survival,
    (degrees),
    "Chi-squared distribution with positive degrees of freedom `degrees`.",
    chi_tails,
    chi_inverse
);

/// Validates positive finite beta shapes `a, b` against the supported numerical limit.
fn beta_shapes(a: f64, b: f64) -> Result<(), Error> {
    positive(a)?;
    positive(b)?;
    if a > 1e8 || b > 1e8 {
        Err(Error::UnsupportedRange)
    } else {
        Ok(())
    }
}

/// Returns beta lower and upper tails at `x` for positive shapes `a, b`, including support limits.
fn beta_tails(x: f64, a: f64, b: f64) -> Result<(f64, f64), Error> {
    observation(x)?;
    beta_shapes(a, b)?;
    let x = x.clamp(0.0, 1.0);
    if a <= 1e6 && b <= 1e6 {
        return Ok((
            incomplete_beta::regularized(a, b, x)?,
            incomplete_beta::regularized_complement(a, b, x)?,
        ));
    }
    if x == 0.0 {
        return Ok((0.0, 1.0));
    }
    if x == 1.0 {
        return Ok((1.0, 0.0));
    }
    // GSL's large/unequal-shape gamma limit (AS 26.5.17). Besides avoiding
    // cancellation of large logs, gamma's small-shape series retains O(a)
    // complements when a is tiny and the other shape exceeds 1e6.
    if b > 1e6 && a < 10.0 {
        return gamma_unit_tails(a, -(b + 0.5 * (a - 1.0)) * (-x).ln_1p());
    }
    if a > 1e6 && b < 10.0 {
        let (q, p) = gamma_unit_tails(b, -(a + 0.5 * (b - 1.0)) * x.ln())?;
        return Ok((p, q));
    }
    let swap = x > (a + 1.0) / (a + b + 2.0);
    let (aa, bb, xx) = if swap { (b, a, 1.0 - x) } else { (a, b, x) };
    let mut c = 1.0;
    let mut d = 1.0 / guard(1.0 - (aa + bb) * xx / (aa + 1.0));
    let mut h = d;
    for m in 1..=65_536 {
        let m = f64::from(m);
        let even = m * (bb - m) * xx / ((aa + 2.0 * m - 1.0) * (aa + 2.0 * m));
        let odd = -(aa + m) * (aa + bb + m) * xx / ((aa + 2.0 * m) * (aa + 2.0 * m + 1.0));
        let mut delta = 0.0;
        for numerator in [even, odd] {
            d = 1.0 / guard(1.0 + numerator * d);
            c = guard(1.0 + numerator / c);
            delta = c * d;
            h *= delta;
        }
        if (delta - 1.0).abs() <= 4.0 * f64::EPSILON {
            let log = a * x.ln() + b * (-x).ln_1p() - ln_beta(a, b)? + h.ln() - aa.ln();
            let small = log.exp();
            let large = -log.exp_m1();
            if !(0.0..=1.0).contains(&small) || !(0.0..=1.0).contains(&large) {
                return Err(Error::NonConvergence);
            }
            return Ok(if swap { (large, small) } else { (small, large) });
        }
    }
    Err(Error::NonConvergence)
}

/// Inverts the beta tail at probability `p` for shapes `a, b`; `upper` selects the survival inverse.
fn beta_inverse(p: f64, upper: bool, a: f64, b: f64) -> Result<f64, Error> {
    probability(p)?;
    beta_shapes(a, b)?;
    if p == 0.0 {
        return Ok(if upper { 1.0 } else { 0.0 });
    }
    if p == 1.0 {
        return Ok(if upper { 0.0 } else { 1.0 });
    }
    if p > 0.5 {
        return beta_inverse(1.0 - p, !upper, a, b);
    }
    if upper {
        return Ok(1.0 - beta_inverse(p, false, b, a)?);
    }
    let log_x = (p.ln() + a.ln() + ln_beta(a, b)?) / a;
    if log_x < -700.0 {
        return Ok(log_x.exp());
    }
    positive_quantile(p, false, 1.0, |x| Ok(beta_tails(x, a, b)?.0))
}
family!(
    beta_cdf,
    beta_survival,
    beta_quantile,
    beta_inverse_survival,
    (a, b),
    "Beta distribution on `[0,1]` with positive shapes `a,b`.",
    beta_tails,
    beta_inverse
);

/// Returns Fisher F lower and upper tails at `x` for numerator and denominator degrees of freedom.
fn fisher_tails(
    x: f64,
    numerator_degrees: f64,
    denominator_degrees: f64,
) -> Result<(f64, f64), Error> {
    observation(x)?;
    beta_shapes(numerator_degrees / 2.0, denominator_degrees / 2.0)?;
    if x <= 0.0 {
        return Ok((0.0, 1.0));
    }
    let ratio = denominator_degrees / numerator_degrees;
    if x < ratio {
        beta_tails(
            x / (ratio + x),
            numerator_degrees / 2.0,
            denominator_degrees / 2.0,
        )
    } else {
        let (q, p) = beta_tails(
            ratio / (ratio + x),
            denominator_degrees / 2.0,
            numerator_degrees / 2.0,
        )?;
        Ok((p, q))
    }
}
/// Inverts a Fisher F tail at `p` for numerator and denominator degrees of freedom; `upper` selects the survival inverse.
fn fisher_inverse(
    p: f64,
    upper: bool,
    numerator_degrees: f64,
    denominator_degrees: f64,
) -> Result<f64, Error> {
    probability(p)?;
    beta_shapes(numerator_degrees / 2.0, denominator_degrees / 2.0)?;
    let (p_small, upper) = if p > 0.5 {
        (1.0 - p, !upper)
    } else {
        (p, upper)
    };
    let y = if upper {
        beta_inverse(
            p_small,
            false,
            denominator_degrees / 2.0,
            numerator_degrees / 2.0,
        )?
    } else {
        beta_inverse(
            p_small,
            false,
            numerator_degrees / 2.0,
            denominator_degrees / 2.0,
        )?
    };
    finite_quantile(
        if upper {
            (denominator_degrees / numerator_degrees) * ((1.0 - y) / y)
        } else {
            (denominator_degrees / numerator_degrees) * (y / (1.0 - y))
        },
        p,
    )
}
family!(
    fisher_f_cdf,
    fisher_f_survival,
    fisher_f_quantile,
    fisher_f_inverse_survival,
    (numerator_degrees, denominator_degrees),
    "Fisher F distribution with positive numerator and denominator degrees of freedom.",
    fisher_tails,
    fisher_inverse
);

/// Returns Student t lower and upper tails at `x` for positive degrees of freedom `degrees`.
fn student_tails(x: f64, degrees: f64) -> Result<(f64, f64), Error> {
    observation(x)?;
    beta_shapes(0.5, degrees / 2.0)?;
    let t = x.abs() / degrees.sqrt();
    let small = if t < 1.0 {
        let u = t * t;
        0.5 * beta_tails(u / (1.0 + u), 0.5, degrees / 2.0)?.1
    } else {
        let v = (1.0 / t).powi(2);
        if v < 1e-16 {
            // Preserve heavy tails even when squaring 1/t underflows.
            let log_t = if t.is_infinite() {
                x.abs().ln() - 0.5 * degrees.ln()
            } else {
                t.ln()
            };
            (-std::f64::consts::LN_2
                - degrees * log_t
                - 0.5 * degrees * v.ln_1p()
                - (degrees / 2.0).ln()
                - ln_beta(degrees / 2.0, 0.5)?)
            .exp()
        } else {
            0.5 * beta_tails(v / (1.0 + v), degrees / 2.0, 0.5)?.0
        }
    };
    Ok(if x < 0.0 {
        (small, 1.0 - small)
    } else {
        (1.0 - small, small)
    })
}
/// Inverts a Student t tail at `p` for `degrees`; `upper` selects the survival inverse.
fn student_inverse(p: f64, upper: bool, degrees: f64) -> Result<f64, Error> {
    probability(p)?;
    beta_shapes(0.5, degrees / 2.0)?;
    if p == 0.0 {
        return Ok(if upper {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        });
    }
    if p == 1.0 {
        return Ok(if upper {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        });
    }
    if p == 0.5 {
        return Ok(0.0);
    }
    let small = p.min(1.0 - p);
    let value = positive_quantile(small, true, f64::MAX, |x| Ok(student_tails(x, degrees)?.1))?;
    Ok(if (p < 0.5) != upper { -value } else { value })
}
family!(
    student_t_cdf,
    student_t_survival,
    student_t_quantile,
    student_t_inverse_survival,
    (degrees),
    "Student t distribution with positive degrees of freedom `degrees`.",
    student_tails,
    student_inverse
);

/// Returns exponential-power lower and upper tails at `x` for positive `scale` and `shape`.
fn exponential_power_tails(x: f64, scale: f64, shape: f64) -> Result<(f64, f64), Error> {
    observation(x)?;
    positive(scale)?;
    positive(shape)?;
    let (p, q) = gamma_unit_tails(1.0 / shape, (x.abs() / scale).powf(shape))?;
    Ok(if x < 0.0 {
        (0.5 * q, 0.5 + 0.5 * p)
    } else {
        (0.5 + 0.5 * p, 0.5 * q)
    })
}
/// Exponential-power CDF with positive `scale,shape` (GSL a,b); gamma shape is `1/shape`.
pub fn exponential_power_cdf(x: f64, scale: f64, shape: f64) -> Result<f64, Error> {
    Ok(exponential_power_tails(x, scale, shape)?.0)
}
/// Exponential-power survival with positive `scale,shape`, retaining small upper tails.
pub fn exponential_power_survival(x: f64, scale: f64, shape: f64) -> Result<f64, Error> {
    Ok(exponential_power_tails(x, scale, shape)?.1)
}

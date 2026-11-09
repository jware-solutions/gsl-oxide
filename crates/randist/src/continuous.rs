// Adapted from GSL 2.8 randist (GPL-3.0-or-later).
// Copyright (C) 1996-2000, 2007 James Theiler, Brian Gough.
use crate::{
    count, density, deviance, finite, log_beta, log_gamma, normal, normal_ziggurat, observation,
    open, positive, standard_normal, stirling, uniform, xlog, Error, RngCore, LOG_SQRT_2_PI,
};
use std::f64::consts::{E, FRAC_PI_2, PI};

/// Draws a gamma-distributed sample from random source `r`.
///
/// Finite positive `shape` and `scale` give mean `shape * scale` and variance
/// `shape * scale^2`. Returns a nonnegative sample using Marsaglia-Tsang;
/// invalid parameters return a domain error and overflow returns a range error.
pub fn gamma<R: RngCore + ?Sized>(r: &mut R, shape: f64, scale: f64) -> Result<f64, Error> {
    positive(shape)?;
    positive(scale)?;
    if shape < 1.0 {
        let u = open(r);
        return finite(gamma(r, shape + 1.0, scale)? * u.powf(1.0 / shape));
    }
    let d = shape - 1.0 / 3.0;
    let c = (1.0 / 3.0) / d.sqrt();
    loop {
        let x = normal_ziggurat(r, 1.0)?;
        let v = 1.0 + c * x;
        if v <= 0.0 {
            continue;
        }
        let v = v * v * v;
        let u = open(r);
        if u < 1.0 - 0.0331 * x.powi(4) || u.ln() < 0.5 * x * x + d * (1.0 - v + v.ln()) {
            return finite(scale * (d * v));
        }
    }
}

/// Draws a unit-scale gamma sample with integer `shape` using `r`.
///
/// For positive shapes, the result has mean and variance equal to `shape`.
/// Zero shape returns zero, following GSL's convention. Small shapes sum
/// exponential variates; larger shapes use Knuth's rejection algorithm.
pub fn gamma_integer<R: RngCore + ?Sized>(r: &mut R, shape: u32) -> Result<f64, Error> {
    if shape < 12 {
        return Ok((0..shape).map(|_| -open(r).ln()).sum());
    }
    gamma_large(r, f64::from(shape))
}

/// Draws a unit-scale gamma sample for a caller-validated shape `a > 1`.
///
/// Uses random source `r` and Ahrens' tangent-rejection method as described
/// by Knuth. Returns a positive sample or `Error::NumericalRange` if the
/// proposal scale or accepted result cannot be represented.
fn gamma_large<R: RngCore + ?Sized>(r: &mut R, a: f64) -> Result<f64, Error> {
    let s = (2.0 * a - 1.0).sqrt();
    if !s.is_finite() {
        return Err(Error::NumericalRange);
    }
    loop {
        let y = (PI * uniform(r)).tan();
        let x = s * y + a - 1.0;
        if x <= 0.0 {
            continue;
        }
        if uniform(r)
            <= (1.0 + y * y) * ((a - 1.0) * ((x - (a - 1.0)) / (a - 1.0)).ln_1p() - s * y).exp()
        {
            return finite(x);
        }
    }
}

/// Draws a unit-scale gamma sample for a caller-validated `0 <= a < 1`.
///
/// Uses random source `r` and Knuth's fractional-shape rejection method.
/// Returns a nonnegative sample; zero shape returns zero and very small
/// positive samples can underflow to zero.
fn gamma_fraction<R: RngCore + ?Sized>(r: &mut R, a: f64) -> f64 {
    if a == 0.0 {
        return 0.0;
    }
    let p = E / (a + E);
    loop {
        let u = uniform(r);
        let v = open(r);
        let (x, q) = if u < p {
            let x = v.powf(1.0 / a);
            (x, (-x).exp())
        } else {
            let x = 1.0 - v.ln();
            (x, x.powf(a - 1.0))
        };
        if uniform(r) < q {
            return x;
        }
    }
}

/// Draws a gamma sample using Knuth's integer/fraction decomposition.
///
/// `r` supplies randomness; `shape` and `scale` must be finite and positive.
/// Returns a nonnegative sample with mean `shape * scale` and variance
/// `shape * scale^2`, or an error for invalid parameters or numerical overflow.
pub fn gamma_knuth<R: RngCore + ?Sized>(r: &mut R, shape: f64, scale: f64) -> Result<f64, Error> {
    positive(shape)?;
    positive(scale)?;
    let n = shape.floor();
    let whole = if n >= f64::from(u32::MAX) {
        gamma_large(r, n)?
    } else {
        gamma_integer(r, n as u32)?
    };
    finite(scale * (whole + gamma_fraction(r, shape - n)))
}

/// Evaluates the gamma density at observation `x`.
///
/// Finite positive `shape` and `scale` specify the distribution as in [`gamma`].
/// Returns zero for negative or infinite `x`. At zero the limiting density is
/// infinite for `shape < 1`, `1/scale` for `shape = 1`, and zero for `shape > 1`.
/// NaN observations and invalid parameters return a domain error.
pub fn gamma_density(x: f64, shape: f64, scale: f64) -> Result<f64, Error> {
    positive(shape)?;
    positive(scale)?;
    observation(x)?;
    if x < 0.0 || x.is_infinite() {
        return Ok(0.0);
    }
    if x == 0.0 {
        return Ok(if shape < 1.0 {
            f64::INFINITY
        } else if shape == 1.0 {
            1.0 / scale
        } else {
            0.0
        });
    }
    let t = x / scale;
    if t.is_infinite() {
        return Ok(0.0);
    }
    let log = if shape >= 16.0 && t > 0.0 {
        -deviance(shape, t) - t.ln() + 0.5 * shape.ln()
            - LOG_SQRT_2_PI
            - stirling(shape)
            - scale.ln()
    } else {
        (shape - 1.0) * (x.ln() - scale.ln()) - t - log_gamma(shape)? - scale.ln()
    };
    density(log.exp())
}

/// Draws a beta-distributed sample on `[0, 1]` from random source `r`.
///
/// `a` and `b` are finite positive shape parameters; the mean is `a/(a+b)`.
/// Returns a sample from a density proportional to `x^(a-1) * (1-x)^(b-1)`.
/// Roundoff or underflow can produce either endpoint for small shapes.
pub fn beta<R: RngCore + ?Sized>(r: &mut R, a: f64, b: f64) -> Result<f64, Error> {
    positive(a)?;
    positive(b)?;
    if a <= 1.0 && b <= 1.0 {
        loop {
            let lx = open(r).ln() / a;
            let ly = open(r).ln() / b;
            if lx.exp() + ly.exp() <= 1.0 {
                return finite(if lx > ly {
                    1.0 / (1.0 + (ly - lx).exp())
                } else {
                    let t = (lx - ly).exp();
                    t / (1.0 + t)
                });
            }
        }
    }
    let x = gamma(r, a, 1.0)?;
    let y = gamma(r, b, 1.0)?;
    if x > y {
        Ok(1.0 / (1.0 + y / x))
    } else {
        finite((x / y) / (1.0 + x / y))
    }
}

/// Evaluates a beta density with shape parameters `a` and `b` at `x`.
///
/// Both shapes must be finite and positive. Returns
/// `x^(a-1) * (1-x)^(b-1) / Beta(a,b)` on `[0,1]`, zero outside that interval,
/// and the one-sided limiting density at either endpoint, possibly infinity.
/// NaN observations or invalid shapes return a domain error.
pub fn beta_density(x: f64, a: f64, b: f64) -> Result<f64, Error> {
    positive(a)?;
    positive(b)?;
    observation(x)?;
    if !(0.0..=1.0).contains(&x) {
        return Ok(0.0);
    }
    density(
        (xlog(a - 1.0, x)
            + if b == 1.0 {
                0.0
            } else {
                (b - 1.0) * (-x).ln_1p()
            }
            - log_beta(a, b)?)
        .exp(),
    )
}

/// Draws an exponential sample with positive mean `mean` using `r`.
///
/// `mean` must be finite and is also the scale parameter. Returns a
/// nonnegative sample with variance `mean^2`; invalid means and overflowing
/// samples return domain and range errors, respectively.
pub fn exponential<R: RngCore + ?Sized>(r: &mut R, mean: f64) -> Result<f64, Error> {
    positive(mean)?;
    finite(-mean * open(r).ln())
}
/// Evaluates the exponential density at observation `x`.
///
/// Finite positive `mean` is the distribution's mean and scale. Returns
/// `exp(-x/mean)/mean` for nonnegative `x`, and zero for negative or infinite
/// observations. NaN or invalid `mean` returns a domain error.
pub fn exponential_density(x: f64, mean: f64) -> Result<f64, Error> {
    positive(mean)?;
    observation(x)?;
    if x < 0.0 {
        Ok(0.0)
    } else {
        density((-x / mean - mean.ln()).exp())
    }
}
/// Draws a symmetric exponential-power sample using random source `r`.
///
/// Finite positive `a` is the scale and finite positive `b` is the exponent
/// controlling tail shape. Returns a real sample centered at zero: `b = 1`
/// gives Laplace scale `a`, while `b = 2` gives normal deviation `a/sqrt(2)`.
pub fn exponential_power<R: RngCore + ?Sized>(r: &mut R, a: f64, b: f64) -> Result<f64, Error> {
    positive(a)?;
    positive(b)?;
    let sign = if uniform(r) < 0.5 { -1.0 } else { 1.0 };
    finite(sign * a * gamma(r, 1.0 / b, 1.0)?.powf(1.0 / b))
}
/// Evaluates the symmetric exponential-power density at observation `x`.
///
/// `a` is a finite positive scale and `b` a finite positive shape exponent.
/// Returns `b * exp(-(|x|/a)^b) / (2*a*Gamma(1/b))` on the real line.
/// Infinite `x` gives zero; NaN or invalid parameters gives a domain error.
pub fn exponential_power_density(x: f64, a: f64, b: f64) -> Result<f64, Error> {
    positive(a)?;
    positive(b)?;
    observation(x)?;
    density((b.ln() - 2.0_f64.ln() - a.ln() - log_gamma(1.0 / b)? - (x.abs() / a).powf(b)).exp())
}
/// Draws a zero-location Cauchy sample from random source `r`.
///
/// Finite positive `scale` controls the width of the distribution. Returns
/// one real sample; this distribution has no finite mean or variance.
/// Invalid scale or an overflowing sample returns an error.
pub fn cauchy<R: RngCore + ?Sized>(r: &mut R, scale: f64) -> Result<f64, Error> {
    positive(scale)?;
    finite(scale * (PI * (open(r) - 0.5)).tan())
}
/// Evaluates the zero-location Cauchy density at observation `x`.
///
/// For finite positive `scale`, returns
/// `scale / (pi * (x^2 + scale^2))`. Infinite observations return zero;
/// NaN observations or invalid scale return a domain error.
pub fn cauchy_density(x: f64, scale: f64) -> Result<f64, Error> {
    positive(scale)?;
    observation(x)?;
    // hypot preserves finite tails when squaring would overflow.
    density((-PI.ln() + scale.ln() - 2.0 * x.hypot(scale).ln()).exp())
}
/// Draws a chi-squared sample with `degrees` degrees of freedom using `r`.
///
/// `degrees` must be finite and positive, and may be noninteger. Returns
/// a nonnegative sample with mean `degrees` and variance `2*degrees`,
/// equivalently a gamma sample with shape `degrees/2` and scale two.
pub fn chi_squared<R: RngCore + ?Sized>(r: &mut R, degrees: f64) -> Result<f64, Error> {
    gamma(r, degrees / 2.0, 2.0)
}
/// Evaluates the chi-squared density at observation `x`.
///
/// `degrees` is a finite positive, possibly noninteger number of degrees
/// of freedom. Returns the gamma density with shape `degrees/2` and scale
/// two, including its limiting value at zero and zero for negative `x`.
pub fn chi_squared_density(x: f64, degrees: f64) -> Result<f64, Error> {
    gamma_density(x, degrees / 2.0, 2.0)
}
/// Draws an Erlang sample with `scale` and `shape` using random source `r`.
///
/// Both parameters must be finite and positive. Following GSL, `shape`
/// may be noninteger. Returns a nonnegative gamma sample with mean
/// `shape * scale`; the argument order places scale before shape.
pub fn erlang<R: RngCore + ?Sized>(r: &mut R, scale: f64, shape: f64) -> Result<f64, Error> {
    gamma(r, shape, scale)
}
/// Evaluates the Erlang density at observation `x`.
///
/// Finite positive `scale` precedes finite positive `shape`; GSL permits
/// noninteger shapes. Returns the gamma density with these parameters,
/// including its limiting value at zero and zero for negative `x`.
pub fn erlang_density(x: f64, scale: f64, shape: f64) -> Result<f64, Error> {
    gamma_density(x, shape, scale)
}
/// Draws a Fisher F sample from random source `r`.
///
/// `numerator` and `denominator` are finite positive degrees of freedom.
/// Returns a nonnegative ratio of two independent chi-squared samples,
/// each divided by its degrees of freedom; overflow returns a range error.
pub fn fisher_f<R: RngCore + ?Sized>(
    r: &mut R,
    numerator: f64,
    denominator: f64,
) -> Result<f64, Error> {
    positive(numerator)?;
    positive(denominator)?;
    finite((chi_squared(r, numerator)? / numerator) / (chi_squared(r, denominator)? / denominator))
}
/// Evaluates the Fisher F density at observation `x`.
///
/// `numerator` and `denominator` are finite positive degrees of freedom.
/// Returns zero for negative or infinite `x`, and the limiting density at
/// zero (infinite, one, or zero as `numerator` is below, equal to, or above two).
/// Logarithmic evaluation avoids the intermediate overflow in GSL bug #28500.
pub fn fisher_f_density(x: f64, numerator: f64, denominator: f64) -> Result<f64, Error> {
    positive(numerator)?;
    positive(denominator)?;
    observation(x)?;
    if x < 0.0 || x.is_infinite() {
        return Ok(0.0);
    }
    if x == 0.0 {
        return Ok(if numerator < 2.0 {
            f64::INFINITY
        } else if numerator == 2.0 {
            1.0
        } else {
            0.0
        });
    }
    let a = numerator / 2.0;
    let b = denominator / 2.0;
    let l = numerator.ln() - denominator.ln() + x.ln();
    let log = if l > 0.0 {
        -b * l - (a + b) * (-l).exp().ln_1p()
    } else {
        a * l - (a + b) * l.exp().ln_1p()
    };
    density((log - x.ln() - log_beta(a, b)?).exp())
}
/// Draws uniformly from `[lower, upper)` using random source `r`.
///
/// Both bounds must be finite and satisfy `lower < upper`. Returns one
/// value including the lower endpoint and excluding the upper endpoint,
/// even when rounding would otherwise produce the upper bound.
pub fn uniform_interval<R: RngCore + ?Sized>(
    r: &mut R,
    lower: f64,
    upper: f64,
) -> Result<f64, Error> {
    interval(lower, upper)?;
    let u = uniform(r);
    let x = lower * (1.0 - u) + upper * u;
    // Rounding may otherwise return the excluded upper endpoint.
    Ok(x.max(lower).min(previous(upper)))
}
/// Returns the immediately preceding representable value below finite `x`.
///
/// The caller supplies a finite interval endpoint. Either signed zero maps
/// to the smallest negative subnormal, preserving an excluded upper bound.
fn previous(x: f64) -> f64 {
    if x == 0.0 {
        -f64::from_bits(1)
    } else {
        f64::from_bits(if x > 0.0 {
            x.to_bits() - 1
        } else {
            x.to_bits() + 1
        })
    }
}
/// Validates finite interval endpoints `a` and `b` with `a < b`.
///
/// Returns `Ok(())` for ordered finite bounds and `Error::Domain` otherwise;
/// the width itself need not fit in a finite `f64`.
fn interval(a: f64, b: f64) -> Result<(), Error> {
    if a.is_finite() && b.is_finite() && a < b {
        Ok(())
    } else {
        Err(Error::Domain)
    }
}
/// Evaluates the uniform density on `[lower, upper)` at observation `x`.
///
/// Finite bounds must satisfy `lower < upper`. Returns `1/(upper-lower)`
/// inside the interval and zero outside it, including at the upper endpoint.
/// NaN observations and invalid bounds return a domain error.
pub fn uniform_interval_density(x: f64, lower: f64, upper: f64) -> Result<f64, Error> {
    interval(lower, upper)?;
    observation(x)?;
    if x < lower || x >= upper {
        Ok(0.0)
    } else {
        let width = upper - lower;
        if width.is_infinite() {
            Ok(0.5 / (0.5 * upper - 0.5 * lower))
        } else {
            Ok(1.0 / width)
        }
    }
}
/// Draws a type-I Gumbel sample with GSL parameters `a` and `b` using `r`.
///
/// Both parameters must be finite and positive. `a` is an inverse scale;
/// `b` shifts the location by `ln(b)/a`. Returns a real sample whose CDF is
/// `exp(-b*exp(-a*x))`; invalid parameters or overflow returns an error.
pub fn gumbel_type_1<R: RngCore + ?Sized>(r: &mut R, a: f64, b: f64) -> Result<f64, Error> {
    positive(a)?;
    positive(b)?;
    finite((b.ln() - (-open(r).ln()).ln()) / a)
}
/// Evaluates the type-I Gumbel density at observation `x`.
///
/// Finite positive `a` is the inverse scale and finite positive `b`
/// sets location `ln(b)/a`. Returns `a*b*exp(-a*x-b*exp(-a*x))` on the
/// real line, or zero at either infinity. NaN returns a domain error.
pub fn gumbel_type_1_density(x: f64, a: f64, b: f64) -> Result<f64, Error> {
    positive(a)?;
    positive(b)?;
    observation(x)?;
    if x.is_infinite() {
        return Ok(0.0);
    }
    let z = b.ln() - a * x;
    if z == f64::INFINITY {
        Ok(0.0)
    } else {
        density((a.ln() + z - z.exp()).exp())
    }
}
/// Draws a type-II Gumbel (Frechet) sample with GSL parameters using `r`.
///
/// Finite positive `a` is the shape; finite positive `b` gives scale
/// `b^(1/a)`. Returns a nonnegative sample with CDF `exp(-b*x^(-a))`
/// for positive `x`; underflow may round a small sample to zero.
pub fn gumbel_type_2<R: RngCore + ?Sized>(r: &mut R, a: f64, b: f64) -> Result<f64, Error> {
    positive(a)?;
    positive(b)?;
    finite(((b.ln() - (-open(r).ln()).ln()) / a).exp())
}
/// Evaluates the type-II Gumbel (Frechet) density at observation `x`.
///
/// Finite positive `a` is the shape and finite positive `b` gives scale
/// `b^(1/a)`. Returns `a*b*x^(-a-1)*exp(-b*x^(-a))` for positive finite `x`,
/// and zero for nonpositive or infinite `x`. NaN returns a domain error.
pub fn gumbel_type_2_density(x: f64, a: f64, b: f64) -> Result<f64, Error> {
    positive(a)?;
    positive(b)?;
    observation(x)?;
    if x <= 0.0 || x.is_infinite() {
        return Ok(0.0);
    }
    let z = b.ln() - a * x.ln();
    density((a.ln() - x.ln() + z - z.exp()).exp())
}
/// Draws a zero-location logistic sample from random source `r`.
///
/// Finite positive `scale` controls spread and gives standard deviation
/// `pi*scale/sqrt(3)`. Returns a real sample; invalid scale and sample
/// overflow return domain and range errors, respectively.
pub fn logistic<R: RngCore + ?Sized>(r: &mut R, scale: f64) -> Result<f64, Error> {
    positive(scale)?;
    let u = open(r);
    finite(scale * (u.ln() - (-u).ln_1p()))
}
/// Evaluates a zero-location logistic density at observation `x`.
///
/// `scale` must be finite and positive. Returns the symmetric density
/// `exp(-|x|/scale) / (scale*(1+exp(-|x|/scale))^2)`, with zero at infinity.
/// NaN observations or invalid scale return a domain error.
pub fn logistic_density(x: f64, scale: f64) -> Result<f64, Error> {
    positive(scale)?;
    observation(x)?;
    let z = -x.abs() / scale;
    density((z - 2.0 * z.exp().ln_1p() - scale.ln()).exp())
}
/// Draws a lognormal sample using random source `r`.
///
/// Finite `location` is the mean of the logarithm, and finite positive
/// `sigma` is its standard deviation. Returns `exp(location + sigma*Z)`
/// for standard normal `Z`; overflow is an error and underflow may give zero.
pub fn lognormal<R: RngCore + ?Sized>(r: &mut R, location: f64, sigma: f64) -> Result<f64, Error> {
    if !location.is_finite() {
        return Err(Error::Domain);
    }
    finite((location + normal(r, sigma)?).exp())
}
/// Evaluates the lognormal density at observation `x`.
///
/// Finite `location` and finite positive `sigma` are the mean and standard
/// deviation of `ln(X)`. Returns the corresponding density for positive `x`,
/// and zero for nonpositive or infinite `x`; NaN returns a domain error.
pub fn lognormal_density(x: f64, location: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    observation(x)?;
    if !location.is_finite() {
        return Err(Error::Domain);
    }
    if x <= 0.0 || x.is_infinite() {
        return Ok(0.0);
    }
    density(
        (-0.5 * ((x.ln() - location) / sigma).powi(2) - LOG_SQRT_2_PI - sigma.ln() - x.ln()).exp(),
    )
}
/// Draws a Pareto sample from random source `r`.
///
/// Finite positive `shape` controls the power-law tail and finite positive
/// `minimum` is the lower support bound. Returns a sample on
/// `[minimum, infinity)`, or a range error if a tail sample overflows.
pub fn pareto<R: RngCore + ?Sized>(r: &mut R, shape: f64, minimum: f64) -> Result<f64, Error> {
    positive(shape)?;
    positive(minimum)?;
    finite((minimum.ln() - open(r).ln() / shape).exp())
}
/// Evaluates the Pareto density at observation `x`.
///
/// `shape` and the lower bound `minimum` must be finite and positive.
/// Returns `shape * minimum^shape / x^(shape+1)` for `x >= minimum`,
/// and zero below that bound or at infinity. NaN returns a domain error.
pub fn pareto_density(x: f64, shape: f64, minimum: f64) -> Result<f64, Error> {
    positive(shape)?;
    positive(minimum)?;
    observation(x)?;
    if x < minimum || x.is_infinite() {
        Ok(0.0)
    } else {
        density((shape.ln() - x.ln() + shape * (minimum.ln() - x.ln())).exp())
    }
}
/// Draws a Rayleigh sample with finite positive scale `sigma` using `r`.
///
/// Returns a nonnegative radius distributed like the norm of two independent
/// zero-mean normals with standard deviation `sigma`. Invalid scales and
/// sample overflow return domain and range errors, respectively.
pub fn rayleigh<R: RngCore + ?Sized>(r: &mut R, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    finite(sigma * (-2.0 * open(r).ln()).sqrt())
}
/// Evaluates the Rayleigh density at observation `x`.
///
/// Finite positive `sigma` is the scale of the underlying normal components.
/// Returns `x*exp(-x^2/(2*sigma^2))/sigma^2` for positive finite `x`, and zero
/// for nonpositive or infinite `x`; NaN returns a domain error.
pub fn rayleigh_density(x: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    observation(x)?;
    if x <= 0.0 || x.is_infinite() {
        Ok(0.0)
    } else {
        density((x.ln() - 2.0 * sigma.ln() - 0.5 * (x / sigma).powi(2)).exp())
    }
}
/// Draws a Rayleigh sample conditioned to be at least `lower` using `r`.
///
/// `lower` is a finite nonnegative cutoff and `sigma` is a finite positive
/// Rayleigh scale. Returns a sample in `[lower, infinity)`, or an error for
/// invalid parameters or an overflowing result.
pub fn rayleigh_tail<R: RngCore + ?Sized>(r: &mut R, lower: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    if lower < 0.0 || !lower.is_finite() {
        return Err(Error::Domain);
    }
    finite(lower.hypot(sigma * (-2.0 * open(r).ln()).sqrt()))
}
/// Evaluates the Rayleigh density conditioned to be at least `lower`.
///
/// `x` is the observation, `lower` a finite nonnegative cutoff, and `sigma`
/// a finite positive scale. Returns the normalized conditional density on
/// `[lower, infinity)` and zero below the cutoff or at infinity; at zero
/// the density is zero. Invalid parameters or NaN return a domain error.
pub fn rayleigh_tail_density(x: f64, lower: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    observation(x)?;
    if lower < 0.0 || !lower.is_finite() {
        return Err(Error::Domain);
    }
    if x < lower || x == 0.0 || x.is_infinite() {
        Ok(0.0)
    } else {
        density(
            (x.ln()
                - 2.0 * sigma.ln()
                - 0.5 * ((x - lower) / sigma) * ((x / sigma) + (lower / sigma)))
                .exp(),
        )
    }
}
/// Draws a standard Student t sample from random source `r`.
///
/// Finite positive `degrees` is the possibly noninteger number of degrees
/// of freedom. Returns a real sample with location zero and unit scale,
/// using a standard normal divided by a scaled chi-squared square root.
pub fn student_t<R: RngCore + ?Sized>(r: &mut R, degrees: f64) -> Result<f64, Error> {
    positive(degrees)?;
    finite(standard_normal(r) / (chi_squared(r, degrees)? / degrees).sqrt())
}
/// Evaluates the standard Student t density at observation `x`.
///
/// Finite positive `degrees` sets the degrees of freedom; location is zero
/// and scale is one. Returns a nonnegative density on the real line, zero
/// at infinity, and a domain error for NaN or invalid degrees of freedom.
pub fn student_t_density(x: f64, degrees: f64) -> Result<f64, Error> {
    positive(degrees)?;
    observation(x)?;
    if x.is_infinite() {
        return Ok(0.0);
    }
    let z = x.abs() / degrees.sqrt();
    let l = if z > 1.0 {
        2.0 * z.ln() + (1.0 / z).powi(2).ln_1p()
    } else {
        (z * z).ln_1p()
    };
    density((-0.5 * degrees.ln() - log_beta(0.5, degrees / 2.0)? - 0.5 * (degrees + 1.0) * l).exp())
}
/// Draws a zero-location Laplace sample from random source `r`.
///
/// Finite positive `scale` gives standard deviation `sqrt(2)*scale`.
/// Returns a real sample with symmetric exponential tails; invalid scale
/// or sample overflow returns an error.
pub fn laplace<R: RngCore + ?Sized>(r: &mut R, scale: f64) -> Result<f64, Error> {
    positive(scale)?;
    let sign = if uniform(r) < 0.5 { -1.0 } else { 1.0 };
    finite(sign * scale * open(r).ln())
}
/// Evaluates a zero-location Laplace density at observation `x`.
///
/// For finite positive `scale`, returns `exp(-|x|/scale)/(2*scale)`.
/// Infinite observations return zero; NaN observations and invalid scales
/// return a domain error.
pub fn laplace_density(x: f64, scale: f64) -> Result<f64, Error> {
    positive(scale)?;
    observation(x)?;
    density((-x.abs() / scale - 2.0_f64.ln() - scale.ln()).exp())
}
/// Draws a Weibull sample from random source `r`.
///
/// Finite positive `scale` and `shape` define CDF
/// `1-exp(-(x/scale)^shape)` for nonnegative `x`. Returns a nonnegative
/// sample; invalid parameters or sample overflow returns an error.
pub fn weibull<R: RngCore + ?Sized>(r: &mut R, scale: f64, shape: f64) -> Result<f64, Error> {
    positive(scale)?;
    positive(shape)?;
    finite(scale * (-open(r).ln()).powf(1.0 / shape))
}
/// Evaluates the Weibull density at observation `x`.
///
/// `scale` and `shape` must be finite and positive. Returns zero for negative
/// or infinite `x`. At zero the limiting density is infinite for `shape < 1`,
/// `1/scale` for `shape = 1`, and zero for `shape > 1`; NaN is a domain error.
pub fn weibull_density(x: f64, scale: f64, shape: f64) -> Result<f64, Error> {
    positive(scale)?;
    positive(shape)?;
    observation(x)?;
    if x < 0.0 || x.is_infinite() {
        return Ok(0.0);
    }
    if x == 0.0 {
        return Ok(if shape < 1.0 {
            f64::INFINITY
        } else if shape == 1.0 {
            1.0 / scale
        } else {
            0.0
        });
    }
    let l = x.ln() - scale.ln();
    density((shape.ln() - scale.ln() + (shape - 1.0) * l - (shape * l).exp()).exp())
}
/// Draws a symmetric Levy stable sample from random source `r`.
///
/// Finite positive `c` is the scale and `0 < alpha <= 2` is the stability
/// index. Returns a real sample with characteristic function
/// `exp(-|c*t|^alpha)`: `alpha = 1` gives Cauchy scale `c`, and `alpha = 2`
/// gives normal deviation `sqrt(2)*c`. Unrepresentable samples are errors.
pub fn stable_symmetric<R: RngCore + ?Sized>(r: &mut R, c: f64, alpha: f64) -> Result<f64, Error> {
    stable_skewed(r, c, alpha, 0.0)
}
/// Draws a skewed Levy stable sample from random source `r`.
///
/// Finite positive `c` is the scale, `0 < alpha <= 2` the stability index,
/// and `-1 <= skew <= 1` the skewness parameter. Returns a real sample in
/// GSL's convention; at `alpha = 1`, this includes the scale-dependent
/// location shift `2*skew*c*ln(c)/pi`. Zero skew gives [`stable_symmetric`].
pub fn stable_skewed<R: RngCore + ?Sized>(
    r: &mut R,
    c: f64,
    alpha: f64,
    skew: f64,
) -> Result<f64, Error> {
    positive(c)?;
    if !(0.0 < alpha && alpha <= 2.0 && (-1.0..=1.0).contains(&skew)) {
        return Err(Error::Domain);
    }
    let v = PI * (open(r) - 0.5);
    let w = -open(r).ln();
    let x = if alpha == 2.0 {
        2.0 * v.sin() * w.sqrt()
    } else if alpha == 1.0 {
        if skew == 0.0 {
            v.tan()
        } else {
            ((FRAC_PI_2 + skew * v) * v.tan()
                - skew * (FRAC_PI_2 * w * v.cos() / (FRAC_PI_2 + skew * v)).ln()
                + skew * c.ln())
                / FRAC_PI_2
        }
    } else {
        let t = skew * (FRAC_PI_2 * alpha).tan();
        let b = t.atan() / alpha;
        (1.0 + t * t).powf(0.5 / alpha) * (alpha * (v + b)).sin() / v.cos().powf(1.0 / alpha)
            * ((v - alpha * (v + b)).cos() / w).powf((1.0 - alpha) / alpha)
    };
    finite(c * x)
}

/// Floors inverse-transform value `x` and converts it to a `u32` count.
///
/// Returns zero when the floored value is zero. Nonfinite values or floored
/// values outside `[0, u32::MAX]` return `Error::NumericalRange`.
pub(crate) fn positive_count(x: f64) -> Result<u32, Error> {
    count(x.floor())
}

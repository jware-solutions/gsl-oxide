// Adapted from GSL 2.8 gauss.c, gausstail.c, gausszig.c (GPL-3.0-or-later).
// Copyright (C) 1996-2000, 2006, 2007 James Theiler, Brian Gough.
// Copyright (C) 2006 Charles Karney; (C) 2005 Jochen Voss.
use crate::{density, finite, observation, open, positive, uniform, Error, RngCore, LOG_SQRT_2_PI};
use std::f64::consts::{PI, SQRT_2};

/// GSL's precomputed Ziggurat heights, widths, and integer acceptance thresholds.
mod tables;

/// Draws a zero-mean normal sample using random source `r` and polar rejection.
///
/// `sigma` is the finite positive standard deviation. Returns one real sample
/// from `N(0, sigma^2)`, or a domain/range error if its parameters or result
/// cannot be represented.
pub fn normal<R: RngCore + ?Sized>(r: &mut R, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    finite(sigma * standard_normal(r))
}

/// Draws a standard normal sample from random source `r`.
///
/// Returns a real value with mean zero and variance one, using the polar
/// Box-Muller method. No cached second sample or shared state is retained.
pub fn standard_normal<R: RngCore + ?Sized>(r: &mut R) -> f64 {
    loop {
        let x = 2.0 * open(r) - 1.0;
        let y = 2.0 * open(r) - 1.0;
        let s = x * x + y * y;
        if s > 0.0 && s <= 1.0 {
            return y * (-2.0 * s.ln() / s).sqrt();
        }
    }
}

/// Draws a zero-mean normal sample using Leva's ratio-of-uniforms method.
///
/// `r` supplies the randomness and finite positive `sigma` sets the standard
/// deviation. Returns one sample from `N(0, sigma^2)` with the same checked
/// parameter and result contract as [`normal`].
pub fn normal_ratio<R: RngCore + ?Sized>(r: &mut R, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    finite(sigma * standard_normal_ratio(r))
}

/// Draws a standard normal sample using Leva's ratio-of-uniforms method.
///
/// Uses random source `r` and returns one real value with mean zero and
/// variance one. This is the unit-deviation counterpart of [`normal_ratio`].
pub fn standard_normal_ratio<R: RngCore + ?Sized>(r: &mut R) -> f64 {
    loop {
        let u = 1.0 - uniform(r);
        let v = 1.7156 * (uniform(r) - 0.5);
        let x = u - 0.449871;
        let y = v.abs() + 0.386595;
        let q = x * x + y * (0.196 * y - 0.25472 * x);
        if q < 0.27597 || (q <= 0.27846 && v * v <= -4.0 * u * u * u.ln()) {
            return v / u;
        }
    }
}

/// Draws a zero-mean normal sample using GSL's 128-strip Ziggurat method.
///
/// `r` supplies random bits and finite positive `sigma` sets the standard
/// deviation. Returns one sample from `N(0, sigma^2)`, reporting invalid
/// parameters or an overflowing scaled sample as an error.
pub fn normal_ziggurat<R: RngCore + ?Sized>(r: &mut R, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    const RIGHT: f64 = 3.44428647676;
    loop {
        let bits = r.next_u32();
        let i = (bits & 127) as usize;
        let sign = if bits & 128 == 0 { -1.0 } else { 1.0 };
        let j = bits >> 8;
        let mut x = f64::from(j) * tables::WTAB[i];
        if j < tables::KTAB[i] {
            return finite(sign * sigma * x);
        }
        let y = if i < 127 {
            tables::YTAB[i + 1] + (tables::YTAB[i] - tables::YTAB[i + 1]) * uniform(r)
        } else {
            x = RIGHT - (1.0 - uniform(r)).ln() / RIGHT;
            (-RIGHT * (x - 0.5 * RIGHT)).exp() * uniform(r)
        };
        if y < (-0.5 * x * x).exp() {
            return finite(sign * sigma * x);
        }
    }
}

/// Evaluates the zero-mean normal density at observation `x`.
///
/// `sigma` is the finite positive standard deviation. Returns
/// `exp(-(x/sigma)^2/2) / (sqrt(2*pi)*sigma)`. Infinite observations return
/// zero; NaN or invalid `sigma` returns `Error::Domain`.
pub fn normal_density(x: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    observation(x)?;
    density((-0.5 * (x / sigma).powi(2) - LOG_SQRT_2_PI - sigma.ln()).exp())
}

/// Evaluates the standard normal density at observation `x`.
///
/// Returns `exp(-x^2/2) / sqrt(2*pi)`, with mean zero and variance one.
/// Infinite `x` returns zero and NaN returns `Error::Domain`.
pub fn standard_normal_density(x: f64) -> Result<f64, Error> {
    normal_density(x, 1.0)
}

/// Draws a zero-mean normal sample conditioned to be at least `lower`.
///
/// `r` supplies randomness, `lower` is a finite cutoff, and `sigma` is a
/// finite positive standard deviation. Returns a sample in `[lower, infinity)`;
/// rejection sampling also supports cutoffs far into the upper tail.
pub fn normal_tail<R: RngCore + ?Sized>(r: &mut R, lower: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    if !lower.is_finite() {
        return Err(Error::Domain);
    }
    let a = lower / sigma;
    if !a.is_finite() {
        return Err(Error::NumericalRange);
    }
    if a < 1.0 {
        loop {
            let x = standard_normal(r);
            if x >= a {
                return finite(x * sigma);
            }
        }
    }
    loop {
        // hypot avoids overflowing a*a for large cutoffs.
        let x = a.hypot((-2.0 * open(r).ln()).sqrt());
        if x * uniform(r) <= a {
            return finite(x * sigma);
        }
    }
}

/// Draws a standard normal sample conditioned to be at least `lower`.
///
/// Uses random source `r` and a finite cutoff `lower`. Returns a real sample
/// in `[lower, infinity)` with the unconditioned standard deviation set to one.
pub fn standard_normal_tail<R: RngCore + ?Sized>(r: &mut R, lower: f64) -> Result<f64, Error> {
    normal_tail(r, lower, 1.0)
}

/// Evaluates a zero-mean normal density conditioned to be at least `lower`.
///
/// `x` is the observation, `lower` is a finite cutoff, and `sigma` is the
/// finite positive standard deviation before conditioning. Returns zero below
/// the cutoff and at infinity; otherwise returns the normal density divided
/// by its upper-tail probability. A Mills-ratio expansion preserves far tails.
pub fn normal_tail_density(x: f64, lower: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    observation(x)?;
    if !lower.is_finite() {
        return Err(Error::Domain);
    }
    if x < lower || x.is_infinite() {
        return Ok(0.0);
    }
    let a = lower / sigma;
    let z = x / sigma;
    if !a.is_finite() || !z.is_finite() {
        return Err(Error::NumericalRange);
    }
    if a > 8.0 {
        // Mills ratio Q(a)/phi(a) = (1/a)(1 - 1/a² + 3/a⁴ - ...).
        let inv = 1.0 / a;
        let mut sum = 1.0;
        let mut term: f64 = 1.0;
        for k in 1..=100 {
            let next = -term * f64::from(2 * k - 1) * inv * inv;
            if next.abs() >= term.abs() {
                break;
            }
            sum += next;
            term = next;
            if term.abs() < f64::EPSILON * sum {
                break;
            }
        }
        density((-0.5 * (z - a) * (z + a) + a.ln() - sum.ln() - sigma.ln()).exp())
    } else {
        let q = 0.5 * gsl_oxide_special::error_functions::erfc(a / SQRT_2);
        density((-0.5 * z * z - LOG_SQRT_2_PI - sigma.ln() - q.ln()).exp())
    }
}

/// Evaluates a standard normal density conditioned to be at least `lower`.
///
/// `x` is the observation and `lower` is a finite cutoff. Returns zero for
/// `x < lower` or infinite `x`, and the normalized conditional density otherwise.
/// NaN observations and invalid cutoffs return `Error::Domain`.
pub fn standard_normal_tail_density(x: f64, lower: f64) -> Result<f64, Error> {
    normal_tail_density(x, lower, 1.0)
}

/// Draws a correlated pair of zero-mean normal samples from `r`.
///
/// Finite positive `sigma_x` and `sigma_y` are the component standard
/// deviations, and `rho` is their correlation in `[-1, 1]`. Returns `(x, y)`;
/// endpoint correlations give perfectly correlated or anticorrelated samples.
pub fn bivariate_normal<R: RngCore + ?Sized>(
    r: &mut R,
    sigma_x: f64,
    sigma_y: f64,
    rho: f64,
) -> Result<(f64, f64), Error> {
    positive(sigma_x)?;
    positive(sigma_y)?;
    if !(-1.0..=1.0).contains(&rho) {
        return Err(Error::Domain);
    }
    let x = standard_normal(r);
    let y = rho * x + ((1.0 - rho) * (1.0 + rho)).sqrt() * standard_normal(r);
    Ok((finite(sigma_x * x)?, finite(sigma_y * y)?))
}

/// Evaluates the joint zero-mean normal density at `(x, y)`.
///
/// Finite positive `sigma_x` and `sigma_y` are the component standard
/// deviations. Correlation `rho` must satisfy `-1 < rho < 1` so a joint density
/// exists. Returns zero for infinite observations and a domain error for NaN,
/// invalid deviations, or a degenerate correlation.
pub fn bivariate_normal_density(
    x: f64,
    y: f64,
    sigma_x: f64,
    sigma_y: f64,
    rho: f64,
) -> Result<f64, Error> {
    positive(sigma_x)?;
    positive(sigma_y)?;
    observation(x)?;
    observation(y)?;
    if !(-1.0..1.0).contains(&rho) || rho == -1.0 {
        return Err(Error::Domain);
    }
    if x.is_infinite() || y.is_infinite() {
        return Ok(0.0);
    }
    let u = x / sigma_x;
    let v = y / sigma_y;
    let s = (1.0 - rho) * (1.0 + rho);
    let quadratic = u * u + ((v - rho * u) / s.sqrt()).powi(2);
    density((-0.5 * quadratic - (2.0 * PI).ln() - sigma_x.ln() - sigma_y.ln() - 0.5 * s.ln()).exp())
}

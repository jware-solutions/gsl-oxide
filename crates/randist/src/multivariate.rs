// Adapted from GSL 2.8 dirichlet.c, mvgauss.c, wishart.c (GPL-3.0-or-later).
// Copyright (C) 2002 Gavin E. Crooks; (C) 2007 Brian Gough.
// Copyright (C) 2016 Timothee Flutre, Patrick Alken; (C) 2017 Timothee Flutre.
use crate::{
    density, finite, gamma, log_gamma, open, positive, standard_normal, xlog, Error, RngCore,
    LOG_SQRT_2_PI,
};

/// Draws a Dirichlet probability vector from random source `r`.
///
/// `concentrations` must be nonempty with finite positive entries. Returns
/// one nonnegative coordinate per entry, in the same order, summing to one
/// up to roundoff. Mean coordinate `i` is `concentrations[i] / sum(concentrations)`.
/// Logarithmic rescaling protects small concentrations; underflow can give zeros.
pub fn dirichlet<R: RngCore + ?Sized>(
    r: &mut R,
    concentrations: &[f64],
) -> Result<Vec<f64>, Error> {
    if concentrations.is_empty() {
        return Err(Error::EmptyInput);
    }
    for &a in concentrations {
        positive(a)?;
    }
    let mut logs = Vec::with_capacity(concentrations.len());
    for &a in concentrations {
        let log = if a < 1.0 {
            gamma(r, a + 1.0, 1.0)?.ln() + open(r).ln() / a
        } else {
            gamma(r, a, 1.0)?.ln()
        };
        logs.push(log);
    }
    let max = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return Err(Error::NumericalRange);
    }
    for x in &mut logs {
        *x = (*x - max).exp();
    }
    let sum: f64 = logs.iter().sum();
    for x in &mut logs {
        *x /= sum;
    }
    Ok(logs)
}

/// Returns the natural logarithm of the Dirichlet density at `coordinates`.
///
/// The two slices must have equal nonzero length; `concentrations` must be
/// finite and positive. Coordinates lie in `[0,1]` and must sum to one within
/// `32 * len * EPSILON`. Points outside that simplex return negative infinity;
/// NaN coordinates return a domain error. Boundary logarithms may be infinite,
/// or give a range error when competing singularities are indeterminate.
pub fn dirichlet_log_density(concentrations: &[f64], coordinates: &[f64]) -> Result<f64, Error> {
    if coordinates.iter().any(|x| x.is_nan()) {
        return Err(Error::Domain);
    }
    if concentrations.len() != coordinates.len() {
        return Err(Error::DimensionMismatch);
    }
    if concentrations.is_empty() {
        return Err(Error::EmptyInput);
    }
    for &a in concentrations {
        positive(a)?;
    }
    if coordinates
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0 || *x > 1.0)
    {
        return Ok(f64::NEG_INFINITY);
    }
    let sum: f64 = coordinates.iter().sum();
    if (sum - 1.0).abs() > 32.0 * coordinates.len() as f64 * f64::EPSILON {
        return Ok(f64::NEG_INFINITY);
    }
    let a: f64 = concentrations.iter().sum();
    let mut log = log_gamma(a)?;
    for (&a, &x) in concentrations.iter().zip(coordinates) {
        log += xlog(a - 1.0, x) - log_gamma(a)?;
    }
    density(log)
}

/// Evaluates a Dirichlet density at the simplex point `coordinates`.
///
/// `concentrations` supplies finite positive shapes in coordinate order;
/// both slices must have equal nonzero length. Returns a density with respect
/// to the first `len - 1` coordinates, zero outside the unit simplex, and
/// possibly infinity at a boundary. See [`dirichlet_log_density`] for the
/// normalization tolerance and logarithmic alternative.
pub fn dirichlet_density(concentrations: &[f64], coordinates: &[f64]) -> Result<f64, Error> {
    density(dirichlet_log_density(concentrations, coordinates)?.exp())
}

/// Checks that `matrix` stores an `n` by `n` flattened square.
///
/// Returns `Ok(())` when `n > 0` and its checked square matches the slice
/// length. Zero dimension returns `Error::EmptyInput`, incompatible or
/// overflowing dimensions `Error::DimensionMismatch`. Values are not inspected.
fn matrix(matrix: &[f64], n: usize) -> Result<(), Error> {
    if n == 0 {
        return Err(Error::EmptyInput);
    }
    if n.checked_mul(n) != Some(matrix.len()) {
        return Err(Error::DimensionMismatch);
    }
    Ok(())
}
/// Validates a row-major lower Cholesky factor `l` of dimension `n`.
///
/// Requires the square dimensions checked by [`matrix`], finite lower
/// entries, and finite strictly positive diagonal entries. Upper-triangular
/// entries are ignored. Returns `Ok(())` or a dimension/domain error.
fn factor(l: &[f64], n: usize) -> Result<(), Error> {
    matrix(l, n)?;
    for i in 0..n {
        positive(l[i * n + i])?;
        if l[i * n..i * n + i].iter().any(|x| !x.is_finite()) {
            return Err(Error::Domain);
        }
    }
    Ok(())
}
/// Solves `l * solution = b` in place by lower-triangular forward substitution.
///
/// The caller supplies a validated row-major factor `l` with dimension
/// `b.len()`. On success, `b` contains the solution and `Ok(())` is returned.
/// An unrepresentable component returns `Error::NumericalRange`, potentially
/// after earlier components of `b` have been replaced.
fn solve(l: &[f64], b: &mut [f64]) -> Result<(), Error> {
    let n = b.len();
    for i in 0..n {
        let sum: f64 = (0..i).map(|j| l[i * n + j] * b[j]).sum();
        b[i] = finite((b[i] - sum) / l[i * n + i])?;
    }
    Ok(())
}

/// Draws a multivariate normal vector with location `mean` using source `r`.
///
/// `mean` is a nonempty finite vector of length `n`; `cholesky` is a flattened
/// row-major `n` by `n` lower factor `L` with finite lower entries and positive
/// diagonal. Covariance is `L * L^T`; upper entries are ignored. Returns an
/// `n`-element sample in mean-coordinate order, with dimension/domain/range
/// errors for invalid inputs or unrepresentable calculations.
pub fn multivariate_normal<R: RngCore + ?Sized>(
    r: &mut R,
    mean: &[f64],
    cholesky: &[f64],
) -> Result<Vec<f64>, Error> {
    let n = mean.len();
    factor(cholesky, n)?;
    if mean.iter().any(|x| !x.is_finite()) {
        return Err(Error::Domain);
    }
    let z: Vec<_> = (0..n).map(|_| standard_normal(r)).collect();
    (0..n)
        .map(|i| finite(mean[i] + (0..=i).map(|j| cholesky[i * n + j] * z[j]).sum::<f64>()))
        .collect()
}
/// Returns the log multivariate normal density at finite vector `x`.
///
/// `x` and finite `mean` must have equal nonzero length `n`. `cholesky` is the
/// row-major `n` by `n` lower factor of the covariance, with finite lower entries
/// and positive diagonal; its upper entries are ignored. Returns the natural
/// logarithm of the joint density, preserving values whose density would
/// underflow. Invalid dimensions, nonfinite vectors, or arithmetic give an error.
pub fn multivariate_normal_log_density(
    x: &[f64],
    mean: &[f64],
    cholesky: &[f64],
) -> Result<f64, Error> {
    if x.len() != mean.len() {
        return Err(Error::DimensionMismatch);
    }
    let n = mean.len();
    factor(cholesky, n)?;
    if x.iter().chain(mean).any(|x| !x.is_finite()) {
        return Err(Error::Domain);
    }
    let mut z: Vec<_> = x.iter().zip(mean).map(|(x, m)| x - m).collect();
    solve(cholesky, &mut z)?;
    let det: f64 = (0..n).map(|i| cholesky[i * n + i].ln()).sum();
    density(-(n as f64) * LOG_SQRT_2_PI - det - 0.5 * z.iter().map(|v| v * v).sum::<f64>())
}
/// Evaluates a multivariate normal density at finite vector `x`.
///
/// `x` and finite `mean` must have equal nonzero length. `cholesky` is the
/// flattened row-major lower covariance factor described by
/// [`multivariate_normal_log_density`]. Returns the nonnegative joint density;
/// underflow gives zero. Invalid dimensions, parameters, or arithmetic give an error.
pub fn multivariate_normal_density(
    x: &[f64],
    mean: &[f64],
    cholesky: &[f64],
) -> Result<f64, Error> {
    density(multivariate_normal_log_density(x, mean, cholesky)?.exp())
}

/// Validates a row-major observation array and returns its number of rows.
///
/// `dimensions` is the nonzero column count and `minimum` the required row
/// count. `data` must be nonempty, finite, and contain complete rows. Returns
/// the row count, or an empty-input, dimension, or domain error on failure.
fn samples(data: &[f64], dimensions: usize, minimum: usize) -> Result<usize, Error> {
    if dimensions == 0 || data.is_empty() {
        return Err(Error::EmptyInput);
    }
    if data.len() % dimensions != 0 || data.len() / dimensions < minimum {
        return Err(Error::DimensionMismatch);
    }
    if data.iter().any(|x| !x.is_finite()) {
        return Err(Error::Domain);
    }
    Ok(data.len() / dimensions)
}
/// Computes the column means of row-major observation array `data`.
///
/// `dimensions` is the nonzero number of columns. `data` must be finite and
/// contain at least one complete row. Returns `dimensions` means in column
/// order, or an empty-input/dimension/domain error for invalid data and a
/// range error for unrepresentable arithmetic.
pub fn multivariate_mean(data: &[f64], dimensions: usize) -> Result<Vec<f64>, Error> {
    let rows = samples(data, dimensions, 1)?;
    let mut means = vec![0.0; dimensions];
    for j in 0..dimensions {
        let scale = data
            .chunks_exact(dimensions)
            .map(|row| row[j].abs())
            .fold(0.0, f64::max);
        if scale > 0.0 {
            means[j] = finite(
                data.chunks_exact(dimensions)
                    .map(|row| row[j] / scale / rows as f64)
                    .sum::<f64>()
                    * scale,
            )?;
        }
    }
    Ok(means)
}
/// Computes the unbiased sample covariance matrix of row-major `data`.
///
/// `dimensions` is the nonzero column count; finite `data` must contain at
/// least two complete observations. Returns a flattened row-major square:
/// diagonal entries are sample variances and off-diagonal entries covariances,
/// all divided by `rows - 1`. Invalid input or unrepresentable arithmetic
/// returns an error.
pub fn multivariate_sample_covariance(data: &[f64], dimensions: usize) -> Result<Vec<f64>, Error> {
    let rows = samples(data, dimensions, 2)?;
    let means = multivariate_mean(data, dimensions)?;
    let size = dimensions
        .checked_mul(dimensions)
        .ok_or(Error::DimensionMismatch)?;
    let mut covariance = vec![0.0; size];
    for i in 0..dimensions {
        for j in 0..=i {
            let sum: f64 = data
                .chunks_exact(dimensions)
                .map(|row| (row[i] - means[i]) * ((row[j] - means[j]) / (rows - 1) as f64))
                .sum();
            covariance[i * dimensions + j] = finite(sum)?;
            covariance[j * dimensions + i] = sum;
        }
    }
    Ok(covariance)
}

/// Draws a Wishart matrix with Bartlett decomposition using random source `r`.
///
/// `dimensions = n` must be positive, and finite `degrees` must exceed `n - 1`.
/// `scale_cholesky` is a row-major `n` by `n` lower factor `L`, with finite lower
/// entries and positive diagonal; upper entries are ignored. Returns a flattened
/// row-major symmetric matrix whose expectation is `degrees * L * L^T`, or an
/// error for invalid inputs or unrepresentable calculations.
pub fn wishart<R: RngCore + ?Sized>(
    r: &mut R,
    degrees: f64,
    scale_cholesky: &[f64],
    dimensions: usize,
) -> Result<Vec<f64>, Error> {
    let n = dimensions;
    factor(scale_cholesky, n)?;
    if !degrees.is_finite() || degrees <= (n - 1) as f64 {
        return Err(Error::Domain);
    }
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        a[i * n + i] = crate::chi_squared(r, degrees - i as f64)?.sqrt();
        for j in 0..i {
            a[i * n + j] = standard_normal(r);
        }
    }
    let mut b = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            b[i * n + j] = finite(
                (j..=i)
                    .map(|k| scale_cholesky[i * n + k] * a[k * n + j])
                    .sum(),
            )?;
        }
    }
    let mut x = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let v = finite((0..=j).map(|k| b[i * n + k] * b[j * n + k]).sum())?;
            x[i * n + j] = v;
            x[j * n + i] = v;
        }
    }
    Ok(x)
}

/// Returns a Wishart log density from the observation and scale factors.
///
/// `sample_cholesky` factors observation `X`, and `scale_cholesky` factors scale
/// `V`. Both are row-major `dimensions` by `dimensions` lower Cholesky factors,
/// with finite lower entries and positive diagonals; upper entries are ignored.
/// `dimensions` must be positive and finite `degrees > dimensions - 1`.
/// Returns the natural logarithm of the density of `X` under `Wishart(degrees,V)`;
/// invalid dimensions, factors, degrees, or arithmetic return an error.
pub fn wishart_log_density(
    sample_cholesky: &[f64],
    degrees: f64,
    scale_cholesky: &[f64],
    dimensions: usize,
) -> Result<f64, Error> {
    let n = dimensions;
    factor(sample_cholesky, n)?;
    factor(scale_cholesky, n)?;
    if !degrees.is_finite() || degrees <= (n - 1) as f64 {
        return Err(Error::Domain);
    }
    let mut trace = 0.0;
    for j in 0..n {
        let mut column: Vec<_> = (0..n)
            .map(|i| {
                if i >= j {
                    sample_cholesky[i * n + j]
                } else {
                    0.0
                }
            })
            .collect();
        solve(scale_cholesky, &mut column)?;
        trace += column.iter().map(|v| v * v).sum::<f64>();
    }
    let log_det_scale: f64 = (0..n).map(|i| 2.0 * scale_cholesky[i * n + i].ln()).sum();
    let log_det_sample: f64 = (0..n).map(|i| 2.0 * sample_cholesky[i * n + i].ln()).sum();
    let log_multigamma = (n * (n - 1)) as f64 * 0.25 * std::f64::consts::PI.ln()
        + (0..n)
            .map(|i| log_gamma(0.5 * (degrees - i as f64)))
            .collect::<Result<Vec<_>, _>>()?
            .iter()
            .sum::<f64>();
    density(
        -0.5 * degrees * n as f64 * 2.0_f64.ln() - 0.5 * degrees * log_det_scale - log_multigamma
            + 0.5 * (degrees - n as f64 - 1.0) * log_det_sample
            - 0.5 * trace,
    )
}
/// Evaluates a Wishart density using observation and scale Cholesky factors.
///
/// `sample_cholesky` and `scale_cholesky` are flattened row-major lower factors
/// of `dimensions` square matrices; both require finite lower entries and
/// positive diagonals. `dimensions > 0` and finite `degrees > dimensions - 1`.
/// Returns a nonnegative density, possibly zero from underflow; see
/// [`wishart_log_density`] for the logarithmic result and error behavior.
pub fn wishart_density(
    sample_cholesky: &[f64],
    degrees: f64,
    scale_cholesky: &[f64],
    dimensions: usize,
) -> Result<f64, Error> {
    density(wishart_log_density(sample_cholesky, degrees, scale_cholesky, dimensions)?.exp())
}

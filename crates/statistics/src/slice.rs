use crate::arithmetic::{
    average, deviation_scale, difference, finite, observation, rescale, Center, Sum,
};
use crate::Error;

/// Arithmetic mean `sum(data) / n` of at least one finite observation.
///
/// Returns a finite `f64` without computing a second moment. Empty data returns
/// [`Error::EmptyInput`]; NaN/infinity returns [`Error::NonFiniteData`]. See the
/// [crate numerical contracts](crate) and the [GSL mean definition](https://www.gnu.org/software/gsl/doc/html/statistics.html#mean-standard-deviation-and-variance).
pub fn mean(data: &[f64]) -> Result<f64, Error> {
    if data.is_empty() {
        return Err(Error::EmptyInput);
    }
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &value in data {
        observation(value)?;
        lo = lo.min(value);
        hi = hi.max(value);
    }
    if lo == hi {
        return Ok(lo);
    }
    average(data.iter().copied(), data.len(), data.len() as f64, lo, hi)
}

fn center(data: &[f64]) -> Result<Center, Error> {
    let mean = mean(data)?;
    let scale = deviation_scale(data, mean);
    let mut residual = Sum::default();
    for &value in data {
        residual.add(difference(value, mean, scale))?;
    }
    Ok(Center {
        mean,
        scale,
        residual: residual.total()? / data.len() as f64,
    })
}

// Recenter in normalized coordinates before squaring. This is a corrected
// centered two-pass calculation without subtracting nearly equal raw moments.
fn variance(data: &[f64], sample: bool) -> Result<f64, Error> {
    let c = center(data)?;
    if sample && data.len() < 2 {
        return Err(Error::InsufficientData);
    }
    let mut squares = Sum::default();
    for &value in data {
        let d = c.at(value);
        squares.add(d * d)?;
    }
    let denominator = (data.len() - usize::from(sample)) as f64;
    rescale(squares.total()? / denominator, c.scale, c.scale)
}

/// Population variance of finite `data`, centered about its computed mean.
///
/// Returns `M2 / n`; requires one observation (a singleton returns positive
/// zero). Errors: [`Error::EmptyInput`], [`Error::NonFiniteData`], or
/// [`Error::NumericalRange`] if the variance cannot be represented. See the
/// [crate contracts](crate) and [GSL definitions](https://www.gnu.org/software/gsl/doc/html/statistics.html#mean-standard-deviation-and-variance).
/// This does not accept a known fixed population mean.
pub fn population_variance(data: &[f64]) -> Result<f64, Error> {
    variance(data, false)
}

/// Sample variance `M2 / (n - 1)` of finite `data`, as in GSL's estimated variance.
///
/// Requires two observations; a nonempty singleton returns
/// [`Error::InsufficientData`]. Other errors and range limits are those of
/// [`population_variance`]. Returns a finite, nonnegative value.
/// [GSL reference](https://www.gnu.org/software/gsl/doc/html/statistics.html#mean-standard-deviation-and-variance).
pub fn sample_variance(data: &[f64]) -> Result<f64, Error> {
    variance(data, true)
}

/// Square root of [`population_variance(data)`](population_variance).
///
/// Requires one finite observation and returns a finite, nonnegative deviation.
/// The same errors apply, including [`Error::NumericalRange`] when the required
/// variance overflows even if its square root would fit.
/// [GSL reference](https://www.gnu.org/software/gsl/doc/html/statistics.html#mean-standard-deviation-and-variance).
pub fn population_standard_deviation(data: &[f64]) -> Result<f64, Error> {
    Ok(population_variance(data)?.sqrt())
}

/// Square root of [`sample_variance(data)`](sample_variance), using `n - 1`.
///
/// Requires two finite observations; returns a finite, nonnegative deviation.
/// All sample-variance errors apply, including overflow of the required variance.
/// [GSL reference](https://www.gnu.org/software/gsl/doc/html/statistics.html#mean-standard-deviation-and-variance).
pub fn sample_standard_deviation(data: &[f64]) -> Result<f64, Error> {
    Ok(sample_variance(data)?.sqrt())
}

fn paired(x: &[f64], y: &[f64]) -> Result<(Center, Center), Error> {
    if x.len() != y.len() {
        return Err(Error::LengthMismatch);
    }
    if x.len() < 2 {
        return Err(Error::InsufficientData);
    }
    Ok((center(x)?, center(y)?))
}

/// Sample covariance of equal-length finite `x` and `y`.
///
/// Returns `sum((x - mean_x)*(y - mean_y)) / (n - 1)`, requiring at least two
/// pairs. A constant input gives zero. Errors: [`Error::LengthMismatch`],
/// [`Error::InsufficientData`], [`Error::NonFiniteData`], or
/// [`Error::NumericalRange`] for an unrepresentable result. Centered products are
/// scaled, so individual public variances need not be representable.
/// [GSL covariance](https://www.gnu.org/software/gsl/doc/html/statistics.html#covariance).
pub fn sample_covariance(x: &[f64], y: &[f64]) -> Result<f64, Error> {
    let (cx, cy) = paired(x, y)?;
    let mut cross = Sum::default();
    for (&x, &y) in x.iter().zip(y) {
        cross.add(cx.at(x) * cy.at(y))?;
    }
    rescale(cross.total()? / (x.len() - 1) as f64, cx.scale, cy.scale)
}

/// Pearson correlation of equal-length finite `x` and `y`, in `[-1, 1]`.
///
/// Requires two pairs and nonzero centered variation in both inputs. Computes
/// `C / sqrt(M2_x) / sqrt(M2_y)` in scaled coordinates without forming a product
/// of moments or requiring representable variances. Errors are those of
/// [`sample_covariance`], plus [`Error::ZeroVariance`] for a constant input.
/// [GSL correlation](https://www.gnu.org/software/gsl/doc/html/statistics.html#correlation).
pub fn correlation(x: &[f64], y: &[f64]) -> Result<f64, Error> {
    let (cx, cy) = paired(x, y)?;
    let mut cross = Sum::default();
    let mut sx = Sum::default();
    let mut sy = Sum::default();
    for (&x, &y) in x.iter().zip(y) {
        let dx = cx.at(x);
        let dy = cy.at(y);
        cross.add(dx * dy)?;
        sx.add(dx * dx)?;
        sy.add(dy * dy)?;
    }
    let (sx, sy) = (sx.total()?, sy.total()?);
    if sx == 0.0 || sy == 0.0 {
        return Err(Error::ZeroVariance);
    }
    let r = finite((cross.total()? / sx.sqrt()) / sy.sqrt())?;
    // Compensated dot products, recentering, two roots and divisions introduce
    // a few ulps. Permit only 16 eps beyond Cauchy-Schwarz, not arbitrary clamping.
    if r.abs() > 1.0 + 16.0 * f64::EPSILON {
        Err(Error::NumericalRange)
    } else {
        Ok(r.clamp(-1.0, 1.0))
    }
}

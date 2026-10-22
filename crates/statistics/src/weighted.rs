use crate::arithmetic::{average, difference, observation, power_scale, rescale, Center, Sum};
use crate::Error;

struct Prepared {
    maximum: f64,
    total: f64,
    positive: usize,
    mean: f64,
}

fn prepare(data: &[f64], weights: &[f64]) -> Result<Prepared, Error> {
    if data.len() != weights.len() {
        return Err(Error::LengthMismatch);
    }
    let mut maximum: f64 = 0.0;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    let mut positive = 0;
    for (&value, &weight) in data.iter().zip(weights) {
        if !weight.is_finite() || weight < 0.0 {
            return Err(Error::InvalidWeight);
        }
        if weight == 0.0 {
            continue;
        }
        observation(value)?;
        maximum = maximum.max(weight);
        lo = lo.min(value);
        hi = hi.max(value);
        positive += 1; // Bounded by slice length, hence cannot overflow usize.
    }
    if positive == 0 {
        return Err(Error::ZeroTotalWeight);
    }
    let mut total = Sum::default();
    for &weight in weights {
        if weight == 0.0 {
            continue;
        }
        let q = weight / maximum;
        if q < f64::MIN_POSITIVE {
            return Err(Error::WeightRange);
        }
        total.add(q)?;
    }
    let total = total.total()?;
    let mean = if lo == hi {
        lo
    } else {
        average(
            data.iter()
                .zip(weights)
                .filter(|(_, w)| **w > 0.0)
                .map(|(&x, &w)| (w / maximum) * x),
            positive,
            total,
            lo,
            hi,
        )?
    };
    Ok(Prepared {
        maximum,
        total,
        positive,
        mean,
    })
}

/// Weighted mean `sum(w*x) / W` of equally long `data` and `weights`.
///
/// Weights must be finite and nonnegative; zero removes the observation before
/// validating its value. Requires one positive weight; returns its value for a
/// single participating observation. Errors: [`Error::LengthMismatch`],
/// [`Error::InvalidWeight`], [`Error::NonFiniteData`], [`Error::ZeroTotalWeight`],
/// [`Error::WeightRange`], or [`Error::NumericalRange`]. Does not compute moments.
/// See the [crate normalization/range contracts](crate) and
/// [GSL weighted samples](https://www.gnu.org/software/gsl/doc/html/statistics.html#weighted-samples).
pub fn weighted_mean(data: &[f64], weights: &[f64]) -> Result<f64, Error> {
    Ok(prepare(data, weights)?.mean)
}

fn variance(data: &[f64], weights: &[f64], reliability: bool) -> Result<f64, Error> {
    let p = prepare(data, weights)?;
    if reliability && p.positive < 2 {
        return Err(Error::InsufficientData);
    }
    let mut magnitude: f64 = 0.0;
    let mut delta: f64 = 0.0;
    for (&value, &weight) in data.iter().zip(weights) {
        if weight == 0.0 {
            continue;
        }
        magnitude = magnitude.max(value.abs());
        delta = delta.max((value - p.mean).abs());
    }
    let scale = power_scale(if delta.is_finite() { delta } else { magnitude });
    let mut residual = Sum::default();
    for (&value, &weight) in data.iter().zip(weights) {
        if weight == 0.0 {
            continue;
        }
        residual.add((weight / p.maximum) * difference(value, p.mean, scale))?;
    }
    let center = Center {
        mean: p.mean,
        scale,
        residual: residual.total()? / p.total,
    };
    let mut squares = Sum::default();
    let mut prefix = Sum::default();
    let mut pairs = Sum::default();
    for (&value, &weight) in data.iter().zip(weights) {
        if weight == 0.0 {
            continue;
        }
        let q = weight / p.maximum;
        let d = center.at(value);
        squares.add(q * (d * d))?;
        // W - W2/W = 2 * sum_{i<j}(q_i*q_j) / sum(q_i).
        // Positive cross terms avoid subtracting nearly equal totals when one
        // weight dominates. GSL's reliability correction, by algebra alone.
        pairs.add(q * prefix.total()?)?;
        prefix.add(q)?;
    }
    let denominator = if reliability {
        2.0 * (pairs.total()? / p.total)
    } else {
        p.total
    };
    if denominator <= 0.0 {
        return Err(Error::NumericalRange);
    }
    rescale(squares.total()? / denominator, scale, scale)
}

/// Weighted population-style variance `M2_w / W` about the computed weighted mean.
///
/// Requires equally long `data` and `weights`, with at least one positive
/// weight; one participating observation gives positive zero. Returns a finite,
/// nonnegative value. Zero weights remove even non-finite observations. Input,
/// normalization, and error contracts are those of [`weighted_mean`]; also returns
/// [`Error::NumericalRange`] if the variance is unrepresentable. This estimates
/// spread about the computed mean, not a supplied fixed population mean.
/// [GSL weighted definitions](https://www.gnu.org/software/gsl/doc/html/statistics.html#weighted-samples).
pub fn weighted_population_variance(data: &[f64], weights: &[f64]) -> Result<f64, Error> {
    variance(data, weights, false)
}

/// Reliability-weight variance `M2_w / (W - W2/W)`, where `W2 = sum(w*w)`.
///
/// This is GSL's weighted sample-variance correction, **not** frequency weighting
/// (`W - 1`), even for integer weights. Equal positive weights reproduce sample
/// variance; common positive weight rescaling preserves the estimator within
/// roundoff and the [supported normalized-weight range](crate).
///
/// Requires at least two positive weights, otherwise [`Error::InsufficientData`]
/// (or [`Error::ZeroTotalWeight`] if there are none). Returns a finite nonnegative
/// value. All input/overflow contracts of [`weighted_population_variance`] apply.
/// The denominator uses positive pair products to retain strongly unequal weights.
/// [GSL weighted samples](https://www.gnu.org/software/gsl/doc/html/statistics.html#weighted-samples).
pub fn weighted_reliability_variance(data: &[f64], weights: &[f64]) -> Result<f64, Error> {
    variance(data, weights, true)
}

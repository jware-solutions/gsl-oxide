use crate::Error;

/// Returns `value` when finite, or `NumericalRange` for overflow or NaN.
pub(crate) fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::NumericalRange)
    }
}

/// Validates that an observed `value` is finite, returning `NonFiniteData` otherwise.
pub(crate) fn observation(value: f64) -> Result<(), Error> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::NonFiniteData)
    }
}

// Neumaier summation, independently expressed from its compensated-sum identity.
// Also keeps the low part of a running mean when its increment is below one ulp.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Sum {
    pub(crate) value: f64,
    pub(crate) correction: f64,
}

impl Sum {
    /// Adds `value` to the compensated sum; reports unrepresentable intermediates without updating the sum.
    pub(crate) fn add(&mut self, value: f64) -> Result<(), Error> {
        let next = finite(self.value + value)?;
        let residual = if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        let correction = finite(self.correction + residual)?;
        self.value = next;
        self.correction = correction;
        Ok(())
    }

    /// Returns the compensated total, or a numerical-range error if it is not finite.
    pub(crate) fn total(self) -> Result<f64, Error> {
        if self.correction == 0.0 {
            finite(self.value)
        } else {
            finite(self.value + self.correction)
        }
    }
}

// First try the original units, preserving very small cancellation residuals.
// On overflow retry with a count-based power-of-two divisor: sum(abs(x))/divisor
// then fits whenever each x fits. Scaling by the *data maximum* instead would
// unnecessarily erase e.g. the 1e-300 term in [1e308, 1e-300, -1e308].
/// Sums `values` and divides by `denominator`, using `count` to choose a fallback scale and clamping roundoff to `[lo, hi]`.
pub(crate) fn average(
    values: impl Iterator<Item = f64> + Clone,
    count: usize,
    denominator: f64,
    lo: f64,
    hi: f64,
) -> Result<f64, Error> {
    let mut sum = Sum::default();
    let direct = values.clone().try_for_each(|x| sum.add(x));
    if direct.is_ok() {
        if let Ok(total) = sum.total() {
            return finite((total / denominator).clamp(lo, hi));
        }
    }
    let divisor = 2.0_f64.powi((usize::BITS - (count - 1).leading_zeros()) as i32);
    let mut sum = Sum::default();
    for value in values {
        sum.add(value / divisor)?;
    }
    // Divide by the (possibly noninteger) denominator before multiplying back
    // when necessary, otherwise apply the combined factor in one operation.
    let value = if divisor >= denominator {
        sum.total()? * (divisor / denominator)
    } else {
        (sum.total()? / denominator) * divisor
    };
    // A mean lies in the convex hull. This only repairs final rounding outside
    // that hull, including spurious overflow at the largest finite boundary.
    finite(value.clamp(lo, hi))
}

// Largest power of two <= magnitude, including subnormal powers.
/// Returns the largest power of two no greater than nonnegative `magnitude`, or one when it is zero.
pub(crate) fn power_scale(magnitude: f64) -> f64 {
    if magnitude == 0.0 {
        return 1.0;
    }
    let bits = magnitude.to_bits();
    let exponent = bits & (0x7ff_u64 << 52);
    if exponent != 0 {
        f64::from_bits(exponent)
    } else {
        f64::from_bits(1_u64 << (63 - bits.leading_zeros()))
    }
}

/// Returns the binary exponent of a positive power-of-two `scale`, including subnormal scales.
fn exponent(scale: f64) -> i32 {
    let bits = scale.to_bits();
    let field = (bits >> 52) as i32;
    if field == 0 {
        -1074 + (63 - bits.leading_zeros()) as i32
    } else {
        field - 1023
    }
}

// Combine binary exponents before scaling: neither scale product nor an
// intermediate square must fit. Monotone scaling allows ordinary underflow.
/// Multiplies `value` by power-of-two scales `a` and `b`, returning a numerical-range error on overflow.
pub(crate) fn rescale(mut value: f64, a: f64, b: f64) -> Result<f64, Error> {
    let mut exp = exponent(a) + exponent(b);
    while exp > 512 {
        value = finite(value * 2.0_f64.powi(512))?;
        exp -= 512;
    }
    while exp < -512 {
        value *= 2.0_f64.powi(-512);
        exp += 512;
    }
    finite(value * 2.0_f64.powi(exp))
}

/// Returns `(value - mean) / scale`, rearranging the subtraction when it would overflow.
pub(crate) fn difference(value: f64, mean: f64, scale: f64) -> f64 {
    let delta = value - mean;
    if delta.is_finite() {
        delta / scale
    } else {
        value / scale - mean / scale
    }
}

/// Finds a power-of-two scale for deviations of `data` from `mean`, avoiding overflowing differences.
pub(crate) fn deviation_scale(data: &[f64], mean: f64) -> f64 {
    let mut largest: f64 = 0.0;
    for &value in data {
        let delta = (value - mean).abs();
        if !delta.is_finite() {
            return power_scale(data.iter().fold(mean.abs(), |a, x| a.max(x.abs())));
        }
        largest = largest.max(delta);
    }
    power_scale(largest)
}

#[derive(Clone, Copy)]
pub(crate) struct Center {
    pub(crate) mean: f64,
    pub(crate) scale: f64,
    pub(crate) residual: f64,
}

impl Center {
    /// Returns the centered, scaled deviation of `value`, corrected for the stored mean residual.
    pub(crate) fn at(self, value: f64) -> f64 {
        difference(value, self.mean, self.scale) - self.residual
    }
}

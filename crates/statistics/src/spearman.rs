use crate::arithmetic::observation;
use crate::{pearson_correlation, Error};

/// Spearman rank correlation of equal-length finite `x` and `y`, in `[-1, 1]`.
///
/// Computes [`pearson_correlation`] of the ascending, one-based ranks of each
/// input, as in GSL 2.8's `gsl_stats_spearman`. Equal values receive the mean of
/// their occupied ranks; `-0.0` and `0.0` are tied. This measures monotonic
/// association and does not compute a significance level or p-value.
///
/// Returns [`Error::LengthMismatch`] for unequal lengths,
/// [`Error::InsufficientData`] for fewer than two pairs,
/// [`Error::NonFiniteData`] for NaN/infinite observations, or
/// [`Error::ZeroVariance`] for a constant input. These checks follow the crate's
/// checked Pearson contract; numerical range errors from [`pearson_correlation`]
/// also propagate.
///
/// Uses O(n log n) time and allocates O(n) auxiliary storage without mutating
/// either input. Unlike GSL's C API, no caller-supplied workspace is required.
/// [GSL definition](https://www.gnu.org/software/gsl/doc/html/statistics.html#correlation).
///
/// ```
/// use gsl_oxide_statistics::spearman_correlation;
/// let rho = spearman_correlation(&[1., 1., 2., 3.], &[1., 2., 2., 3.])?;
/// assert!((rho - 5. / 6.).abs() < 1e-14);
/// # Ok::<(), gsl_oxide_statistics::Error>(())
/// ```
pub fn spearman_correlation(x: &[f64], y: &[f64]) -> Result<f64, Error> {
    if x.len() != y.len() {
        return Err(Error::LengthMismatch);
    }
    if x.len() < 2 {
        return Err(Error::InsufficientData);
    }
    for &value in x.iter().chain(y) {
        observation(value)?;
    }
    pearson_correlation(&ranks(x), &ranks(y))
}

// Sort indices so ranks can be written back in the original pairing order.
// All values have been validated; numerical equality deliberately ties zeros
// even though total_cmp orders their signs separately.
/// Returns the ranks of `data` in input order, assigning tied values their average rank.
fn ranks(data: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..data.len()).collect();
    order.sort_unstable_by(|&a, &b| data[a].total_cmp(&data[b]));
    let mut ranks = vec![0.0; data.len()];
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && data[order[start]] == data[order[end]] {
            end += 1;
        }
        // The occupied one-based positions run from start + 1 through end.
        let rank = 0.5 * (start as f64 + 1.0 + end as f64);
        for &index in &order[start..end] {
            ranks[index] = rank;
        }
        start = end;
    }
    ranks
}

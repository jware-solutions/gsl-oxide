// Adapted from GSL 2.8 randist (GPL-3.0-or-later).
// Copyright (C) 1996-2000, 2007, 2009 James Theiler, Brian Gough.
use crate::{
    beta, continuous::positive_count, count, density, gamma, gamma_integer, log_beta, log_gamma,
    open, positive, probability, uniform, xlog, Error, RngCore,
};

/// Draws a Bernoulli sample with success probability `p` using source `r`.
///
/// `p` must lie in `[0,1]`. Returns one on success and zero otherwise;
/// endpoint probabilities produce deterministic samples. Invalid probabilities
/// return `Error::Domain`.
pub fn bernoulli<R: RngCore + ?Sized>(r: &mut R, p: f64) -> Result<u32, Error> {
    probability(p)?;
    Ok(u32::from(uniform(r) < p))
}
/// Returns the probability that a Bernoulli sample equals `k`.
///
/// For success probability `p` in `[0,1]`, returns `1-p` at zero, `p` at one,
/// and zero for other counts. Invalid `p` returns `Error::Domain`.
pub fn bernoulli_pdf(k: u32, p: f64) -> Result<f64, Error> {
    probability(p)?;
    Ok(match k {
        0 => 1.0 - p,
        1 => p,
        _ => 0.0,
    })
}

/// Draws a binomial success count from random source `r`.
///
/// `p` is the success probability in `[0,1]` and `trials` the number of
/// independent trials. Returns a count in `0..=trials`, with mean `trials*p`.
/// Uses Knuth's beta splitting, with geometric waiting times for small means.
pub fn binomial<R: RngCore + ?Sized>(r: &mut R, p: f64, trials: u32) -> Result<u32, Error> {
    probability(p)?;
    if p == 0.0 {
        return Ok(0);
    }
    if p == 1.0 {
        return Ok(trials);
    }
    let flipped = p > 0.5;
    let p = if flipped { 1.0 - p } else { p };
    let k = if f64::from(trials) * p < 30.0 {
        // Waiting times between successes avoid an underflowing cumulative PMF.
        let log_q = (-p).ln_1p();
        let mut remaining = trials;
        let mut k = 0;
        loop {
            let gap = (open(r).ln() / log_q).floor();
            if gap >= f64::from(remaining) {
                break;
            }
            remaining -= gap as u32 + 1;
            k += 1;
        }
        k
    } else {
        binomial_knuth(r, p, trials)?
    };
    Ok(if flipped { trials - k } else { k })
}

/// Draws a binomial success count using Knuth's beta-splitting algorithm.
///
/// `r` supplies randomness, `p` must lie in `[0,1]`, and `trials` is the
/// number of independent trials. Returns a count in `0..=trials`; zero trials
/// or zero success probability gives zero, while `p = 1` gives `trials`.
pub fn binomial_knuth<R: RngCore + ?Sized>(
    r: &mut R,
    mut p: f64,
    mut trials: u32,
) -> Result<u32, Error> {
    probability(p)?;
    if p == 0.0 {
        return Ok(0);
    }
    if p == 1.0 {
        return Ok(trials);
    }
    let mut k = 0;
    while trials > 10 {
        let a = 1 + trials / 2;
        let b = 1 + (trials - a);
        let x = beta(r, f64::from(a), f64::from(b))?;
        if x >= p {
            trials = a - 1;
            p /= x;
        } else {
            k += a;
            trials = b - 1;
            p = (p - x) / (1.0 - x);
        }
    }
    for _ in 0..trials {
        k += u32::from(uniform(r) < p);
    }
    Ok(k)
}

/// Returns the natural logarithm of the binomial coefficient `n choose k`.
///
/// An impossible choice `k > n` returns negative infinity. Callers supply
/// `n < u64::MAX`; small `k` uses a short logarithmic sum and larger `k`
/// uses log-beta. Special-function failures return a numerical range error.
pub(crate) fn log_choose(n: u64, k: u64) -> Result<f64, Error> {
    if k > n {
        return Ok(f64::NEG_INFINITY);
    }
    let k = k.min(n - k);
    if k < 128 {
        Ok((1..=k)
            .map(|j| ((n - k + j) as f64).ln() - (j as f64).ln())
            .sum())
    } else {
        Ok(-((n + 1) as f64).ln() - log_beta((k + 1) as f64, (n - k + 1) as f64)?)
    }
}

/// Returns the probability of `k` successes in `trials` independent trials.
///
/// The success probability `p` must lie in `[0,1]`. Returns
/// `choose(trials,k)*p^k*(1-p)^(trials-k)`, or zero for `k > trials`.
/// Endpoint probabilities are handled as point masses; invalid `p` is an error.
pub fn binomial_pdf(k: u32, p: f64, trials: u32) -> Result<f64, Error> {
    probability(p)?;
    if k > trials {
        return Ok(0.0);
    }
    if p == 0.0 {
        return Ok(f64::from(k == 0));
    }
    if p == 1.0 {
        return Ok(f64::from(k == trials));
    }
    density(
        (log_choose(u64::from(trials), u64::from(k))?
            + f64::from(k) * p.ln()
            + f64::from(trials - k) * (-p).ln_1p())
        .exp(),
    )
}

/// Draws a Poisson count with mean `mean` from random source `r`.
///
/// `mean` must be finite and nonnegative. Returns a nonnegative integer with
/// mean and variance `mean`; zero mean gives zero. Invalid means return a
/// domain error, and unrepresentable counts or intermediates a range error.
pub fn poisson<R: RngCore + ?Sized>(r: &mut R, mut mean: f64) -> Result<u32, Error> {
    if mean < 0.0 || !mean.is_finite() {
        return Err(Error::Domain);
    }
    let mut k: u32 = 0;
    while mean > 10.0 {
        let m = count((mean * 0.875).floor())?;
        let x = gamma_integer(r, m)?;
        if x >= mean {
            return k
                .checked_add(binomial(r, mean / x, m - 1)?)
                .ok_or(Error::NumericalRange);
        }
        k = k.checked_add(m).ok_or(Error::NumericalRange)?;
        mean -= x;
    }
    let threshold = (-mean).exp();
    let mut product = uniform(r);
    while product > threshold {
        k = k.checked_add(1).ok_or(Error::NumericalRange)?;
        product *= uniform(r);
    }
    Ok(k)
}
/// Fills `output` with independent Poisson counts drawn from source `r`.
///
/// `mean` must be finite and nonnegative and is checked even for an empty
/// output slice. Returns `Ok(())` after filling every element. A sampling
/// range error may leave an already-filled prefix in `output`.
pub fn fill_poisson<R: RngCore + ?Sized>(
    r: &mut R,
    output: &mut [u32],
    mean: f64,
) -> Result<(), Error> {
    if mean < 0.0 || !mean.is_finite() {
        return Err(Error::Domain);
    }
    for x in output {
        *x = poisson(r, mean)?;
    }
    Ok(())
}
/// Returns the probability that a Poisson count with mean `mean` equals `k`.
///
/// `mean` must be finite and nonnegative. Returns `mean^k*exp(-mean)/k!`,
/// with the zero-mean point mass handled explicitly. Invalid means return
/// a domain error; a very small probability may underflow to zero.
pub fn poisson_pdf(k: u32, mean: f64) -> Result<f64, Error> {
    if mean < 0.0 || !mean.is_finite() {
        return Err(Error::Domain);
    }
    if mean == 0.0 {
        return Ok(f64::from(k == 0));
    }
    if k == 0 {
        return Ok((-mean).exp());
    }
    let x = f64::from(k);
    let log = if k >= 16 {
        -crate::deviance(x, mean) - 0.5 * x.ln() - crate::LOG_SQRT_2_PI - crate::stirling(x)
    } else {
        x * mean.ln() - mean - log_gamma(x + 1.0)?
    };
    density(log.exp())
}
/// Draws the number of trials through the first success using source `r`.
///
/// The success probability satisfies `0 < p <= 1`. Returns a positive integer
/// with mean `1/p`, including the successful trial; `p = 1` gives one.
/// Invalid probabilities or counts beyond `u32::MAX` return an error.
pub fn geometric<R: RngCore + ?Sized>(r: &mut R, p: f64) -> Result<u32, Error> {
    probability(p)?;
    positive(p)?;
    if p == 1.0 {
        Ok(1)
    } else {
        positive_count(1.0 + open(r).ln() / (-p).ln_1p())
    }
}
/// Returns the probability of first success on trial `k`.
///
/// The success probability satisfies `0 < p <= 1`. Returns
/// `p*(1-p)^(k-1)` for `k >= 1` and zero at zero. `p = 1` is a point mass
/// at one; invalid probabilities return a domain error.
pub fn geometric_pdf(k: u32, p: f64) -> Result<f64, Error> {
    probability(p)?;
    positive(p)?;
    if k == 0 {
        Ok(0.0)
    } else if p == 1.0 {
        Ok(f64::from(k == 1))
    } else {
        Ok((p.ln() + f64::from(k - 1) * (-p).ln_1p()).exp())
    }
}
/// Draws a count of tagged items sampled without replacement using `r`.
///
/// `tagged` and `untagged` are the two population counts; `draws` is clamped
/// to their sum, following GSL. Returns the number of tagged items selected,
/// between zero and `min(tagged, draws)`; an empty population gives zero.
pub fn hypergeometric<R: RngCore + ?Sized>(
    r: &mut R,
    tagged: u32,
    untagged: u32,
    draws: u32,
) -> Result<u32, Error> {
    let n = u64::from(tagged) + u64::from(untagged);
    let t = u64::from(draws).min(n);
    let complement = t >= n / 2;
    let take = if complement { n - t } else { t };
    let mut a = tagged;
    let mut b = n;
    let mut k = 0;
    for _ in 0..take {
        if (b as f64) * uniform(r) < f64::from(a) {
            k += 1;
            a -= 1;
            if a == 0 {
                break;
            }
        }
        b -= 1;
    }
    Ok(if complement { tagged - k } else { k })
}
/// Returns the probability of selecting exactly `k` tagged items.
///
/// `tagged` and `untagged` specify the two population counts. As for the
/// sampler, `draws` is clamped to the population size. Returns the
/// hypergeometric probability, or zero for an impossible count; logarithmic
/// special-function failures return a range error.
pub fn hypergeometric_pdf(k: u32, tagged: u32, untagged: u32, draws: u32) -> Result<f64, Error> {
    let n = u64::from(tagged) + u64::from(untagged);
    let t = u64::from(draws).min(n);
    let k = u64::from(k);
    if k > u64::from(tagged) || k > t || t - k > u64::from(untagged) {
        return Ok(0.0);
    }
    density(
        (log_choose(u64::from(tagged), k)? + log_choose(u64::from(untagged), t - k)?
            - log_choose(n, t)?)
        .exp(),
    )
}
/// Draws a logarithmic-series count from random source `r`.
///
/// The distribution parameter satisfies `0 < p < 1`. Returns a positive
/// integer with probability `-p^k/(k*ln(1-p))`. Invalid probabilities or
/// samples beyond the `u32` count range return an error.
pub fn logarithmic_series<R: RngCore + ?Sized>(r: &mut R, p: f64) -> Result<u32, Error> {
    if !(p > 0.0 && p < 1.0) {
        return Err(Error::Domain);
    }
    let v = open(r);
    if v >= p {
        return Ok(1);
    }
    let q = -(((-p).ln_1p()) * open(r)).exp_m1();
    if v <= q * q {
        positive_count(1.0 + v.ln() / q.ln())
    } else if v <= q {
        Ok(2)
    } else {
        Ok(1)
    }
}
/// Returns the logarithmic-series probability at integer count `k`.
///
/// For `0 < p < 1`, returns `-p^k/(k*ln(1-p))` when `k >= 1`, and zero
/// at zero. Invalid `p` returns a domain error; small probabilities may
/// underflow to zero.
pub fn logarithmic_series_pdf(k: u32, p: f64) -> Result<f64, Error> {
    if !(p > 0.0 && p < 1.0) {
        return Err(Error::Domain);
    }
    if k == 0 {
        Ok(0.0)
    } else {
        density((f64::from(k) * p.ln() - f64::from(k).ln() - (-(-p).ln_1p()).ln()).exp())
    }
}
/// Draws failures before a given number of successes using source `r`.
///
/// `p` is the success probability with `0 < p <= 1`; finite positive
/// `successes` is the shape and may be noninteger. Returns a nonnegative
/// count with mean `successes*(1-p)/p`. `p = 1` gives zero; invalid parameters
/// or unrepresentable counts return an error.
pub fn negative_binomial<R: RngCore + ?Sized>(
    r: &mut R,
    p: f64,
    successes: f64,
) -> Result<u32, Error> {
    probability(p)?;
    positive(p)?;
    positive(successes)?;
    if p == 1.0 {
        return Ok(0);
    }
    let mean = gamma(r, successes, (1.0 - p) / p)?;
    poisson(r, mean)
}
/// Returns the negative-binomial probability of exactly `k` failures.
///
/// `p` is the success probability with `0 < p <= 1`; finite positive
/// `successes` is the possibly noninteger success-count shape. Returns
/// `Gamma(successes+k)/(Gamma(successes)*k!)*p^successes*(1-p)^k`.
/// For `p = 1`, all probability is at zero.
pub fn negative_binomial_pdf(k: u32, p: f64, successes: f64) -> Result<f64, Error> {
    probability(p)?;
    positive(p)?;
    positive(successes)?;
    if p == 1.0 {
        return Ok(f64::from(k == 0));
    }
    let kf = f64::from(k);
    density(
        (-log_beta(successes, kf + 1.0)? - (successes + kf).ln()
            + successes * p.ln()
            + kf * (-p).ln_1p())
        .exp(),
    )
}
/// Draws failures before an integer number of successes using source `r`.
///
/// `p` must satisfy `0 < p <= 1` and `successes` must be a positive integer.
/// Returns a nonnegative negative-binomial count with mean
/// `successes*(1-p)/p`; invalid parameters or count overflow returns an error.
pub fn pascal<R: RngCore + ?Sized>(r: &mut R, p: f64, successes: u32) -> Result<u32, Error> {
    negative_binomial(r, p, f64::from(successes))
}
/// Returns the Pascal probability of `k` failures before `successes` successes.
///
/// `p` is the success probability with `0 < p <= 1` and `successes` a
/// positive integer. Returns the corresponding negative-binomial probability;
/// `p = 1` puts all probability at zero.
pub fn pascal_pdf(k: u32, p: f64, successes: u32) -> Result<f64, Error> {
    negative_binomial_pdf(k, p, f64::from(successes))
}

/// Walker alias table for repeated sampling from nonnegative category weights.
#[derive(Clone, Debug)]
pub struct DiscreteDistribution {
    /// Normalized category probabilities, in the original weight order.
    probabilities: Vec<f64>,
    /// Per-column thresholds for selecting the column instead of its alias.
    cutoffs: Vec<f64>,
    /// Fallback category for each column when its threshold test fails.
    aliases: Vec<usize>,
}
impl DiscreteDistribution {
    /// Builds an alias table from a nonempty slice of category `weights`.
    ///
    /// Weights must be finite and nonnegative, with at least one positive entry;
    /// they need not sum to one. Returns normalized probabilities and an O(n)
    /// preprocessing table for reuse. Empty input gives `Error::EmptyInput`, and
    /// invalid or all-zero weights give `Error::Domain`.
    pub fn new(weights: &[f64]) -> Result<Self, Error> {
        let probabilities = normalized_weights(weights)?;
        let n = weights.len();
        let mut cutoffs: Vec<_> = probabilities.iter().map(|p| p * n as f64).collect();
        let mut aliases: Vec<_> = (0..n).collect();
        let (mut small, mut large): (Vec<_>, Vec<_>) = (0..n).partition(|&i| cutoffs[i] < 1.0);
        while !small.is_empty() && !large.is_empty() {
            let s = small.pop().unwrap();
            let l = large.pop().unwrap();
            aliases[s] = l;
            cutoffs[l] = (cutoffs[l] + cutoffs[s]) - 1.0;
            if cutoffs[l] < 1.0 {
                small.push(l);
            } else {
                large.push(l);
            }
        }
        for i in small.into_iter().chain(large) {
            cutoffs[i] = 1.0;
        }
        Ok(Self {
            probabilities,
            cutoffs,
            aliases,
        })
    }
    /// Draws a category index from this alias table using random source `r`.
    ///
    /// Returns an index in `0..weights.len()` with the normalized probability
    /// established by [`Self::new`]. Sampling takes O(1) expected time and does
    /// not rebuild or modify the table.
    pub fn sample<R: RngCore + ?Sized>(&self, r: &mut R) -> usize {
        let i = crate::sampling::index(r, self.probabilities.len());
        if uniform(r) < self.cutoffs[i] {
            i
        } else {
            self.aliases[i]
        }
    }
    /// Returns the probability of the given category index.
    ///
    /// `category` indexes the weights passed to [`Self::new`]. Returns its
    /// normalized probability in `[0,1]`, or zero for an out-of-range index.
    /// This query does not draw from a random source.
    pub fn pdf(&self, category: usize) -> f64 {
        self.probabilities.get(category).copied().unwrap_or(0.0)
    }
}

/// Converts finite nonnegative `weights` into category probabilities.
///
/// Returns a nonempty vector with the same length and a sum approximately
/// one, scaling by the maximum weight to avoid sum overflow. Empty input
/// returns `Error::EmptyInput`; invalid or all-zero weights give a domain error.
/// Extremely small normalized weights may underflow to zero.
pub(crate) fn normalized_weights(weights: &[f64]) -> Result<Vec<f64>, Error> {
    if weights.is_empty() {
        return Err(Error::EmptyInput);
    }
    if weights.iter().any(|x| !x.is_finite() || *x < 0.0) {
        return Err(Error::Domain);
    }
    let max = weights.iter().copied().fold(0.0, f64::max);
    if max == 0.0 {
        return Err(Error::Domain);
    }
    let total: f64 = weights.iter().map(|x| x / max).sum();
    Ok(weights.iter().map(|x| (x / max) / total).collect())
}

/// Draws category counts for `trials` independent weighted choices using `r`.
///
/// `weights` must be nonempty, finite, nonnegative, and have a positive sum;
/// they need not be normalized. Returns one count per category, in weight
/// order, whose sum is exactly `trials`. Zero-weight categories receive zero.
pub fn multinomial<R: RngCore + ?Sized>(
    r: &mut R,
    trials: u32,
    weights: &[f64],
) -> Result<Vec<u32>, Error> {
    let p = normalized_weights(weights)?;
    let mut output = vec![0; p.len()];
    let mut remaining = trials;
    // Reverse sums avoid subtracting almost equal floating-point totals.
    let mut tails = p.clone();
    for i in (0..p.len() - 1).rev() {
        tails[i] += tails[i + 1];
    }
    for i in 0..p.len() {
        let n = if i + 1 == p.len() {
            remaining
        } else if p[i] == 0.0 {
            0
        } else {
            binomial(r, (p[i] / tails[i]).min(1.0), remaining)?
        };
        output[i] = n;
        remaining -= n;
    }
    Ok(output)
}
/// Returns the natural logarithm of a multinomial count probability.
///
/// `weights` gives nonnegative category weights and `counts` the observed
/// counts in the same order. Slices must have equal nonzero length and the
/// weights must be finite with at least one positive entry. The trial count
/// is the sum of `counts`. Impossible observations return negative infinity;
/// invalid dimensions, weights, or numerical range return an error.
pub fn multinomial_log_pdf(weights: &[f64], counts: &[u32]) -> Result<f64, Error> {
    if weights.len() != counts.len() {
        return Err(Error::DimensionMismatch);
    }
    let p = normalized_weights(weights)?;
    let total = counts
        .iter()
        .try_fold(0_u64, |s, &n| s.checked_add(u64::from(n)))
        .ok_or(Error::NumericalRange)?;
    let mut log = log_gamma(total as f64 + 1.0)?;
    for (&n, &p) in counts.iter().zip(&p) {
        log += xlog(f64::from(n), p) - log_gamma(f64::from(n) + 1.0)?;
    }
    density(log)
}
/// Returns the multinomial probability of the observed category `counts`.
///
/// `weights` and `counts` must have equal nonzero length; weights must be
/// finite, nonnegative, and not all zero, but need not sum to one. The trial
/// count is inferred from `counts`. Returns zero for impossible observations
/// or probability underflow; see [`multinomial_log_pdf`] for the logarithm.
pub fn multinomial_pdf(weights: &[f64], counts: &[u32]) -> Result<f64, Error> {
    density(multinomial_log_pdf(weights, counts)?.exp())
}

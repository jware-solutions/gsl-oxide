use crate::{beta_cdf, beta_survival, positive, probability, shaped::gamma_unit_tails, Error};
use gsl_oxide_special::beta::ln_beta;

/// Returns binomial lower and upper tails at count `k` for success probability `p` and `trials` trials.
fn binomial_tails(k: u32, p: f64, trials: u32) -> Result<(f64, f64), Error> {
    probability(p)?;
    if k >= trials {
        return Ok((1.0, 0.0));
    }
    let a = f64::from(k) + 1.0;
    let b = f64::from(trials - k);
    Ok((beta_survival(p, a, b)?, beta_cdf(p, a, b)?))
}
/// Binomial `Pr(X <= k)` for success probability `p` in `[0,1]` and `trials` trials.
/// The incomplete-beta shape limit in the crate documentation applies.
pub fn binomial_cdf(k: u32, p: f64, trials: u32) -> Result<f64, Error> {
    Ok(binomial_tails(k, p, trials)?.0)
}
/// Binomial `Pr(X > k)`, retaining small upper tails; parameters follow [`binomial_cdf`].
pub fn binomial_survival(k: u32, p: f64, trials: u32) -> Result<f64, Error> {
    Ok(binomial_tails(k, p, trials)?.1)
}

/// Returns Poisson lower and upper tails at count `k` for the nonnegative `mean`.
fn poisson_tails(k: u32, mean: f64) -> Result<(f64, f64), Error> {
    if !mean.is_finite() || mean < 0.0 {
        return Err(Error::Domain);
    }
    if mean == 0.0 {
        return Ok((1.0, 0.0));
    }
    let (q, p) = gamma_unit_tails(f64::from(k) + 1.0, mean)?;
    Ok((p, q))
}
/// Poisson `Pr(X <= k)` with finite nonnegative `mean`; zero is a point mass at zero.
/// The gamma shape limit applies to `k+1` when `mean>0`.
pub fn poisson_cdf(k: u32, mean: f64) -> Result<f64, Error> {
    Ok(poisson_tails(k, mean)?.0)
}
/// Poisson `Pr(X > k)`, retaining small upper tails; parameters follow [`poisson_cdf`].
pub fn poisson_survival(k: u32, mean: f64) -> Result<f64, Error> {
    Ok(poisson_tails(k, mean)?.1)
}

/// Returns geometric lower and upper tails at trial count `k` for success probability `p`.
fn geometric_tails(k: u32, p: f64) -> Result<(f64, f64), Error> {
    probability(p)?;
    if k == 0 {
        return Ok((0.0, 1.0));
    }
    let log = f64::from(k) * (-p).ln_1p();
    Ok((-log.exp_m1(), log.exp()))
}
/// Geometric `Pr(X <= k)` counting trials through the first success (support starts at one).
/// `p` is in `[0,1]`; at zero all finite-count CDFs are zero, as in GSL.
pub fn geometric_cdf(k: u32, p: f64) -> Result<f64, Error> {
    Ok(geometric_tails(k, p)?.0)
}
/// Geometric `Pr(X > k)`; parameters and counting convention follow [`geometric_cdf`].
pub fn geometric_survival(k: u32, p: f64) -> Result<f64, Error> {
    Ok(geometric_tails(k, p)?.1)
}

/// Returns negative-binomial tails at failure count `k` for probability `p` and required `successes`.
fn negative_binomial_tails(k: u32, p: f64, successes: f64) -> Result<(f64, f64), Error> {
    probability(p)?;
    positive(successes)?;
    Ok((
        beta_cdf(p, successes, f64::from(k) + 1.0)?,
        beta_survival(p, successes, f64::from(k) + 1.0)?,
    ))
}
/// Negative-binomial `Pr(X <= k)` counting failures before positive real `successes` successes.
/// `p` is in `[0,1]`; incomplete-beta shape limits apply to `successes` and `k+1`.
pub fn negative_binomial_cdf(k: u32, p: f64, successes: f64) -> Result<f64, Error> {
    Ok(negative_binomial_tails(k, p, successes)?.0)
}
/// Negative-binomial `Pr(X > k)`; parameters follow [`negative_binomial_cdf`].
pub fn negative_binomial_survival(k: u32, p: f64, successes: f64) -> Result<f64, Error> {
    Ok(negative_binomial_tails(k, p, successes)?.1)
}
/// Pascal CDF: negative binomial with a positive integer number of `successes`.
pub fn pascal_cdf(k: u32, p: f64, successes: u32) -> Result<f64, Error> {
    negative_binomial_cdf(k, p, f64::from(successes))
}
/// Pascal survival: negative binomial with a positive integer number of `successes`.
pub fn pascal_survival(k: u32, p: f64, successes: u32) -> Result<f64, Error> {
    negative_binomial_survival(k, p, f64::from(successes))
}

/// Returns the natural logarithm of the binomial coefficient `n choose k` for `k <= n`.
fn log_choose(n: u64, k: u64) -> Result<f64, Error> {
    let k = k.min(n - k);
    if k < 128 {
        Ok((1..=k)
            .map(|i| ((n - k + i) as f64).ln() - (i as f64).ln())
            .sum())
    } else {
        Ok(-((n + 1) as f64).ln() - ln_beta((k + 1) as f64, (n - k + 1) as f64)?)
    }
}

/// Returns hypergeometric tails at `k` for `tagged` and `untagged` population counts and `draws` draws, or an error for invalid inputs or nonconvergence.
fn hypergeometric_tails(
    k: u32,
    tagged: u32,
    untagged: u32,
    draws: u32,
) -> Result<(f64, f64), Error> {
    let n1 = u64::from(tagged);
    let n2 = u64::from(untagged);
    let t = u64::from(draws);
    if t > n1 + n2 {
        return Err(Error::Domain);
    }
    let minimum = t.saturating_sub(n2);
    let maximum = t.min(n1);
    let k = u64::from(k);
    if k < minimum {
        return Ok((0.0, 1.0));
    }
    if k >= maximum {
        return Ok((1.0, 0.0));
    }
    let upper = k as f64 >= (t as f64) * (n1 as f64) / ((n1 + n2) as f64);
    let mut i = if upper { k + 1 } else { k };
    let mut term = (log_choose(n1, i)? + log_choose(n2, t - i)? - log_choose(n1 + n2, t)?).exp();
    let mut sum = term;
    for _ in 0..1_000_000 {
        if (upper && i == maximum) || (!upper && i == minimum) || term == 0.0 {
            let sum = sum.clamp(0.0, 1.0);
            return Ok(if upper {
                (1.0 - sum, sum)
            } else {
                (sum, 1.0 - sum)
            });
        }
        let factor = if upper {
            ((n1 - i) as f64) / ((i + 1) as f64) * ((t - i) as f64) / ((n2 + i + 1 - t) as f64)
        } else {
            (i as f64) / ((n1 - i + 1) as f64) * ((n2 + i - t) as f64) / ((t - i + 1) as f64)
        };
        term *= factor;
        sum += term;
        if term <= f64::EPSILON * sum {
            return Ok(if upper {
                (1.0 - sum, sum)
            } else {
                (sum, 1.0 - sum)
            });
        }
        if upper {
            i += 1;
        } else {
            i -= 1;
        }
    }
    Err(Error::NonConvergence)
}
/// Hypergeometric `Pr(X <= k)` for `tagged,untagged` population counts and `draws` draws.
/// Requires `draws <= tagged+untagged` (unlike randist's clamped sampler). Population
/// addition is widened to avoid `u32` overflow; impossible counts return support limits.
pub fn hypergeometric_cdf(k: u32, tagged: u32, untagged: u32, draws: u32) -> Result<f64, Error> {
    Ok(hypergeometric_tails(k, tagged, untagged, draws)?.0)
}
/// Hypergeometric `Pr(X > k)`, summing the smaller tail directly; parameters follow [`hypergeometric_cdf`].
pub fn hypergeometric_survival(
    k: u32,
    tagged: u32,
    untagged: u32,
    draws: u32,
) -> Result<f64, Error> {
    Ok(hypergeometric_tails(k, tagged, untagged, draws)?.1)
}

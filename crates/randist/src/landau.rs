// Adapted from GSL 2.8 landau.c (GPL-3.0-or-later).
// Copyright (C) 2001, 2004 David Morrison; CERN DENLAN/RANLAN by K. S. Kolbig.
use crate::{density, observation, open, Error, RngCore};
/// DENLAN rational coefficients and RANLAN inverse-distribution lookup values.
mod tables;
use tables::*;

/// Evaluates coefficients `c` at `x` with Horner's rule.
///
/// Coefficients are ordered from the constant term upward; an empty slice
/// returns zero. Callers supply the DENLAN approximation's coefficients.
fn polynomial(c: &[f64], x: f64) -> f64 {
    c.iter().rev().fold(0.0, |s, &c| s * x + c)
}
/// Returns the ratio of the polynomials with coefficients `p` and `q` at `x`.
///
/// Both slices use ascending powers. Callers choose an approximation region
/// where the denominator is nonzero; floating-point division is unchecked.
fn rational(p: &[f64], q: &[f64], x: f64) -> f64 {
    polynomial(p, x) / polynomial(q, x)
}

/// Evaluates the standard Landau density at observation `x`.
///
/// Returns a nonnegative density using GSL/CERN DENLAN's piecewise
/// approximation, whose coefficients have about nine significant digits.
/// Infinite `x` returns zero; NaN returns `Error::Domain`.
pub fn landau_density(x: f64) -> Result<f64, Error> {
    observation(x)?;
    if x.is_infinite() {
        return Ok(0.0);
    }
    let y = if x < -5.5 {
        let u = (x + 1.0).exp();
        if u == 0.0 {
            return Ok(0.0);
        }
        0.3989422803 * (-1.0 / u).exp() / u.sqrt() * (1.0 + u * polynomial(&A1, u))
    } else if x < -1.0 {
        let u = (-x - 1.0).exp();
        (-u).exp() * u.sqrt() * rational(&P1, &Q1, x)
    } else if x < 1.0 {
        rational(&P2, &Q2, x)
    } else if x < 5.0 {
        rational(&P3, &Q3, x)
    } else if x < 12.0 {
        let u = 1.0 / x;
        u * u * rational(&P4, &Q4, u)
    } else if x < 50.0 {
        let u = 1.0 / x;
        u * u * rational(&P5, &Q5, u)
    } else if x < 300.0 {
        let u = 1.0 / x;
        u * u * rational(&P6, &Q6, u)
    } else {
        let u = (1.0 / x) / (1.0 - x.ln() / (x + 1.0));
        u * u * (1.0 + u * polynomial(&A2, u))
    };
    density(y)
}

/// Draws a standard Landau sample using random source `r`.
///
/// Returns one real value in the GSL/CERN standard Landau parameterization,
/// using RANLAN's inverse table and tail fits. The distribution has a long
/// right tail and no finite mean or variance.
pub fn landau<R: RngCore + ?Sized>(r: &mut R) -> f64 {
    let x = open(r);
    let v = 1000.0 * x;
    let i = v as usize;
    let u = v - i as f64;
    if (70..=800).contains(&i) {
        F[i] + u * (F[i + 1] - F[i])
    } else if (7..=980).contains(&i) {
        F[i] + u * (F[i + 1] - F[i] - 0.25 * (1.0 - u) * (F[i + 2] - F[i + 1] - F[i] + F[i - 1]))
    } else if i < 7 {
        let v = x.ln();
        let u = 1.0 / v;
        (0.99858950 + (34.5213058 + 17.0854528 * u) * u) / (1.0 + (34.1760202 + 4.01244582 * u) * u)
            * (-(-0.91893853 - v).ln() - 1.0)
    } else {
        let u = 1.0 - x;
        let v = u * u;
        if x <= 0.999 {
            (1.00060006 + 263.991156 * u + 4373.20068 * v)
                / ((1.0 + 257.368075 * u + 3414.48018 * v) * u)
        } else {
            (1.00001538 + 6075.14119 * u + 734266.409 * v)
                / ((1.0 + 6065.11919 * u + 694021.044 * v) * u)
        }
    }
}

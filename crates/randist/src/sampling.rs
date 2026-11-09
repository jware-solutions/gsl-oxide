// Adapted from GSL 2.8 shuffle.c and sphere.c (GPL-3.0-or-later).
// Copyright (C) 1996-2000, 2007 James Theiler, Brian Gough.
use crate::{finite, open, standard_normal, Error, RngCore};
use std::f64::consts::TAU;

/// Draws an unbiased integer index in `0..n` from random source `r`.
///
/// The caller must supply `n > 0`. Rejection of the incomplete remainder
/// of the `u64` range avoids modulo bias; zero `n` would panic.
pub(crate) fn index<R: RngCore + ?Sized>(r: &mut R, n: usize) -> usize {
    let n = n as u64;
    let threshold = n.wrapping_neg() % n;
    loop {
        let x = r.next_u64();
        if x >= threshold {
            return (x % n) as usize;
        }
    }
}

/// Randomly permutes the elements of `values` in place using source `r`.
///
/// Fisher-Yates gives each permutation equal probability. The function
/// returns `()`; empty and one-element slices require no random draws.
pub fn shuffle<R: RngCore + ?Sized, T>(r: &mut R, values: &mut [T]) {
    for i in (1..values.len()).rev() {
        let j = index(r, i + 1);
        values.swap(i, j);
    }
}

/// Chooses `count` items from `source` without replacement using `r`.
///
/// Returns cloned items in their original source order, so no source position
/// is selected twice. `count = 0` returns an empty vector; counts larger than
/// `source.len()` return `Error::SampleTooLarge`.
pub fn choose<R: RngCore + ?Sized, T: Clone>(
    r: &mut R,
    source: &[T],
    count: usize,
) -> Result<Vec<T>, Error> {
    if count > source.len() {
        return Err(Error::SampleTooLarge);
    }
    let mut output = Vec::with_capacity(count);
    for (i, x) in source.iter().enumerate() {
        if output.len() == count {
            break;
        }
        if index(r, source.len() - i) < count - output.len() {
            output.push(x.clone());
        }
    }
    Ok(output)
}

/// Draws `count` items independently from `source` using random source `r`.
///
/// Returns a vector of cloned items in draw order; repeated source positions
/// are allowed. An empty source is valid only for `count = 0`, which returns
/// an empty vector; otherwise it returns `Error::EmptyInput`.
pub fn sample_with_replacement<R: RngCore + ?Sized, T: Clone>(
    r: &mut R,
    source: &[T],
    count: usize,
) -> Result<Vec<T>, Error> {
    if count > 0 && source.is_empty() {
        return Err(Error::EmptyInput);
    }
    Ok((0..count)
        .map(|_| source[index(r, source.len())].clone())
        .collect())
}

/// Draws a uniform direction on the unit circle using random source `r`.
///
/// Returns Cartesian coordinates `[x, y]` with `x^2 + y^2` approximately one.
/// Uses disk rejection and a rational map to avoid trigonometric functions.
pub fn direction_2d<R: RngCore + ?Sized>(r: &mut R) -> [f64; 2] {
    loop {
        let x = 2.0 * open(r) - 1.0;
        let y = 2.0 * open(r) - 1.0;
        let s = x * x + y * y;
        if s > 0.0 && s <= 1.0 {
            return [(x * x - y * y) / s, 2.0 * x * y / s];
        }
    }
}
/// Draws a uniform direction on the unit circle using random source `r`.
///
/// Returns `[cos(angle), sin(angle)]` for a uniform angle in `(0, 2*pi)`.
/// The result is a unit vector up to floating-point roundoff.
pub fn direction_2d_trigonometric<R: RngCore + ?Sized>(r: &mut R) -> [f64; 2] {
    let angle = TAU * open(r);
    [angle.cos(), angle.sin()]
}
/// Draws a uniform direction on the unit sphere using random source `r`.
///
/// Returns Cartesian coordinates `[x, y, z]` whose squared norm is
/// approximately one, using Marsaglia's disk-rejection method.
pub fn direction_3d<R: RngCore + ?Sized>(r: &mut R) -> [f64; 3] {
    loop {
        let x = 2.0 * open(r) - 1.0;
        let y = 2.0 * open(r) - 1.0;
        let s = x * x + y * y;
        if s <= 1.0 {
            let t = 2.0 * (1.0 - s).sqrt();
            return [x * t, y * t, 1.0 - 2.0 * s];
        }
    }
}
/// Draws a uniform unit direction in `dimensions` dimensions using `r`.
///
/// Returns a vector of `dimensions` normalized independent standard normals.
/// Its Euclidean norm is approximately one. Zero dimensions returns
/// `Error::EmptyInput`; an unrepresentable normalization returns a range error.
pub fn direction_nd<R: RngCore + ?Sized>(r: &mut R, dimensions: usize) -> Result<Vec<f64>, Error> {
    if dimensions == 0 {
        return Err(Error::EmptyInput);
    }
    loop {
        let mut x: Vec<_> = (0..dimensions).map(|_| standard_normal(r)).collect();
        let norm = x.iter().fold(0.0_f64, |s, &v| s.hypot(v));
        if norm == 0.0 {
            continue;
        }
        finite(norm)?;
        for v in &mut x {
            *v /= norm;
        }
        return Ok(x);
    }
}

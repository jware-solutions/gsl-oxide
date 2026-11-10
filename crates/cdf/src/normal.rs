// Cody's rational normal tails and Wichura's AS 241 inverse, from GSL 2.8.
#![allow(clippy::excessive_precision)]
use crate::{finite_quantile, observation, positive, probability, Error};

/// Returns a normal tail from positive `x` and the supplied `rational` factor, splitting its exponential for accuracy.
fn scaled_tail(x: f64, rational: f64) -> f64 {
    let rounded = (x * 16.0).floor() / 16.0;
    let remainder = 0.5 * (x - rounded) * (x + rounded);
    (-0.5 * rounded * rounded).exp() * (-remainder).exp() * rational
}
/// Returns the normal CDF offset from one half for a small argument `x`.
fn small(x: f64) -> f64 {
    const A: [f64; 5] = [
        2.2352520354606839287,
        161.02823106855587881,
        1067.6894854603709582,
        18154.981253343561249,
        0.065682337918207449113,
    ];
    const B: [f64; 4] = [
        47.20258190468824187,
        976.09855173777669322,
        10260.932208618978205,
        45507.789335026729956,
    ];
    let square = x * x;
    let mut num = A[4] * square;
    let mut den = square;
    for i in 0..3 {
        num = (num + A[i]) * square;
        den = (den + B[i]) * square;
    }
    x * (num + A[3]) / (den + B[3])
}
/// Returns the standard normal survival probability at nonnegative `x`.
fn tail(x: f64) -> f64 {
    if x >= 39.0 {
        return 0.0;
    }
    if x < 0.66291 {
        return 0.5 - small(x);
    }
    if x < 32.0_f64.sqrt() {
        const C: [f64; 9] = [
            0.39894151208813466764,
            8.8831497943883759412,
            93.506656132177855979,
            597.27027639480026226,
            2494.5375852903726711,
            6848.1904505362823326,
            11602.651437647350124,
            9842.7148383839780218,
            1.0765576773720192317e-8,
        ];
        const D: [f64; 8] = [
            22.266688044328115691,
            235.38790178262499861,
            1519.377599407554805,
            6485.558298266760755,
            18615.571640885098091,
            34900.952721145977266,
            38912.003286093271411,
            19685.429676859990727,
        ];
        let mut num = C[8] * x;
        let mut den = x;
        for i in 0..7 {
            num = (num + C[i]) * x;
            den = (den + D[i]) * x;
        }
        scaled_tail(x, (num + C[7]) / (den + D[7]))
    } else {
        const P: [f64; 6] = [
            0.21589853405795699,
            0.1274011611602473639,
            0.022235277870649807,
            0.001421619193227893466,
            2.9112874951168792e-5,
            0.02307344176494017303,
        ];
        const Q: [f64; 5] = [
            1.28426009614491121,
            0.468238212480865118,
            0.0659881378689285515,
            0.00378239633202758244,
            7.29751555083966205e-5,
        ];
        let square = 1.0 / (x * x);
        let mut num = P[5] * square;
        let mut den = square;
        for i in 0..4 {
            num = (num + P[i]) * square;
            den = (den + Q[i]) * square;
        }
        let rational = (0.39894228040143267794 - square * (num + P[4]) / (den + Q[4])) / x;
        scaled_tail(x, rational)
    }
}
/// Evaluates coefficient arrays `a` and `b` as a rational polynomial at `x`.
fn rational(a: &[f64], b: &[f64], x: f64) -> f64 {
    a.iter().rev().fold(0.0, |v, c| v * x + c) / b.iter().rev().fold(0.0, |v, c| v * x + c)
}
/// Returns the central standard-normal quantile for centered probability `q`.
fn inverse_small(q: f64) -> f64 {
    const A: [f64; 8] = [
        3.387132872796366608,
        133.14166789178437745,
        1971.5909503065514427,
        13731.693765509461125,
        45921.953931549871457,
        67265.770927008700853,
        33430.575583588128105,
        2509.0809287301226727,
    ];
    const B: [f64; 8] = [
        1.0,
        42.313330701600911252,
        687.1870074920579083,
        5394.1960214247511077,
        21213.794301586595867,
        39307.89580009271061,
        28729.085735721942674,
        5226.495278852854561,
    ];
    q * rational(&A, &B, 0.180625 - q * q)
}
/// Returns the intermediate-tail normal quantile magnitude for transformed probability `r`.
fn inverse_intermediate(r: f64) -> f64 {
    const A: [f64; 8] = [
        1.42343711074968357734,
        4.6303378461565452959,
        5.7694972214606914055,
        3.64784832476320460504,
        1.27045825245236838258,
        0.24178072517745061177,
        0.0227238449892691845833,
        7.7454501427834140764e-4,
    ];
    const B: [f64; 8] = [
        1.0,
        2.05319162663775882187,
        1.6763848301838038494,
        0.68976733498510000455,
        0.14810397642748007459,
        0.0151986665636164571966,
        5.475938084995344946e-4,
        1.05075007164441684324e-9,
    ];
    rational(&A, &B, r - 1.6)
}
/// Returns the extreme-tail normal quantile magnitude for transformed probability `r`.
fn inverse_tail(r: f64) -> f64 {
    const A: [f64; 8] = [
        6.6579046435011037772,
        5.4637849111641143699,
        1.7848265399172913358,
        0.29656057182850489123,
        0.026532189526576123093,
        0.0012426609473880784386,
        2.71155556874348757815e-5,
        2.01033439929228813265e-7,
    ];
    const B: [f64; 8] = [
        1.0,
        0.59983220655588793769,
        0.13692988092273580531,
        0.0148753612908506148525,
        7.868691311456132591e-4,
        1.8463183175100546818e-5,
        1.4215117583164458887e-7,
        2.04426310338993978564e-15,
    ];
    rational(&A, &B, r - 5.0)
}

/// Standard normal `Pr(X <= x)`. Infinite observations return the limiting probability.
pub fn standard_normal_cdf(x: f64) -> Result<f64, Error> {
    observation(x)?;
    Ok(if x <= 0.0 { tail(-x) } else { 1.0 - tail(x) })
}
/// Standard normal `Pr(X > x)`, retaining small upper tails.
pub fn standard_normal_survival(x: f64) -> Result<f64, Error> {
    standard_normal_cdf(-x)
}
/// Standard normal CDF inverse for `p` in `[0,1]`, with infinite endpoint limits.
pub fn standard_normal_quantile(p: f64) -> Result<f64, Error> {
    probability(p)?;
    if p == 0.0 {
        return Ok(f64::NEG_INFINITY);
    }
    if p == 1.0 {
        return Ok(f64::INFINITY);
    }
    let q = p - 0.5;
    if q.abs() <= 0.425 {
        return Ok(inverse_small(q));
    }
    let r = (-p.min(1.0 - p).ln()).sqrt();
    let value = if r <= 5.0 {
        inverse_intermediate(r)
    } else {
        inverse_tail(r)
    };
    Ok(if p < 0.5 { -value } else { value })
}
/// Standard normal survival inverse, evaluated without subtracting `p` from one.
pub fn standard_normal_inverse_survival(p: f64) -> Result<f64, Error> {
    Ok(-standard_normal_quantile(p)?)
}
/// Zero-mean normal CDF with finite positive standard deviation `sigma`.
pub fn normal_cdf(x: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    standard_normal_cdf(x / sigma)
}
/// Zero-mean normal survival with finite positive standard deviation `sigma`.
pub fn normal_survival(x: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    standard_normal_survival(x / sigma)
}
/// Zero-mean normal quantile with finite positive standard deviation `sigma`.
pub fn normal_quantile(p: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    finite_quantile(sigma * standard_normal_quantile(p)?, p)
}
/// Zero-mean normal survival inverse with finite positive standard deviation `sigma`.
pub fn normal_inverse_survival(p: f64, sigma: f64) -> Result<f64, Error> {
    positive(sigma)?;
    finite_quantile(sigma * standard_normal_inverse_survival(p)?, p)
}

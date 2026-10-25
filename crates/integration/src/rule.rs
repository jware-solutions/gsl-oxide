//! Private embedded rule with mathematically reproducible coefficients.

use crate::{finite, Error, Sum};

// P7 = (429*x^7 - 693*x^5 + 315*x^3 - 35*x)/16.
// The monic even Stieltjes polynomial E8 is determined by
// integral_-1^1 P7(x)*E8(x)*x^k dx = 0, k=1,3,5,7. Exact rational moments give
// E8 = x^8 - (36/17)x^6 + (7794/5491)x^4 - (202548/653429)x^2
//      + 52932681/4854324041.
// Rounded f64 nodes and weights follow from exact rational elimination,
// root bisection, and integration of Lagrange basis polynomials.
// Gauss nodes are roots of P7; added nodes are roots of E8. Moment and
// polynomial tests below check this construction.
// Mathematical basis: https://dlmf.nist.gov/3.5.iv and /3.5.v;
// Embedded extension: https://www.gnu.org/software/gsl/doc/html/integration.html
// ("Integrands without weight functions").
const X: [f64; 8] = [
    0.0,
    0.207_784_955_007_898_48,
    0.405_845_151_377_397_2,
    0.586_087_235_467_691_1,
    0.741_531_185_599_394_5,
    0.864_864_423_359_769_1,
    0.949_107_912_342_758_5,
    0.991_455_371_120_812_6,
];
const K: [f64; 8] = [
    0.209_482_141_084_727_82,
    0.204_432_940_075_298_89,
    0.190_350_578_064_785_42,
    0.169_004_726_639_267_9,
    0.140_653_259_715_525_92,
    0.104_790_010_322_250_19,
    0.063_092_092_629_978_56,
    0.022_935_322_010_529_224,
];
const G: [f64; 4] = [
    0.417_959_183_673_469_4,
    0.381_830_050_505_118_9,
    0.279_705_391_489_276_64,
    0.129_484_966_168_869_7,
];

pub(super) fn weighted(y: f64, half: f64, weight: f64) -> Result<f64, Error> {
    // Multiply by interval size first to preserve subnormal widths. Fall back
    // to weighting first when the unweighted intermediate would overflow.
    let scaled = y * half;
    finite(if scaled.is_finite() {
        scaled * weight
    } else {
        y * (half * weight)
    })
}

pub(super) struct Estimate {
    pub(super) value: f64,
    pub(super) error: f64,
    pub(super) absolute: f64,
}

pub(super) fn geometry(a: f64, b: f64) -> Result<(f64, f64), Error> {
    let width = b - a;
    let (mid, half) = if width.is_finite() {
        (a + width * 0.5, width * 0.5)
    } else {
        (a * 0.5 + b * 0.5, b * 0.5 - a * 0.5)
    };
    if !(a < mid && mid < b && half > 0.0) {
        Err(Error::Subdivision)
    } else {
        Ok((mid, half))
    }
}

pub(super) fn evaluate(
    f: &mut impl FnMut(f64) -> Result<f64, Error>,
    a: f64,
    b: f64,
) -> Result<Estimate, Error> {
    let (mid, half) = geometry(a, b)?;
    // Use the full width when representable: a subnormal width can be an odd
    // number of ulps, whose exact half is not representable. Keep its scale.
    let width = b - a;
    let (scale, weight_scale) = if width.is_finite() {
        (width, 0.5)
    } else {
        (half, 1.0)
    };
    let mut nodes = [(mid, mid); 8];
    let mut previous = (mid, mid);
    // Validate every mapped node before sampling. No endpoints or duplicate nodes.
    for i in 1..8 {
        let left = (-half).mul_add(X[i], mid);
        let right = half.mul_add(X[i], mid);
        if !(a < left && left < previous.0 && previous.1 < right && right < b) {
            return Err(Error::Subdivision);
        }
        nodes[i] = (left, right);
        previous = (left, right);
    }
    let center = finite(f(mid)?)?;
    let mut samples = [(center, center); 8];
    for i in 1..8 {
        samples[i] = (finite(f(nodes[i].0)?)?, finite(f(nodes[i].1)?)?);
    }
    let mut kronrod = Sum::default();
    let mut gauss = Sum::default();
    let mut absolute = Sum::default();
    let mut mean = Sum::default();
    for i in 0..8 {
        let count = if i == 0 { 1 } else { 2 };
        for y in [samples[i].0, samples[i].1].into_iter().take(count) {
            // Weight before adding a symmetric pair: f(left)+f(right) may overflow.
            let term = weighted(y, scale, K[i] * weight_scale)?;
            kronrod.add(term)?;
            absolute.add(term.abs())?;
            mean.add(y * (K[i] * 0.5))?;
            if i % 2 == 0 {
                gauss.add(weighted(y, scale, G[i / 2] * weight_scale)?)?;
            }
        }
    }
    let value = kronrod.total()?;
    let absolute = absolute.total()?;
    let mean = mean.total()?;
    let mut variation = Sum::default();
    for i in 0..8 {
        for y in [samples[i].0, samples[i].1]
            .into_iter()
            .take(if i == 0 { 1 } else { 2 })
        {
            variation.add(
                (weighted(y, scale, K[i] * weight_scale)?
                    - weighted(mean, scale, K[i] * weight_scale)?)
                .abs(),
            )?;
        }
    }
    let variation = variation.total()?;
    let difference = finite(value - gauss.total()?)?.abs();
    // A safety factor grows from 1 to 9 with sampled relative variation.
    // This heuristic inflates the embedded discrepancy without amplifying
    // roundoff-sized discrepancies into a spurious sqrt(epsilon) error floor.
    let factor = if absolute > 0.0 {
        1.0 + 8.0 * (variation / absolute).min(1.0).sqrt()
    } else {
        1.0
    };
    let scaled = finite(difference * factor)?;
    let resolution =
        ((f64::EPSILON * a.abs().max(b.abs())) / half + f64::from_bits(1) / half).min(1.0);
    // resolution already allows a full epsilon (twice the usual rounding unit).
    let mapping_error = finite(variation * resolution)?;
    let floor = (32.0 * f64::EPSILON * absolute)
        .max(mapping_error)
        // Absolute rounding of the 15 weighted samples also matters when all
        // terms are subnormal and epsilon*absolute has underflowed to zero.
        .max(16.0 * f64::from_bits(1));
    Ok(Estimate {
        value,
        error: finite(scaled.max(floor))?,
        absolute,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Budget;

    #[test]
    fn independent_polynomial_construction_and_moments() {
        assert!(X.windows(2).all(|w| w[0] < w[1]));
        assert!(K.iter().chain(G.iter()).all(|&w| w > 0.0));
        for (i, &x) in X.iter().enumerate().skip(1) {
            let z = x * x;
            let polynomial = if i % 2 == 0 {
                x * (((429.0 * z - 693.0) * z + 315.0) * z - 35.0) / 16.0
            } else {
                (((z - 36.0 / 17.0) * z + 7794.0 / 5491.0) * z - 202548.0 / 653429.0) * z
                    + 52932681.0 / 4854324041.0
            };
            assert!(polynomial.abs() < 5e-15, "root {i}: {polynomial}");
        }
        // Exactness: Gauss degree 13; Kronrod degree 23 (odd moments by symmetry).
        for degree in 0..=23 {
            let expected = if degree % 2 == 0 {
                2.0 / f64::from(degree + 1)
            } else {
                0.0
            };
            let mut km = K[0] * 0.0_f64.powi(degree);
            let mut gm = G[0] * 0.0_f64.powi(degree);
            for i in 1..8 {
                let moment = X[i].powi(degree) + (-X[i]).powi(degree);
                km += K[i] * moment;
                if i % 2 == 0 {
                    gm += G[i / 2] * moment;
                }
            }
            assert!((km - expected).abs() < 7e-16, "K degree {degree}");
            if degree <= 13 {
                assert!((gm - expected).abs() < 7e-16, "G degree {degree}");
            }
            let result = evaluate(&mut |x| Ok(x.powi(degree)), -1.0, 1.0).unwrap();
            assert!((result.value - expected).abs() < 7e-16);
            assert!(result.error.is_finite() && result.error >= 0.0);
        }
    }

    #[test]
    fn mapping_symmetry_reuse_and_roundoff_floor() {
        let mut nodes = Vec::new();
        let result = evaluate(
            &mut |x| {
                nodes.push(x);
                Ok((x - 4.0).powi(2))
            },
            2.0,
            6.0,
        )
        .unwrap();
        assert!((result.value - 16.0 / 3.0).abs() < 2e-15);
        assert_eq!(nodes.len(), 15);
        nodes.sort_by(f64::total_cmp);
        assert!(nodes.windows(2).all(|x| x[0] < x[1]));
        for i in 0..15 {
            assert!((nodes[i] + nodes[14 - i] - 8.0).abs() < 2e-15);
        }
        let result = evaluate(&mut |_| Ok(3.0), -2.0, 5.0).unwrap();
        assert!((result.value - 21.0).abs() < 4e-15);
        assert!(result.error > 0.0);
        // Normalization and zero/reversal belong to the public layer.
        for (a, b) in [
            (1.0, 1.0),
            (2.0, 1.0),
            (1.0, f64::from_bits(1.0_f64.to_bits() + 1)),
        ] {
            assert!(matches!(
                evaluate(&mut |_| panic!("no nodes"), a, b),
                Err(Error::Subdivision)
            ));
        }
    }

    #[test]
    fn arithmetic_and_sample_failures() {
        for y in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut budget = Budget::new(|_| y, 15);
            assert!(matches!(
                evaluate(&mut |x| budget.eval(x), 0.0, 1.0),
                Err(Error::NonFiniteIntegrand)
            ));
        }
        assert!(matches!(
            evaluate(&mut |_| Ok(f64::MAX), -1.0, 1.0),
            Err(Error::NumericalBreakdown)
        ));
        // Avoid overflow in pair sums, midpoint, and interval subtraction.
        let r = evaluate(&mut |_| Ok(1e308), 0.0, 1.0).unwrap();
        assert!((r.value / 1e308 - 1.0).abs() < 5e-16);
        let r = evaluate(&mut |_| Ok(1.0 / f64::MAX), -f64::MAX, f64::MAX).unwrap();
        assert!((r.value - 2.0).abs() < 2e-15);
        // Regression: weighting a subnormal width before the callback value
        // loses several digits, even when the final integral is normal.
        for bits in [10_000, 10_001] {
            let width = f64::from_bits(bits);
            let r = evaluate(&mut |_| Ok(1e308), 0.0, width).unwrap();
            assert!((r.value / (width * 1e308) - 1.0).abs() < 5e-16);
            let r = evaluate(&mut |x| Ok(x / width), 0.0, width).unwrap();
            assert!(
                r.error > width * 1e-6,
                "subnormal node mapping needs an error floor"
            );
        }
        let tiny = f64::from_bits(2);
        let r = evaluate(&mut |_| Ok(tiny), 0.0, 1.0).unwrap();
        assert!((r.value - tiny).abs() <= r.error);
    }
}

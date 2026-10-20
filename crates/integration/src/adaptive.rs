use crate::{finite, rule, target, Error, Sum};

struct Interval {
    a: f64,
    b: f64,
    part: u8,
    estimate: rule::Estimate,
    stalled: usize,
}

// Recompute compensated totals from active leaves, rather than subtracting a
// parent's error from an accumulating total. No negative bookkeeping residuals
// can hide outstanding error, and cancellation between segments remains visible.
pub(super) fn run(
    f: &mut impl FnMut(f64, u8) -> Result<f64, Error>,
    seeds: impl IntoIterator<Item = (f64, f64, u8)>,
    abs: f64,
    rel: f64,
) -> Result<(f64, f64), Error> {
    let mut intervals = Vec::new();
    for (a, b, part) in seeds {
        let estimate = rule::evaluate(&mut |x| f(x, part), a, b)?;
        intervals.push(Interval {
            a,
            b,
            part,
            estimate,
            stalled: 0,
        });
    }
    loop {
        let mut values = Sum::default();
        let mut errors = Sum::default();
        let mut absolute = Sum::default();
        let mut largest = 0;
        for (i, interval) in intervals.iter().enumerate() {
            values.add(interval.estimate.value)?;
            errors.add(interval.estimate.error)?;
            absolute.add(interval.estimate.absolute)?;
            if interval.estimate.error > intervals[largest].estimate.error {
                largest = i;
            }
        }
        let value = values.total()?;
        let error = errors.total()?;
        let goal = target(abs, rel, value);
        if error <= goal {
            return Ok((value, error));
        }
        if 32.0 * f64::EPSILON * absolute.total()? > goal {
            return Err(Error::Roundoff);
        }
        let parent = &intervals[largest];
        let (mid, _) = rule::geometry(parent.a, parent.b)?;
        let left = rule::evaluate(&mut |x| f(x, parent.part), parent.a, mid)?;
        let right = rule::evaluate(&mut |x| f(x, parent.part), mid, parent.b)?;
        let child_value = finite(left.value + right.value)?;
        let child_error = finite(left.error + right.error)?;
        let child_absolute = finite(left.absolute + right.absolute)?;
        let ineffective = child_error >= parent.estimate.error * 0.99
            && (child_value - parent.estimate.value).abs() <= 64.0 * f64::EPSILON * child_absolute;
        let stalled = if ineffective { parent.stalled + 1 } else { 0 };
        if stalled >= 8 {
            return Err(Error::Roundoff);
        }
        let (a, b, part) = (parent.a, parent.b, parent.part);
        intervals[largest] = Interval {
            a,
            b: mid,
            part,
            estimate: left,
            stalled,
        };
        intervals.push(Interval {
            a: mid,
            b,
            part,
            estimate: right,
            stalled,
        });
        // Every completed iteration spends 30 original callback calls (also on
        // transformed domains). Thus the callback budget bounds this loop too.
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Budget;

    #[test]
    fn repeated_ineffective_subdivisions_take_the_stagnation_path() {
        // Every 15-sample rule sees the same normalized pattern. Both children
        // reproduce the parent's result and error; the discrepancy never shrinks.
        // A controlled stateful callback exercises the actual stagnation policy.
        let mut calls = 0;
        let mut budget = Budget::new(
            |_| {
                let value = if calls % 15 == 0 { 0.0 } else { 1.0 };
                calls += 1;
                value
            },
            20_000,
        );
        let result = run(&mut |x, _| budget.eval(x), [(0.0, 1.0, 0)], 1e-10, 0.0);
        assert_eq!(result, Err(Error::Roundoff));
        assert!(budget.used > 30 && budget.used < 20_000);
    }
}

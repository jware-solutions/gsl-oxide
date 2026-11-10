// Adapted from GSL 2.8 sys/pow_int.c and gsl_pow_int.h.
// Copyright (C) 1996–2000, 2004, 2007 Gerard Jungman and Brian Gough.
// Distributed under GPL-3.0-or-later; see LICENSE.

/// Integer powers of an `f64`, with the same method signatures as `rgsl::Pow`.
///
/// Adapted from GSL 2.8's `gsl_pow_int`, `gsl_pow_uint`, and `gsl_pow_2`
/// through `gsl_pow_9`; see the [GSL reference]. General integer powers use
/// repeated squaring, requiring at most 32 iterations. Negative exponents
/// take the reciprocal of the base before exponentiation.
///
/// Like GSL, these methods return floating-point results directly without
/// checking overflow or underflow. An exponent of zero returns `1.0` for
/// every base, including zero, infinities, and NaN. For nonzero exponents,
/// NaN propagates; zeros and infinities follow IEEE 754 multiplication and
/// division, including their signs.
///
/// [GSL reference]: https://www.gnu.org/software/gsl/doc/html/math.html#small-integer-powers
///
/// # Examples
///
/// ```
/// use gsl_oxide_special::Pow as _;
///
/// assert_eq!(2.0_f64.pow_uint(3), 8.0);
/// assert_eq!(2.0_f64.pow_int(-3), 0.125);
/// assert_eq!((-3.0_f64).pow_2(), 9.0);
/// ```
pub trait Pow {
    /// Returns `self` raised to the signed integer power `n`.
    ///
    /// Supports the full `i32` range, including `i32::MIN`.
    #[doc(alias = "gsl_pow_int")]
    fn pow_int(&self, n: i32) -> Self;

    /// Returns `self` raised to the unsigned integer power `n`.
    ///
    /// Supports the full `u32` range without converting `n` to a signed integer.
    #[doc(alias = "gsl_pow_uint")]
    fn pow_uint(&self, n: u32) -> Self;

    /// Returns `self` squared.
    #[doc(alias = "gsl_pow_2")]
    fn pow_2(&self) -> Self;

    /// Returns `self` cubed.
    #[doc(alias = "gsl_pow_3")]
    fn pow_3(&self) -> Self;

    /// Returns `self` raised to the fourth power.
    #[doc(alias = "gsl_pow_4")]
    fn pow_4(&self) -> Self;

    /// Returns `self` raised to the fifth power.
    #[doc(alias = "gsl_pow_5")]
    fn pow_5(&self) -> Self;

    /// Returns `self` raised to the sixth power.
    #[doc(alias = "gsl_pow_6")]
    fn pow_6(&self) -> Self;

    /// Returns `self` raised to the seventh power.
    #[doc(alias = "gsl_pow_7")]
    fn pow_7(&self) -> Self;

    /// Returns `self` raised to the eighth power.
    #[doc(alias = "gsl_pow_8")]
    fn pow_8(&self) -> Self;

    /// Returns `self` raised to the ninth power.
    #[doc(alias = "gsl_pow_9")]
    fn pow_9(&self) -> Self;
}

impl Pow for f64 {
    /// Returns the base raised to signed integer exponent `n` using reciprocal-first exponentiation for negative exponents.
    #[inline]
    fn pow_int(&self, n: i32) -> Self {
        let base = if n < 0 { 1.0 / *self } else { *self };
        base.pow_uint(n.unsigned_abs())
    }

    /// Returns the base raised to unsigned integer exponent `n` by repeated squaring.
    #[inline]
    fn pow_uint(&self, mut n: u32) -> Self {
        let mut base = *self;
        let mut value = 1.0;
        loop {
            if n & 1 != 0 {
                value *= base;
            }
            n >>= 1;
            base *= base;
            if n == 0 {
                return value;
            }
        }
    }

    /// Returns the base squared using GSL's multiplication order.
    #[inline]
    fn pow_2(&self) -> Self {
        *self * *self
    }

    /// Returns the base cubed using GSL's multiplication order.
    #[inline]
    fn pow_3(&self) -> Self {
        *self * *self * *self
    }

    /// Returns the base raised to the 4th power using GSL’s multiplication order.
    #[inline]
    fn pow_4(&self) -> Self {
        let x2 = *self * *self;
        x2 * x2
    }

    /// Returns the base raised to the 5th power using GSL’s multiplication order.
    #[inline]
    fn pow_5(&self) -> Self {
        let x2 = *self * *self;
        x2 * x2 * *self
    }

    /// Returns the base raised to the 6th power using GSL’s multiplication order.
    #[inline]
    fn pow_6(&self) -> Self {
        let x2 = *self * *self;
        x2 * x2 * x2
    }

    /// Returns the base raised to the 7th power using GSL’s multiplication order.
    #[inline]
    fn pow_7(&self) -> Self {
        let x3 = *self * *self * *self;
        x3 * x3 * *self
    }

    /// Returns the base raised to the 8th power using GSL’s multiplication order.
    #[inline]
    fn pow_8(&self) -> Self {
        let x2 = *self * *self;
        let x4 = x2 * x2;
        x4 * x4
    }

    /// Returns the base raised to the 9th power using GSL’s multiplication order.
    #[inline]
    fn pow_9(&self) -> Self {
        let x3 = *self * *self * *self;
        x3 * x3 * x3
    }
}

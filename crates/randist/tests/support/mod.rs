use gsl_oxide_randist::RngCore;

/// Initial sample count shared by the GSL statistical regressions.
pub(crate) const N: usize = 100_000;

// Fixed SplitMix64 stream keeps statistical regressions reproducible and needs
// no OS entropy. This generator is only test infrastructure, not a public API.
/// Deterministic SplitMix64 random source for reproducible sampling tests.
pub(crate) struct Generator(
    /// Current 64-bit SplitMix state.
    u64,
);
impl Generator {
    /// Returns a test generator initialized to the fixed regression-test seed.
    ///
    /// Each call starts the same SplitMix64 stream, keeping statistical tests
    /// reproducible without accessing operating-system entropy.
    pub(crate) fn new() -> Self {
        Self(0x123456789abcdef0)
    }
}
impl RngCore for Generator {
    /// Advances the SplitMix64 state and returns the next 64 random bits.
    ///
    /// Wrapping arithmetic fixes the stream independently of debug overflow checks.
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    /// Advances the generator and returns the upper 32 bits of its next word.
    fn next_u32(&mut self) -> u32 {
        (self.next_u64() >> 32) as u32
    }
    /// Fills every byte of `dest` using successive generator words.
    ///
    /// Returns `()` after writing the entire slice, including an empty slice.
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        rand_core::impls::fill_bytes_via_next(self, dest);
    }
    /// Fills every byte of `dest` and returns `Ok(())`.
    ///
    /// The deterministic in-memory generator has no fallible entropy source.
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

/// Asserts that observed `x` is within `tolerance` of `expected` relatively.
///
/// The absolute bound is `tolerance * max(abs(expected), 1e-300)`.
/// Returns `()` on success and panics with both values when the bound is exceeded.
pub(crate) fn relative(x: f64, expected: f64, tolerance: f64) {
    assert!(
        (x - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{x} != {expected}"
    );
}

// Port of GSL randist/test.c testPDF: same 100 bins on [-5,5), N,
// five-sigma threshold and two-count guard. Like GSL, extend a trial on an
// outlier; cap at four extensions so a broken sampler fails promptly.
/// Checks samples against the continuous density callback `pdf`.
///
/// `sample` draws from a fixed-seed generator; `pdf` evaluates the reference
/// density. Integrates 100 bins over `[-5,5)` and delegates to
/// [`continuous_masses`]. Returns `()` on agreement and panics on failure.
pub(crate) fn continuous(mut sample: impl FnMut(&mut Generator) -> f64, pdf: impl Fn(f64) -> f64) {
    let mut p = [0.0; 100];
    for (i, p) in p.iter_mut().enumerate() {
        let a = if i == 50 { 0.0 } else { -5.0 + i as f64 * 0.1 };
        let b = if i == 49 {
            0.0
        } else {
            -5.0 + (i + 1) as f64 * 0.1
        };
        *p = gsl_oxide_integration::integrate(&pdf, a, b, 1e-13, 1e-7, 200_000)
            .unwrap()
            .0;
        assert!(p.is_finite() && *p >= 0.0, "invalid bin mass {i}: {p}");
    }
    continuous_masses(&mut sample, p);
}

/// Checks a continuous sampler against 100 expected bin probabilities `p`.
///
/// `sample` draws from a fixed-seed generator. Bins span `[-5,5)` with width
/// 0.1; probabilities outside that range are omitted. Uses GSL's five-sigma
/// and two-count guards, extending the initial 100,000 draws at most four
/// times for outliers. Returns `()` on agreement and panics otherwise.
pub(crate) fn continuous_masses(mut sample: impl FnMut(&mut Generator) -> f64, p: [f64; 100]) {
    let mut r = Generator::new();
    let mut counts = [0_usize; 100];
    let mut previous = 0;
    for attempt in 0..=4 {
        let n = N << attempt;
        let mut edges = [0_usize; 100];
        for _ in previous..n {
            let x = sample(&mut r);
            assert!(x.is_finite(), "nonfinite sample: {x}");
            if x > -5.0 && x < 5.0 {
                let u = (x + 5.0) / 0.1;
                let bin = (u.floor() as usize).min(99);
                if u.fract() == 0.0 {
                    edges[bin] += 1;
                } else {
                    counts[bin] += 1;
                }
            }
        }
        // GSL assigns rounded support endpoints to the occupied adjacent bin.
        for i in 0..100 {
            if i + 1 < 100 && counts[i + 1] == 0 {
                counts[i] += edges[i + 1];
                edges[i + 1] = 0;
            }
            counts[i] += edges[i];
        }
        let failed = (0..100).find(|&i| {
            let expected = n as f64 * p[i];
            let delta = (counts[i] as f64 - expected).abs();
            if expected > 0.0 {
                delta > 5.0 * expected.sqrt() && delta > 2.0
            } else {
                counts[i] != 0
            }
        });
        if let Some(i) = failed {
            assert!(
                attempt < 4,
                "bin {i}: {} observed, {} expected",
                counts[i],
                n as f64 * p[i]
            );
        } else {
            return;
        }
        previous = n;
    }
}

// Port of GSL testDiscretePDF, retaining N, 100 bins and five-sigma criterion.
/// Checks integer sampler `sample` against probability callback `pmf`.
///
/// Draws 100,000 values from a fixed seed and compares frequencies for counts
/// 0 through 99 with GSL's five-sigma and two-count criteria. Returns `()`
/// on agreement; invalid probabilities or discrepant frequencies panic.
pub(crate) fn discrete(mut sample: impl FnMut(&mut Generator) -> u32, pmf: impl Fn(u32) -> f64) {
    let mut r = Generator::new();
    let mut counts = [0_usize; 100];
    for _ in 0..N {
        let x = sample(&mut r) as usize;
        if x < 100 {
            counts[x] += 1;
        }
    }
    for (i, &observed) in counts.iter().enumerate() {
        let expected = N as f64 * pmf(i as u32);
        let delta = (observed as f64 - expected).abs();
        assert!(expected.is_finite() && expected >= 0.0);
        assert!(
            if expected > 0.0 {
                delta <= 5.0 * expected.sqrt() || delta <= 2.0
            } else {
                observed == 0
            },
            "bin {i}: {observed} observed, {expected} expected"
        );
    }
}

/// Checks that `sample` hits the open interval `(lower, upper)` as expected.
///
/// `probability` is the reference probability of that interval. Draws 100,000
/// values from a fixed seed and requires the observed count to be within
/// three Poisson standard deviations. Returns `()` on success and panics on failure.
pub(crate) fn interval(
    mut sample: impl FnMut(&mut Generator) -> f64,
    lower: f64,
    upper: f64,
    probability: f64,
) {
    let mut r = Generator::new();
    let observed = (0..N)
        .filter(|_| {
            let x = sample(&mut r);
            x > lower && x < upper
        })
        .count() as f64;
    let expected = N as f64 * probability;
    let sigma = if expected > 0.0 {
        (observed - expected).abs() / expected.sqrt()
    } else {
        observed.abs()
    };
    assert!(
        sigma <= 3.0,
        "{observed} observed vs {expected} expected ({sigma} sigma)"
    );
}

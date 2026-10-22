# gsl-oxide

Independent native Rust implementations of selected GNU Scientific Library (GSL) 2.8 numerical functionality: special functions, numerical integration, and statistics.

The workspace uses safe Rust, forbids unsafe code in every library, and verifies numerical behavior with ordinary Rust tests. No GSL installation, bindings, or native numerical dependencies are required. Only a selected subset is in scope; gsl-oxide is not a full GSL replacement.

## Crates and coverage

| Crate | Family | Status |
| --- | --- | --- |
| `gsl-oxide-special` | Error functions | `error_functions::{erf, erfc}` |
| | Gamma and signed log gamma | `gamma::{gamma, ln_gamma}` |
| | Beta and log beta | `beta::{beta, ln_beta}` for positive finite shapes |
| | Regularized incomplete gamma | `incomplete_gamma::{p, q}` |
| | Regularized incomplete beta | `incomplete_beta::{regularized, complement}` |
| | Cylindrical Bessel J/Y/I/K | `bessel::{j0, j1, jn, y0, y1, yn, i0, i1, i_n, k0, k1, kn}` and `_scaled` I/K forms |
| `gsl-oxide-integration` | Gauss-Kronrod quadrature | Private Gauss 7 / Kronrod 15 kernel |
| | Adaptive finite-interval integration | `integrate` |
| | Breakpoint-aware integration | `integrate_with_breakpoints` |
| | Infinite-interval integration | `integrate_upper_infinite`, `integrate_lower_infinite`, `integrate_infinite` |
| | Integrable endpoint singularities | `integrate_singular` (tanh–sinh; logarithms and powers with exponent ≥ −1/2) |
| `gsl-oxide-statistics` | Mean, variance, and standard deviation | `mean`, `population_variance`, `sample_variance`, `population_standard_deviation`, `sample_standard_deviation` |
| | Covariance and correlation | `sample_covariance`, `correlation` (Pearson) |
| | Weighted statistics | `weighted_mean`, `weighted_population_variance`, `weighted_reliability_variance` |
| | Streaming and mergeable statistics | `RunningStatistics` with checked `push`, `merge`, queries, and reset |

Each library builds independently and has no dependencies, including on the other workspace crates. Unlisted GSL functionality is unsupported. Statistics quantiles, medians, histograms, higher moments, regression, rank correlations, moving windows, streaming covariance, and weighted streaming remain outside the implemented subset. See the repository's [implementation plan](https://github.com/jware-solutions/gsl-oxide/blob/HEAD/IMPLEMENTATION_PLAN.md) for the work breakdown.

See the Rust API documentation for supported arguments, numerical limits, and error behavior. Build and open it locally with `cargo doc -p gsl-oxide-special --open`, `cargo doc -p gsl-oxide-integration --open`, or `cargo doc -p gsl-oxide-statistics --open`.

Statistics use checked `Result<f64, Error>` returns and allocation-free borrowed slices. Population variance divides the centered moment by `n`; sample variance and covariance use `n - 1`. Weighted population variance divides by total weight `W`; reliability variance uses GSL's `W - sum(w*w)/W` correction, including for integer weights. Finite nonnegative weights are normalized by their maximum; zero weight removes even a non-finite observation. Each positive normalized weight must be at least `f64::MIN_POSITIVE`. Correlation requires nonconstant inputs and uses scaled centered products.

Compensated sums and centered scaling avoid raw-moment cancellation and unnecessary overflow. Successful results are finite; empty inputs, undefined estimators, non-finite participating observations, unsupported weight ratios, and numerical/count overflow return errors. Ordinary floating-point underflow is allowed. Standard-deviation queries require a representable variance even when the mathematical deviation would fit. Streaming additionally requires representable mean differences and an unnormalized centered sum of squares; merge order affects rounding. Counts above `2^53` incur rounding when converted for numerical ratios. These are numerical algorithms, without a correct-rounding or cross-platform bitwise guarantee.

## Getting started

Requires **Rust 1.85 or newer**, using the Rust 2021 edition and Cargo resolver 2. The virtual workspace shares package metadata, lints, `Cargo.lock`, and the root `target/` directory.

Crates have not been published. Until publication, point your application's `Cargo.toml` at a local checkout (adjust the path):

```toml
[dependencies]
gsl-oxide-special = { path = "../gsl-oxide/crates/special" }
gsl-oxide-integration = { path = "../gsl-oxide/crates/integration" }
gsl-oxide-statistics = { path = "../gsl-oxide/crates/statistics" }
```

```rust
use gsl_oxide_special::error_functions::erf;

fn main() {
    assert!((erf(1.0) - 0.8427007929497149).abs() < 2.0e-15);
}
```

```rust
use gsl_oxide_special::{bessel, gamma, incomplete_gamma};

fn main() -> Result<(), gsl_oxide_special::Error> {
    let (log_abs, sign) = gamma::ln_gamma(-0.5)?;
    assert!(log_abs > 1.0 && sign == -1);
    assert!(incomplete_gamma::q(1.0, 100.0)? > 0.0);
    assert!(bessel::i0_scaled(1000.0)? > 0.01);
    Ok(())
}
```

```rust
use gsl_oxide_integration::integrate;

fn main() -> Result<(), gsl_oxide_integration::Error> {
    let (area, estimated_error) = integrate(
        f64::sin, 0.0, std::f64::consts::PI, 1e-12, 1e-12, 10_000,
    )?;
    assert!((area - 2.0).abs() < 1e-12);
    assert!(estimated_error <= 2e-12);
    Ok(())
}
```

```rust
use gsl_oxide_statistics::{mean, sample_variance, weighted_reliability_variance, RunningStatistics};

fn main() -> Result<(), gsl_oxide_statistics::Error> {
    assert_eq!(mean(&[1., 2., 3., 4., 5.])?, 3.);
    assert_eq!(sample_variance(&[1., 2., 3., 4., 5.])?, 2.5);
    assert!((weighted_reliability_variance(&[1., 4.], &[1., 3.])? - 4.5).abs() < 1e-14);
    let mut left = RunningStatistics::new();
    left.push(1.)?;
    let mut right = RunningStatistics::new();
    for x in [2., 3., 4., 5.] { right.push(x)?; }
    left.merge(&right)?;
    assert_eq!(left.count(), 5);
    assert_eq!(left.sample_variance()?, 2.5);
    Ok(())
}
```

Applications will be able to depend on each published crate separately. From this checkout, run `cargo test --workspace`, or `cargo test -p gsl-oxide-special`, `cargo test -p gsl-oxide-integration`, or `cargo test -p gsl-oxide-statistics` for one crate. Each library can also be checked independently with `cargo check -p <package>`. Statistics examples are checked with `cargo test -p gsl-oxide-statistics --doc`.

## Alternatives for modules outside this workspace

These domains are deliberately excluded because established Rust libraries already address them. Applications should use those libraries directly. gsl-oxide provides no wrappers, adapters, conversions, or compatibility layers for them.

| Area outside this workspace's scope | Existing alternatives |
| --- | --- |
| N-dimensional arrays and array operations | [ndarray](https://docs.rs/ndarray/latest/ndarray/) |
| Linear algebra and matrix decompositions | [nalgebra](https://docs.rs/nalgebra/latest/nalgebra/), [faer](https://docs.rs/faer/latest/faer/) |
| FFT | [RustFFT](https://docs.rs/rustfft/latest/rustfft/) |
| Optimization | [argmin](https://docs.rs/argmin/latest/argmin/) |
| Random-number generation and sampling | [rand](https://docs.rs/rand/latest/rand/) |
| Statistical distributions | [statrs](https://docs.rs/statrs/latest/statrs/), [rand_distr](https://docs.rs/rand_distr/latest/rand_distr/) |

## Contributing

Contributions are welcome from the community. Please follow the guidelines in [CONTRIBUTING.md](CONTRIBUTING.md).

## Licensing status

The final license for the new Rust implementation is still under consideration. MIT is preferred, contingent on confirming the necessary rights and provenance; no MIT-licensed release has been declared. Independent development preserves that option, and Cargo license metadata remains unset pending a final decision.

The existing root `LICENSE` is preserved unchanged. Its current text is MIT; this does not finalize the new Rust implementation's licensing.

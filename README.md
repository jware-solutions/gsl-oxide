# gsl-oxide

Native Rust implementations of selected GNU Scientific Library (GSL) 2.8 numerical functionality: special functions, numerical integration, and statistics. Released under **GPL-3.0-or-later**.

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

Special functions use `f64`; `erf`/`erfc` return floating-point values directly (including NaN and infinity limits), while the other families return checked results. Gamma accepts finite nonpole real arguments; beta requires positive finite shapes. Incomplete gamma supports `0 < a <= 1e6`, `x >= 0` including positive infinity; incomplete beta supports finite `0 < a,b <= 1e6`, `0 <= x <= 1`. Complementary tails are computed directly where needed, but may underflow. The largest balanced incomplete-beta shapes can lose about `1e-9` absolute accuracy. Bessel orders are nonnegative integers `0..=1000`; finite arguments satisfy `|x| <= 1e6` for J/I and `0 < x <= 1e6` for Y/K, including scaled variants. Y/K reject zero. Scaled I/K avoid unscaled overflow/underflow; relative accuracy deteriorates near J/Y zeros and recurrence roundoff grows with order. Unsupported ranges and convergence failures are explicit errors.

Integration returns `(integral, estimated_absolute_error)` inside `Result`. Estimates are numerical diagnostics, not rigorous bounds for arbitrary callbacks; finite sampling can miss narrow features. Breakpoints help with known interior difficulties. Tolerances are finite and nonnegative, with at least one positive; pure relative requests require at least `64 * f64::EPSILON`. All segments share a hard original-callback evaluation budget. Infinite integration requires independently convergent tails and sufficiently regular transformed integrands; principal values and specialized oscillatory methods are unsupported. Endpoint integration supports interior-smooth functions with logarithms or algebraic endpoint powers `distance^alpha`, `alpha >= -1/2`, and never evaluates endpoints. Rounded callback arguments near nonzero endpoints can limit inverse-square-root accuracy to roughly `1e-7` times the integral's scale. Unresolved tails, stagnation, non-finite samples, and exhausted budgets return errors.

Statistics use checked `Result<f64, Error>` returns and allocation-free borrowed slices. Population variance divides the centered moment by `n`; sample variance and covariance use `n - 1`. Weighted population variance divides by total weight `W`; reliability variance uses GSL's `W - sum(w*w)/W` correction, including for integer weights. Finite nonnegative weights are normalized by their maximum; zero weight removes even a non-finite observation. Each positive normalized weight must be at least `f64::MIN_POSITIVE`. Correlation requires nonconstant inputs and uses scaled centered products.

Compensated sums and centered scaling avoid raw-moment cancellation and unnecessary overflow. Successful results are finite; empty inputs, undefined estimators, non-finite participating observations, unsupported weight ratios, and numerical/count overflow return errors. Ordinary floating-point underflow is allowed. Standard-deviation queries require a representable variance even when the mathematical deviation would fit. Streaming additionally requires representable mean differences and an unnormalized centered sum of squares; merge order affects rounding. Counts above `2^53` incur rounding when converted for numerical ratios. These are numerical algorithms, without a correct-rounding or cross-platform bitwise guarantee.

## Getting started

Requires **Rust 1.85.0 or newer**, using the Rust 2021 edition and Cargo resolver 2. The virtual workspace shares package metadata, lints, `Cargo.lock`, and the root `target/` directory. Ordinary platform linker/toolchain requirements still apply; GSL and a C/C++ numerical harness are unnecessary.

CI is configured to execute debug/release tests, debug/release doctests, documentation builds, and separate package builds on Linux x86_64 (`ubuntu-24.04`), Windows x86_64 (`windows-2022`), and native macOS ARM64 (`macos-15`). Linux also tests exactly Rust 1.85.0. These are the intended supported release targets; writing the workflow does not establish cross-platform success. Local Phase 6 tests, doctests and independent builds pass on Windows x86_64 with Rust 1.99.1 and exactly 1.85.0; the first hosted runs must confirm the matrix, including native ARM64, before release. Runner access and Actions billing must be enabled for this private repository.

Crates have not been published. Until publication, point your application's `Cargo.toml` at a local checkout (adjust the path):

```toml
[dependencies]
gsl-oxide-special = { path = "../gsl-oxide/crates/special" }
gsl-oxide-integration = { path = "../gsl-oxide/crates/integration" }
gsl-oxide-statistics = { path = "../gsl-oxide/crates/statistics" }
```

The initial shared registry version is **0.1.0**, pending successful CI validation and publication. After the required crates are published, add only the dependencies your application uses:

```toml
[dependencies]
gsl-oxide-special = "0.1.0"
gsl-oxide-integration = "0.1.0"
gsl-oxide-statistics = "0.1.0"
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

Applications will be able to depend on each published crate separately. From this checkout, run `cargo test --workspace --locked`, or `cargo test -p gsl-oxide-special --locked`, `cargo test -p gsl-oxide-integration --locked`, or `cargo test -p gsl-oxide-statistics --locked` for one crate. Each library can also be built independently with `cargo build -p <package> --locked`. Equivalent examples live in each crate's documentation and run with `cargo test --workspace --doc --locked`; root README code blocks are not themselves doctested.

## CI and release

[CI](https://github.com/jware-solutions/gsl-oxide/actions/workflows/ci.yml) runs on pull requests and pushes to `main`, with read-only permissions and no publication token. It checks formatting, Clippy, debug/release tests and doctests, documentation, independent builds, and Cargo package verification. See [CONTRIBUTING.md](https://github.com/jware-solutions/gsl-oxide/blob/main/CONTRIBUTING.md) for contribution instructions.

The owner triggers [Release to crates.io](https://github.com/jware-solutions/gsl-oxide/actions/workflows/release.yml) **from `main`**, with `dry_run` defaulting to `true`. It reuses CI for the dispatch commit, then uses [release-plz](https://release-plz.dev/docs/usage/release) to dry-run all unpublished workspace versions. Setting `dry_run=false` also publishes them. Release-plz handles registry checks, skips versions already published, and waits for registry availability. Versions stay shared, initially `0.1.0`; this workflow does not bump versions, create release PRs, tags, or GitHub Releases. There are no tag/package inputs or custom release scripts. Concurrent release runs cannot cancel an active publication.

Create the GitHub Actions secret CARGO_REGISTRY_TOKEN using your crates.io API token. No username is required.

The token is exposed only to the real publishing step. Dry-runs need no registry token. Confirm the published versions and test registry dependencies before declaring a release complete; docs.rs generation is asynchronous. **P6.7 remains incomplete until hosted validation, publication, and registry verification succeed.**

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

Contributions are welcome from the community. Please follow the guidelines in [CONTRIBUTING.md](https://github.com/jware-solutions/gsl-oxide/blob/main/CONTRIBUTING.md).

## License

This project derives from [GNU GSL](https://www.gnu.org/software/gsl/) and is distributed under **GPLv3 or later** (`GPL-3.0-or-later`); see [LICENSE](https://github.com/jware-solutions/gsl-oxide/blob/main/LICENSE).

Rust adaptation with revised numerical methods and APIs: Copyright (c) 2026 JWare Solutions.

GNU GSL attribution: Copyright (C) 1996–2000, 2007 Brian Gough; Copyright (C) 1996–2000 Gerard Jungman, and other GSL contributors.

Distributed **without any warranty**, including merchantability or fitness for a particular purpose.

# gsl-oxide

Native Rust implementations of selected GNU Scientific Library (GSL) 2.8 numerical functionality: special functions, numerical integration, and statistics. Released under **GPL-3.0-or-later**.

The workspace uses safe Rust, forbids unsafe code in every library, and verifies numerical behavior with ordinary Rust tests. No GSL installation, bindings, or native numerical dependencies are required. Only a selected subset is in scope; gsl-oxide is not a full GSL replacement.

## Crates and coverage

| Crate | Supported functionality |
| --- | --- |
| `gsl-oxide-special` | Error functions; gamma and beta functions, including logarithms and regularized incomplete forms; integer-order Bessel J/Y/I/K functions with scaled I/K variants; real exponential integrals with scaled variants; Legendre and generalized Laguerre polynomials. |
| `gsl-oxide-integration` | Adaptive integration over finite and infinite intervals, user-supplied breakpoints, and endpoint singularities involving logarithms or powers with exponent ≥ −1/2. |
| `gsl-oxide-statistics` | Mean, population/sample variance and standard deviation, sample covariance, Pearson and Spearman correlation, weighted mean and variance, and streaming statistics with merge support. |

Unlisted GSL functionality is unsupported. Consult the Rust API documentation for supported arguments, numerical limits, and error behavior. Generate it locally with `cargo doc --workspace --no-deps`.

See [GSL function name equivalencies](FUNCTION_EQUIVALENCIES.md) for a reference when transitioning from GSL. The guide maps GSL functions to this workspace's Rust API and lists previous Rust names that were replaced to make their purpose clearer.

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
    let (log_abs, sign) = gamma::ln_abs_gamma_with_sign(-0.5)?;
    assert!(log_abs > 1.0 && sign == -1);
    assert!(incomplete_gamma::upper_regularized(1.0, 100.0)? > 0.0);
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
| Kendall rank correlation | [kendalls](https://docs.rs/kendalls/latest/kendalls/) |

## Contributing

Contributions are welcome from the community. Please follow the guidelines in [CONTRIBUTING.md](https://github.com/jware-solutions/gsl-oxide/blob/main/CONTRIBUTING.md).

## License

This project derives from [GNU GSL](https://www.gnu.org/software/gsl/) and is distributed under **GPLv3 or later** (`GPL-3.0-or-later`); see [LICENSE](https://github.com/jware-solutions/gsl-oxide/blob/main/LICENSE).

Rust adaptation with revised numerical methods and APIs: Copyright (c) 2026 JWare Solutions.

GNU GSL attribution: Copyright (C) 1996–2000, 2007 Brian Gough; Copyright (C) 1996–2000 Gerard Jungman, and other GSL contributors.

Distributed **without any warranty**, including merchantability or fitness for a particular purpose.

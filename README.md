# gsl-oxide

Native Rust implementations of selected GNU Scientific Library (GSL) 2.8 numerical functionality: special functions, numerical integration, and statistics. Released under **GPL-3.0-or-later**.

The libraries use safe Rust and require no GSL installation or native numerical dependencies. Only a selected subset is supported; gsl-oxide is not a full GSL replacement.

## Crates and coverage

| Crate | Supported functionality |
| --- | --- |
| [gsl-oxide-special](https://docs.rs/gsl-oxide-special) | Error, gamma, beta, Bessel, and exponential integral functions; Legendre and generalized Laguerre polynomials; integer powers via `Pow`. |
| [gsl-oxide-integration](https://docs.rs/gsl-oxide-integration) | Adaptive integration over finite and infinite intervals, with breakpoints and supported endpoint singularities. |
| [gsl-oxide-statistics](https://docs.rs/gsl-oxide-statistics) | Mean, variance, standard deviation, covariance, correlations, weighted statistics, and streaming statistics. |

See each crate's API documentation for supported arguments, numerical limits, and error behavior.

See [GSL function name equivalencies](FUNCTION_EQUIVALENCIES.md) when transitioning from GSL.

## Getting started

Requires **Rust 1.85.0 or newer**. Add the crates you need from crates.io to your `Cargo.toml`:

```toml
[dependencies]
gsl-oxide-special = "0.1.0"
gsl-oxide-integration = "0.1.0"
gsl-oxide-statistics = "0.1.0"
```

```rust
use gsl_oxide_integration::integrate;
use gsl_oxide_special::error_functions::erf;
use gsl_oxide_statistics::mean;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("erf(1) = {}", erf(1.0));
    println!("Mean = {}", mean(&[1.0, 2.0, 3.0])?);

    let (area, _) = integrate(|x| x * x, 0.0, 1.0, 1e-8, 1e-8, 1_000)?;
    println!("Integral of x² from 0 to 1 = {area}");

    Ok(())
}
```

## Alternatives for modules outside this workspace

For functionality outside this workspace's scope, use these libraries directly:

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

### CI and release

Every push to `main` runs the release workflow, which validates the full CI matrix before publishing unpublished crate versions to crates.io using [release-plz](https://release-plz.dev/). Each published crate also gets a `<package>-v<version>` tag and a GitHub Release with automatically generated release notes (for example, `gsl-oxide-special-v0.1.0`).

The repository must have a `CARGO_REGISTRY_TOKEN` Actions secret with publishing access to the crates. Tags and GitHub Releases use the automatic `GITHUB_TOKEN`, with `contents: write` permission granted to the publish job.

Versions are taken from `Cargo.toml`. Versions already on crates.io are skipped, including tag and GitHub Release creation.

## License

This project derives from [GNU GSL](https://www.gnu.org/software/gsl/) and is distributed under **GPLv3 or later** (`GPL-3.0-or-later`); see [LICENSE](https://github.com/jware-solutions/gsl-oxide/blob/main/LICENSE).

Rust adaptation with revised numerical methods and APIs: Copyright (c) 2026 JWare Solutions.

GNU GSL attribution: Copyright (C) 1996–2000, 2007 Brian Gough; Copyright (C) 1996–2000 Gerard Jungman, and other GSL contributors.

Distributed **without any warranty**, including merchantability or fitness for a particular purpose.

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
| `gsl-oxide-statistics` | Mean, variance, and standard deviation | Planned |
| | Covariance and correlation | Planned |
| | Weighted statistics | Planned |
| | Streaming and mergeable statistics | Planned |

Each library builds independently and has no dependencies, including on the other workspace crates. Statistics currently exposes no numerical routines. Planned functionality is not yet available; unlisted GSL functionality is unsupported. See the repository's [implementation plan](https://github.com/jware-solutions/gsl-oxide/blob/HEAD/IMPLEMENTATION_PLAN.md) for the work breakdown.

See the Rust API documentation for supported arguments, numerical limits, and error behavior. Build and open it locally with `cargo doc -p gsl-oxide-special --open` or `cargo doc -p gsl-oxide-integration --open`.

## Getting started

Requires **Rust 1.85 or newer**, using the Rust 2021 edition and Cargo resolver 2. The virtual workspace shares package metadata, lints, `Cargo.lock`, and the root `target/` directory.

Crates have not been published. Until publication, point your application's `Cargo.toml` at a local checkout (adjust the path):

```toml
[dependencies]
gsl-oxide-special = { path = "../gsl-oxide/crates/special" }
gsl-oxide-integration = { path = "../gsl-oxide/crates/integration" }
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

Applications will be able to depend on each published crate separately. From this checkout, run `cargo test --workspace`, or `cargo test -p gsl-oxide-special` / `cargo test -p gsl-oxide-integration` for one crate. Each library can also be checked independently with `cargo check -p <package>`.

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

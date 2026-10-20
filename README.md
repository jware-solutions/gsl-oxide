# gsl-oxide

Independent native Rust implementations of selected GNU Scientific Library (GSL) 2.8 numerical functionality: special functions, numerical integration, and statistics.

The workspace uses safe Rust, forbids unsafe code in every library, and verifies numerical behavior with ordinary Rust tests. No GSL installation, bindings, or native numerical dependencies are required. Only a selected subset is in scope; gsl-oxide is not a full GSL replacement.

## Crates and coverage

| Crate | Family | Status |
| --- | --- | --- |
| `gsl-oxide-special` | Error functions | `error_functions::erf` implemented and tested; `erfc` planned |
| | Gamma and log gamma | Planned |
| | Beta and log beta | Planned |
| | Regularized incomplete gamma | Planned |
| | Regularized incomplete beta | Planned |
| | Bessel J/Y/I/K, including scaled variants | Planned |
| `gsl-oxide-integration` | Gauss-Kronrod quadrature | Planned |
| | Adaptive finite-interval integration | Planned |
| | Breakpoint-aware integration | Planned |
| | Infinite-interval integration | Planned |
| | Integrable endpoint singularities | Planned |
| `gsl-oxide-statistics` | Mean, variance, and standard deviation | Planned |
| | Covariance and correlation | Planned |
| | Weighted statistics | Planned |
| | Streaming and mergeable statistics | Planned |

Each library builds independently and has no dependencies, including on the other workspace crates. Integration and statistics currently expose no numerical routines. Planned functionality is not yet available; unlisted GSL functionality is unsupported. See the repository's [implementation plan](https://github.com/jware-solutions/gsl-oxide/blob/HEAD/IMPLEMENTATION_PLAN.md) for the work breakdown.

APIs use standard Rust types and conventional mathematical names. `erf(f64) -> f64` uses ordinary floating-point semantics: NaN propagates, infinities return their signed unit limits, and signed zero is preserved. Checked operations will use `Result` when explicit failure handling is needed.

## Getting started

Requires **Rust 1.85 or newer**, using the Rust 2021 edition and Cargo resolver 2. The virtual workspace shares package metadata, lints, `Cargo.lock`, and the root `target/` directory.

Crates have not been published. Until publication, point your application's `Cargo.toml` at a local checkout (adjust the path):

```toml
[dependencies]
gsl-oxide-special = { path = "../gsl-oxide/crates/special" }
```

```rust
use gsl_oxide_special::error_functions::erf;

fn main() {
    assert!((erf(1.0) - 0.8427007929497149).abs() < 2.0e-15);
}
```

Applications will be able to depend on each published crate separately. From this checkout, run `cargo test --workspace`, or `cargo test -p gsl-oxide-special` for a single crate. `cargo check -p gsl-oxide-integration` and `cargo check -p gsl-oxide-statistics` independently check the remaining libraries.

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

## References and testing

The public [GSL 2.8 definitions and documented behavior](https://www.gnu.org/software/gsl/doc/html/) are the behavior reference. Implementations and tests are independently written from mathematical definitions and established numerical methods, without copying or translating GSL implementation code or its coefficient tables. Short mathematical provenance notes live beside the Rust source.

Unit, integration, and documentation tests run through Cargo with expected values stored directly in Rust, justified floating-point tolerances, and explicit edge cases. Tests require no external scripts, numerical fixtures, GSL, or C/C++ numerical compiler.

Reference cases specify both absolute and relative tolerances. Finite comparisons validate the tolerances and use `abs(actual - expected) <= max(abs_tol, rel_tol * abs(expected))`; NaN, infinity, and signed-zero behavior use direct assertions, including zero sign bits. Deterministic identity checks supplement the reference values, with explicit coverage of subnormals and adjacent representable branch-boundary inputs.

## Licensing status

The final license for the new Rust implementation is still under consideration. MIT is preferred, contingent on confirming the necessary rights and provenance; no MIT-licensed release has been declared. Independent development preserves that option, and Cargo license metadata remains unset pending a final decision.

The existing root `LICENSE` is preserved unchanged. Its current text is MIT; this does not finalize the new Rust implementation's licensing.

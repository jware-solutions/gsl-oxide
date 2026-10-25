# gsl oxide implementation plan

Revised on October 25, 2026 (GPL licensing decision; implementation scope unchanged).

Implement the supported numerical behavior of GSL 2.8 in a small Rust workspace: special functions, numerical integration, and selected statistics. Publish each area as a separately usable crate. Implement the routines independently from public mathematical definitions and documented behavior, using short references instead of separate specifications or research documents. Verify the implementation with ordinary Rust tests and committed expected values.

The project adopted **GPL-3.0-or-later** on October 25, 2026. The root `LICENSE` now contains the full GPLv3 text; before this decision it contained MIT, contrary to this plan's earlier assertion. The historical MIT notice is retained in `NOTICE`, without revoking permissions validly granted for earlier versions. All crates inherit GPL metadata and distribute the license and notices. Keep accurate mathematical and third-party provenance, including original notices and dated modifications when adapting GPL material. The upstream `gsl-2.8/` reference retains its own licenses, including separate documentation terms, and stays outside published packages.

## 1 Implementation rules

- Use three workspace members: `gsl-oxide-special`, `gsl-oxide-integration`, and `gsl-oxide-statistics`. Users can depend on just the package they need. Keep family-level functions as modules within these packages; do not create a crate per function, a shared core, an umbrella facade, or a plugin system.
- Share the Cargo lockfile, build output, package metadata, and lint settings through the workspace. Each member has its own public API, tests, and package release. Keep the three crates independent initially; add a dependency between them only when an implementation actually needs it. Cargo's [workspace documentation](https://doc.rust-lang.org/cargo/reference/workspaces.html) describes the shared configuration.
- Implement the numerical routines in Rust. Do not add GSL bindings, FFI, native GSL dependencies, or C/C++ build and comparison programs, including in tests.
- Use the GSL 2.8 manual for mathematical definitions and expected behavior. Prefer original mathematical implementations. Any adaptation of compatible GPL source must preserve original notices, identify its origin, and record dated changes; manual text requires its own license review. Building, executing, or linking the reference source is unnecessary for this plan.
- Mention existing libraries for other domains in `README.md`. Do not implement adapters, wrappers, conversions, or integration layers for those libraries.
- Prefer standard Rust types. Add a custom type only when it represents necessary state or failures that standard types cannot express clearly.
- Keep crate and module scope, implementation status, and alternatives in the root README. Do not create a `selection/` directory, a separate scope document, or a coverage database.
- Keep documentation to the README and concise Rust documentation comments: purpose, inputs, outputs, important edge cases, and a GSL reference. Separate papers, derivations, algorithm documents, and precision-report infrastructure are unnecessary.
- Use `cargo test --workspace` as the test harness, or `cargo test -p <package>` for one crate. Expected results and tolerances belong in Rust tests; external scripts, fixture generators, GSL executables, and fuzzing infrastructure are outside the plan.

## 2 Native types and minimal API

| Need | Preferred representation |
| --- | --- |
| Numerical values and tolerances | `f64` |
| Counts, lengths, evaluation limits | `usize`; use an appropriate integer type for mathematical orders |
| Borrowed data | `&[f64]`, `&mut [f64]` |
| Fixed coefficients and small fixed collections | Rust arrays |
| Owned buffers | `Vec<f64>` or `Vec` of simple tuples |
| Integrands | `FnMut(f64) -> f64` |
| Value with estimated absolute error | `(f64, f64)` with documented tuple order |
| Log gamma with sign | `(f64, i32)` with documented tuple order |
| Checked operations | `Result<T, Error>` using a small error enum local to the relevant crate when needed |

Return `f64` directly where the function's documented floating-point behavior is sufficient. Use `Result` where invalid inputs or failure to converge need explicit handling. Do not add wrapper types such as `Real`, `Tolerance`, `AccuracyEstimate`, or crate-specific vector and matrix types.

Start with normal local buffers for integration and slice-based statistics. Add a reusable workspace or a small streaming-state struct only when the corresponding capability requires it. Such types should hold actual state, rather than wrap existing Rust types without a practical benefit.

## 3 Proposed layout

```text
gsl-oxide/
  Cargo.toml           # virtual workspace manifest
  Cargo.lock
  README.md
  LICENSE
  crates/
    special/           # gsl-oxide-special
      Cargo.toml
      src/lib.rs
      tests/special.rs
    integration/       # gsl-oxide-integration
      Cargo.toml
      src/lib.rs
      tests/integration.rs
    statistics/        # gsl-oxide-statistics
      Cargo.toml
      src/lib.rs
      tests/statistics.rs
  .github/workflows/
    ci.yml
  gsl-2.8/              # upstream reference source
```

Put unit tests next to private implementation details and public API tests in each member's `tests/` directory. Reuse private helpers within a crate when needed. Consider shared implementation infrastructure only after substantial real reuse appears; small local helpers do not justify another crate. Split a large module into family files when that improves readability. The directories above are proposed; this planning edit does not create the Rust implementation.

## 4 Phase 1 Set scope in the README and create the workspace

- [ ] **P1.1** Maintain a root README section listing the three crates and the modules planned and implemented in each, with their current status. Keep planned features clearly distinguished from completed ones.
- [ ] **P1.2** Maintain a README section linking existing alternatives for arrays, linear algebra, FFT, optimization, random sampling, and distributions. State that applications use those libraries directly.
- [ ] **P1.3** Create a virtual Cargo workspace with the three library members. Set an explicit resolver compatible with the chosen edition and minimum Rust version. Share common metadata and lints through workspace inheritance; use only license metadata supported by the project's actual licensing decision. Give each package a description and README that is present in its published archive. Expose only working implementations.
- [ ] **P1.4** Keep dependencies minimal. Use standard Rust facilities first. If an elementary numerical operation needs a dependency to retain a Rust-native calculation path, choose an existing Rust implementation rather than reimplementing it or linking a C numerical library.
- [ ] **P1.5** Add `#![forbid(unsafe_code)]` to each crate's project-owned implementation code. Do not add native build scripts or bindings.
- [ ] **P1.6** Establish naming, direct-return versus checked-return behavior, handling of NaN/infinity, and the tuple conventions through the first function and its tests. Explain these briefly in Rust documentation comments.

**Done when:** all workspace members build, the first numerical test runs, and the root README describes the actual scope and status without extra scope-management files. Each crate can build on its own without pulling in the other two.

## 5 Phase 2 Establish ordinary Rust tests

Start these tests with the first function and extend them alongside each subsequent implementation.

- [ ] **P2.1** Write table-driven Rust tests containing inputs, expected results, and per-case absolute or relative tolerances in the relevant crate. Use documented numerical examples and known mathematical results with short source references. Write the test code independently; review permissions before copying an upstream test suite or other protected material.
- [ ] **P2.2** Add a small test-only comparison helper if repeated assertions need it. For finite values, use an appropriate absolute/relative bound, for example `abs(actual - expected) <= max(abs_tol, rel_tol * abs(expected))`. Handle NaN, infinity, and signed zero explicitly where relevant.
- [ ] **P2.3** Match GSL's numerical results for the supported operations and domains at the precision required by their tests. Use exact assertions where meaningful and justified tolerances for rounded floating-point results. Do not treat arbitrary cross-platform bitwise identity as an implicit requirement.
- [ ] **P2.4** Test normal inputs, endpoints, poles, near-zero values, large magnitudes, overflow/underflow, and invalid inputs where applicable. Use adjacent representable values when testing a numerical branch boundary.
- [ ] **P2.5** Add simple deterministic identity tests where useful: gamma recurrence, complementary probability tails, Bessel recurrences, known integrals, and batch-versus-streaming statistics. Keep these as ordinary Rust loops and assertions.
- [ ] **P2.6** Add explicit failure tests for invalid tolerances, invalid domains, non-finite callbacks, and iteration/evaluation limits. Keep a regression case for every numerical bug fixed.
- [ ] **P2.7** Keep tests runnable with `cargo test --workspace` and `cargo test -p <package>` on a machine without GSL or a C/C++ compiler. Store expected values directly in Rust rather than introducing a JSON/CSV fixture system or an external validation pipeline.

**Done when:** the first implemented function has sourced expected-value tests and edge-case tests, all executed by Cargo without external numerical programs.

## 6 Phase 3 Implement special functions

Reference: [GSL special functions](https://www.gnu.org/software/gsl/doc/html/specfunc.html).

Implement each family in `gsl-oxide-special` from public mathematical definitions and documented numerical methods, using safe Rust code and tests with accurate provenance. Link the relevant reference in its Rust documentation. Preserve any applicable third-party permissions and notices when reusing material; separate algorithm documents are unnecessary.

- [ ] **P3.1** Implement `erf` and `erfc`, including signs, zero, infinities, small arguments, and tails. Preserve the accuracy of the small complementary tail by following the appropriate existing method.
- [ ] **P3.2** Implement log gamma and gamma. Cover positive arguments, supported negative noninteger arguments, the sign of gamma, poles, and overflow. Return log absolute gamma and its sign when that operation is exposed.
- [ ] **P3.3** Implement beta and log beta using the established stable formulas. Test large parameters and cases where intermediate gamma values would overflow.
- [ ] **P3.4** Implement regularized incomplete gamma `P(a, x)` and `Q(a, x)`. Follow the existing series/continued-fraction and region handling as needed, including endpoints, complementary tails, and bounded iteration.
- [ ] **P3.5** Implement regularized incomplete beta and its complement. Test shape parameters, values near zero and one, symmetry, and convergence failures.
- [ ] **P3.6** Implement Bessel J, Y, I, and K for order zero and one first. Include the established scaled I and K forms where necessary to avoid overflow or underflow.
- [ ] **P3.7** Extend Bessel support to integer orders using the existing stable recurrence strategy. State the supported order and argument ranges in the function documentation and test them. Keep unsupported real-order or wider-domain variants out of the public API until implemented and tested.
- [ ] **P3.8** Add focused expected-value and boundary tests with each function. Reuse private arithmetic helpers only where more than one implementation benefits.
- [ ] **P3.9** Update the README status as each family becomes available. Add a short working Rust example once the first family is usable.

**Done when:** each exposed function passes its GSL-reference expected-value tests, edge cases, and applicable identities, with supported domains documented in its API comments.

## 7 Phase 4 Implement numerical integration

Reference: [GSL numerical integration](https://www.gnu.org/software/gsl/doc/html/integration.html).

- [ ] **P4.1** Implement one established Gauss-Kronrod rule pair and its error estimate in `gsl-oxide-integration`. Use coefficients with a documented, compatible origin or compute them independently; test known polynomial integrals, symmetry, interval mapping, and zero-length intervals.
- [ ] **P4.2** Implement adaptive finite-interval integration using the established subdivision and error-estimation approach. Start with ordinary `Vec` storage and a simple interval-management strategy; optimize only when necessary.
- [ ] **P4.3** Accept a Rust closure, endpoints, `f64` absolute/relative tolerances, and a `usize` work limit. Return the integral and estimated absolute error using a tuple inside `Result` when failure reporting is needed.
- [ ] **P4.4** Validate tolerance values and work limits. Handle reversed bounds, exhausted iterations, roundoff stagnation, and intervals that cannot be subdivided further. Bound callback evaluations.
- [ ] **P4.5** Support caller-supplied breakpoints with a slice. Check ordering and bounds, and use a consistent total error budget across subintervals.
- [ ] **P4.6** Add semi-infinite and doubly infinite integration using the existing transformations and their Jacobians. Test Gaussian and heavy-tail integrals with known answers.
- [ ] **P4.7** Add an existing method for integrable endpoint singularities. Keep the supported problem classes explicit; more specialized oscillatory and principal-value methods come later.
- [ ] **P4.8** Test smooth functions, cancellation, discontinuities, narrow peaks, endpoint singularities, infinite intervals, and non-finite callbacks. Check the returned estimate against known answers and test failure behavior separately.
- [ ] **P4.9** Add reusable scratch storage only when repeated-call usage justifies it. Prefer a caller-provided standard collection when practical; use a small workspace struct only if it makes required state management clearer. Test reset and reuse behavior if added.

**Done when:** supported integration methods pass known-answer and failure tests using only Rust closures and Cargo tests.

## 8 Phase 5 Implement selected statistics

Reference: [GSL statistics](https://www.gnu.org/software/gsl/doc/html/statistics.html).

Implement these operations in `gsl-oxide-statistics` using original Rust code and the existing mathematical definitions.

- [ ] **P5.1** Implement mean, population/sample variance, and standard deviation on `&[f64]`. Use established stable accumulation rather than raw-moment subtraction where it loses accuracy.
- [ ] **P5.2** Implement covariance and correlation on paired slices. Validate lengths, sample counts, and zero-variance cases.
- [ ] **P5.3** Implement weighted mean and variance using separate, clearly named behavior for the supported weight definitions. Test zero, negative, non-finite, and all-zero weights.
- [ ] **P5.4** Define and test empty/singleton data, constant data, NaN/infinity, degrees of freedom, and overflow behavior. Prefer straightforward function names over custom estimator or policy types.
- [ ] **P5.5** Add streaming accumulation and merging after the slice functions work. Reuse existing stable update formulas. Introduce a small state struct only for the accumulated count and moments that must persist between calls.
- [ ] **P5.6** Test ordinary samples, alternating scales, extreme magnitudes, and a large offset with small representable variation. Compare streaming and merged results with the batch implementation using suitable tolerances.
- [ ] **P5.7** Keep slice and streaming operations allocation-free when their computation permits it. Add quantiles, histograms, or higher moments only after a concrete need and an existing GSL definition are identified.

**Done when:** exposed statistics return the expected values and handle degenerate inputs consistently, with tests for numerical stability and streaming behavior.

## 9 Phase 6 Add CI and release the implemented subset

- [ ] **P6.1** Add CI for `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace`, and `cargo test --workspace --doc`. Run numerical tests in debug and release mode, and check that each package builds separately.
- [ ] **P6.2** Run ordinary test jobs on Windows, Linux, and macOS where runners are available. Add ARM64 execution when supported rather than building a separate validation system.
- [ ] **P6.3** Confirm that builds and tests use no GSL bindings, C/C++ harness, external expected-value generator, or native GSL installation. Keep reference source outside the published package.
- [ ] **P6.4** Review the public API for unnecessary custom types, duplicated entry points, redundant configuration, and avoidable dependencies. Keep only the minimal useful interface.
- [ ] **P6.5** Update the README with implemented modules, supported families, known restrictions, links to alternatives, and examples that compile as Rust doctests where practical.
- [ ] **P6.6** Inspect `cargo package -p <package> --list` and build each package before release. Ensure every archive includes its applicable license, notices, and README; exclude reference source and unrelated files. Use path and version requirements together if a real dependency between published workspace crates is introduced.
- [ ] **P6.7** Publish the tested crates separately so users can install only the required domain. Start with a consistent version policy to keep maintenance simple; publish packages as their tested implementations become usable. Add expected-value and regression tests with every new function or numerical fix.

**Done when:** the documented subset is installable, its examples work, and its Cargo tests pass on the supported targets.

## 10 Later extensions

After the initial subset works, add exponential integrals and selected orthogonal polynomials, then wider Bessel domains or specialized integration methods when needed. Follow the same short GSL-reference-plus-Rust-tests workflow.

Arrays, linear algebra, FFT, optimization, general RNG infrastructure, and existing distribution ecosystems remain outside this workspace's implementation scope. Maintain the alternative-library links in the README instead of adding adapters or compatibility crates.

Performance work is optional after correctness tests pass: profile the Rust implementation, reduce unnecessary allocations, and add a small Rust benchmark only for an observed bottleneck. A cross-language benchmark suite is unnecessary.

## 11 Implementation order and definition of done

Start with the README and three workspace member skeletons, then implement `erf` and `erfc` in `gsl-oxide-special` with ordinary Rust tests. Continue with log gamma, gamma, beta, incomplete gamma/beta, and Bessel. Integration and statistics use their own package and test structures and do not require the special-functions crate unless a concrete implementation needs it.

A function is complete when its public API uses standard types where practical, its GSL reference and important input behavior are documented briefly, expected-value and edge-case tests pass, and it runs entirely in Rust. Add dedicated custom types, tools, documents, or optimization only to solve a demonstrated need.

# GSL function name equivalencies

This guide helps users moving from GSL find the corresponding gsl-oxide functions. It also records Rust function renames made to clarify the calculation, normalization, or return value. Names that are already clear, including established notation such as `erf` and `bessel::j0`, are retained.

The tables cover the public functions currently implemented in this workspace. They describe corresponding operations within the Rust API's documented domains, not identical signatures, algorithms, accuracy, or error handling. Check the full API documentation for those details.

## Statistics

Rust names in this table are relative to `gsl_oxide_statistics`.

| GSL function | Rust function | Calculation |
| --- | --- | --- |
| `gsl_stats_mean` | `mean` | Arithmetic mean. |
| `gsl_stats_variance` | `sample_variance` | Variance with denominator `n - 1`. |
| `gsl_stats_sd` | `sample_standard_deviation` | Square root of sample variance. |
| `gsl_stats_covariance` | `sample_covariance` | Covariance with denominator `n - 1`. |
| `gsl_stats_correlation` | `pearson_correlation` | Pearson correlation of paired observations. |
| `gsl_stats_spearman` | `spearman_correlation` | Pearson correlation of ranks, using average ranks for ties. |
| `gsl_stats_wmean` | `weighted_mean` | Reliability-weighted mean. |
| `gsl_stats_wvariance` | `weighted_reliability_variance` | Weighted variance with denominator `W - sum(w*w)/W`. |

Rust accepts contiguous `f64` slices without stride or count arguments. Weighted functions take `(data, weights)`, whereas GSL takes weights first. Spearman allocates its ranking storage internally. Non-finite observations and invalid weights follow the Rust crate's checked error contract.

### Running statistics

Rust methods belong to `gsl_oxide_statistics::RunningStatistics`.

| GSL operation | Rust constructor or method |
| --- | --- |
| `gsl_rstat_alloc` | `RunningStatistics::new` or `RunningStatistics::default` |
| `gsl_rstat_add` | `RunningStatistics::push` |
| `gsl_rstat_n` | `RunningStatistics::count` |
| `gsl_rstat_reset` | `RunningStatistics::reset` |
| `gsl_rstat_mean` | `RunningStatistics::mean` |
| `gsl_rstat_variance` | `RunningStatistics::sample_variance` |
| `gsl_rstat_sd` | `RunningStatistics::sample_standard_deviation` |
| `gsl_rstat_free` | Automatic cleanup when the Rust value is dropped. |

### Additional Rust operations

These functions have no direct GSL entry point with the same inputs and meaning.

| Rust function or method | Calculation |
| --- | --- |
| `population_variance` | Variance about the computed mean, with denominator `n`. |
| `population_standard_deviation` | Square root of population variance. |
| `weighted_population_variance` | Variance about the computed weighted mean, with denominator `W`. |
| `RunningStatistics::population_variance` | Running variance with denominator `n`. |
| `RunningStatistics::population_standard_deviation` | Square root of running population variance. |
| `RunningStatistics::merge` | Combine two accumulators. |

GSL's `*_variance_with_fixed_mean` functions use a supplied population mean; they are not direct replacements for these population-variance functions. `gsl_stats_pvariance` computes pooled variance, not population variance.

## Special functions

Rust paths in this section are relative to `gsl_oxide_special`. Where GSL also supplies an `_e` variant, the row identifies the mathematical function only: Rust's `Result` reports failures but does not provide GSL's numerical error estimate.

| GSL function | Rust function |
| --- | --- |
| `gsl_sf_erf` | `error_functions::erf` |
| `gsl_sf_erfc` | `error_functions::erfc` |
| `gsl_sf_gamma` | `gamma::gamma` |
| `gsl_sf_lngamma_sgn_e` | `gamma::ln_abs_gamma_with_sign` |
| `gsl_sf_beta` | `beta::beta` |
| `gsl_sf_lnbeta` | `beta::ln_beta` |
| `gsl_sf_gamma_inc_P` | `incomplete_gamma::lower_regularized` |
| `gsl_sf_gamma_inc_Q` | `incomplete_gamma::upper_regularized` |
| `gsl_sf_beta_inc` | `incomplete_beta::regularized` |
| `gsl_sf_legendre_Pl` | `legendre::polynomial` |
| `gsl_sf_laguerre_n` | `laguerre::generalized_polynomial` |

`gamma::ln_abs_gamma_with_sign(x)` returns `(ln(abs(Gamma(x))), sign)`, where the sign is an `i32` equal to `+1` or `-1`. Use the first tuple component when transitioning from `gsl_sf_lngamma`.

`incomplete_beta::regularized_complement(a, b, x)` evaluates `1 - I_x(a, b)` directly to retain small upper tails. GSL has no separate special-function entry point for this complement; mathematically it equals `gsl_sf_beta_inc(b, a, 1 - x)`, but forming `1 - x` can lose precision.

### Cylindrical Bessel functions

These correspond to GSL's uppercase `J`, `Y`, `I`, and `K` cylindrical families. GSL's lowercase Bessel names denote spherical functions, which this workspace does not implement.

| GSL function | Rust function |
| --- | --- |
| `gsl_sf_bessel_J0` | `bessel::j0` |
| `gsl_sf_bessel_J1` | `bessel::j1` |
| `gsl_sf_bessel_Jn` | `bessel::jn` |
| `gsl_sf_bessel_Y0` | `bessel::y0` |
| `gsl_sf_bessel_Y1` | `bessel::y1` |
| `gsl_sf_bessel_Yn` | `bessel::yn` |
| `gsl_sf_bessel_I0` | `bessel::i0` |
| `gsl_sf_bessel_I1` | `bessel::i1` |
| `gsl_sf_bessel_In` | `bessel::i_n` |
| `gsl_sf_bessel_I0_scaled` | `bessel::i0_scaled` |
| `gsl_sf_bessel_I1_scaled` | `bessel::i1_scaled` |
| `gsl_sf_bessel_In_scaled` | `bessel::i_n_scaled` |
| `gsl_sf_bessel_K0` | `bessel::k0` |
| `gsl_sf_bessel_K1` | `bessel::k1` |
| `gsl_sf_bessel_Kn` | `bessel::kn` |
| `gsl_sf_bessel_K0_scaled` | `bessel::k0_scaled` |
| `gsl_sf_bessel_K1_scaled` | `bessel::k1_scaled` |
| `gsl_sf_bessel_Kn_scaled` | `bessel::kn_scaled` |

Scaled `I` functions return `exp(-abs(x)) * I_n(x)`; scaled `K` functions return `exp(x) * K_n(x)`.

### Exponential integrals

| GSL function | Rust function |
| --- | --- |
| `gsl_sf_expint_E1` | `exponential_integrals::e1` |
| `gsl_sf_expint_E2` | `exponential_integrals::e2` |
| `gsl_sf_expint_Ei` | `exponential_integrals::ei` |
| `gsl_sf_expint_E1_scaled` | `exponential_integrals::e1_scaled` |
| `gsl_sf_expint_E2_scaled` | `exponential_integrals::e2_scaled` |
| `gsl_sf_expint_Ei_scaled` | `exponential_integrals::ei_scaled` |

Scaled `E1` and `E2` return `exp(x)` times the unscaled integral; scaled `Ei` returns `exp(-x) * Ei(x)`.

## Numerical integration

Rust names are relative to `gsl_oxide_integration`. These are counterparts by integration task; they do not reproduce GSL's algorithms or all supported cases.

| GSL function | Rust function | Integration task |
| --- | --- | --- |
| `gsl_integration_qag` | `integrate` | Adaptive integration over finite bounds. |
| `gsl_integration_qagp` | `integrate_with_breakpoints` | Finite interval with supplied breakpoints. |
| `gsl_integration_qagiu` | `integrate_upper_infinite` | From a finite lower bound to positive infinity. |
| `gsl_integration_qagil` | `integrate_lower_infinite` | From negative infinity to a finite upper bound. |
| `gsl_integration_qagi` | `integrate_infinite` | Over the whole real line. |
| `gsl_integration_qags` | `integrate_singular` | Integrable endpoint singularities within the documented Rust limits. |

Rust takes a closure, manages its own storage, and returns `(integral, estimated_absolute_error)` inside `Result`. Its `max_evals` parameter limits callback evaluations; GSL's workspace/subinterval limit is not interchangeable with it. Finite integration uses a fixed Gauss 7 / Kronrod 15 rule, and endpoint-singularity integration uses tanh-sinh. The latter supports endpoint logarithms and powers with exponent at least `-1/2`; consult its API contract before replacing a GSL `qags` call.

## GSL references

The GSL 2.8 definitions are documented in [Statistics](https://www.gnu.org/software/gsl/doc/html/statistics.html), [Running Statistics](https://www.gnu.org/software/gsl/doc/html/rstat.html), [Special Functions](https://www.gnu.org/software/gsl/doc/html/specfunc.html), and [Numerical Integration](https://www.gnu.org/software/gsl/doc/html/integration.html).

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

## Integer powers

Import `gsl_oxide_special::Pow as _` to use these methods on `f64`. The trait uses GSL 2.8's multiplication order. Signed exponents cover the full `i32` range; unsigned exponents cover the full `u32` range.

| GSL function | Rust method |
| --- | --- |
| `gsl_pow_int(x, n)` | `x.pow_int(n)` |
| `gsl_pow_uint(x, n)` | `x.pow_uint(n)` |
| `gsl_pow_2(x)` through `gsl_pow_9(x)` | `x.pow_2()` through `x.pow_9()` |

These methods return `f64` directly, with no overflow or underflow checks, as in GSL. A zero exponent returns `1.0` even for zero, infinite, or NaN bases. Negative exponents take the reciprocal of the base before repeated squaring.

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

## Random distributions and sampling

Rust names in this section are relative to `gsl_oxide_randist`. The tables cover all 102 public function declarations in GSL 2.8's `gsl_randist.h`, including compatibility aliases and memory management.

Sampling takes a caller-owned `&mut RngCore` from `rand_core` 0.6 (also used by `rand` 0.8); the trait is re-exported by the crate. The crate supplies distributions, not RNG implementations or seeding. Scalar samplers generally return `Result`, while standard-normal, Landau, and fixed-dimensional direction samplers return values directly. Densities and probability masses use checked `Result` APIs. See the crate documentation for domains and range errors.

The `_density` suffix identifies a continuous density; `_pdf` identifies a discrete probability; `_log_density` and `_log_pdf` use natural logarithms. At singular support boundaries Rust uses limiting densities, which can be infinite where GSL sometimes returns zero. These are distribution counterparts; RNG sequences and implementation algorithms can differ. In particular, `binomial` uses Knuth splitting and geometric waiting times rather than GSL's BTPE implementation. GSL's `binomial_tpe` and `gamma_mt` compatibility aliases map to the corresponding default Rust samplers.

| GSL function | Rust function or method |
| --- | --- |
| `gsl_ran_bernoulli` | `bernoulli` |
| `gsl_ran_bernoulli_pdf` | `bernoulli_pdf` |
| `gsl_ran_beta` | `beta` |
| `gsl_ran_beta_pdf` | `beta_density` |
| `gsl_ran_binomial` | `binomial` |
| `gsl_ran_binomial_knuth` | `binomial_knuth` |
| `gsl_ran_binomial_tpe` | `binomial` |
| `gsl_ran_binomial_pdf` | `binomial_pdf` |
| `gsl_ran_exponential` | `exponential` |
| `gsl_ran_exponential_pdf` | `exponential_density` |
| `gsl_ran_exppow` | `exponential_power` |
| `gsl_ran_exppow_pdf` | `exponential_power_density` |
| `gsl_ran_cauchy` | `cauchy` |
| `gsl_ran_cauchy_pdf` | `cauchy_density` |
| `gsl_ran_chisq` | `chi_squared` |
| `gsl_ran_chisq_pdf` | `chi_squared_density` |
| `gsl_ran_dirichlet` | `dirichlet` |
| `gsl_ran_dirichlet_pdf` | `dirichlet_density` |
| `gsl_ran_dirichlet_lnpdf` | `dirichlet_log_density` |
| `gsl_ran_erlang` | `erlang` |
| `gsl_ran_erlang_pdf` | `erlang_density` |
| `gsl_ran_fdist` | `fisher_f` |
| `gsl_ran_fdist_pdf` | `fisher_f_density` |
| `gsl_ran_flat` | `uniform_interval` |
| `gsl_ran_flat_pdf` | `uniform_interval_density` |
| `gsl_ran_gamma` | `gamma` |
| `gsl_ran_gamma_int` | `gamma_integer` |
| `gsl_ran_gamma_pdf` | `gamma_density` |
| `gsl_ran_gamma_mt` | `gamma` |
| `gsl_ran_gamma_knuth` | `gamma_knuth` |
| `gsl_ran_gaussian` | `normal` |
| `gsl_ran_gaussian_ratio_method` | `normal_ratio` |
| `gsl_ran_gaussian_ziggurat` | `normal_ziggurat` |
| `gsl_ran_gaussian_pdf` | `normal_density` |
| `gsl_ran_ugaussian` | `standard_normal` |
| `gsl_ran_ugaussian_ratio_method` | `standard_normal_ratio` |
| `gsl_ran_ugaussian_pdf` | `standard_normal_density` |
| `gsl_ran_gaussian_tail` | `normal_tail` |
| `gsl_ran_gaussian_tail_pdf` | `normal_tail_density` |
| `gsl_ran_ugaussian_tail` | `standard_normal_tail` |
| `gsl_ran_ugaussian_tail_pdf` | `standard_normal_tail_density` |
| `gsl_ran_bivariate_gaussian` | `bivariate_normal` |
| `gsl_ran_bivariate_gaussian_pdf` | `bivariate_normal_density` |
| `gsl_ran_multivariate_gaussian` | `multivariate_normal` |
| `gsl_ran_multivariate_gaussian_log_pdf` | `multivariate_normal_log_density` |
| `gsl_ran_multivariate_gaussian_pdf` | `multivariate_normal_density` |
| `gsl_ran_multivariate_gaussian_mean` | `multivariate_mean` |
| `gsl_ran_multivariate_gaussian_vcov` | `multivariate_sample_covariance` |
| `gsl_ran_wishart` | `wishart` |
| `gsl_ran_wishart_log_pdf` | `wishart_log_density` |
| `gsl_ran_wishart_pdf` | `wishart_density` |
| `gsl_ran_landau` | `landau` |
| `gsl_ran_landau_pdf` | `landau_density` |
| `gsl_ran_geometric` | `geometric` |
| `gsl_ran_geometric_pdf` | `geometric_pdf` |
| `gsl_ran_hypergeometric` | `hypergeometric` |
| `gsl_ran_hypergeometric_pdf` | `hypergeometric_pdf` |
| `gsl_ran_gumbel1` | `gumbel_type_1` |
| `gsl_ran_gumbel1_pdf` | `gumbel_type_1_density` |
| `gsl_ran_gumbel2` | `gumbel_type_2` |
| `gsl_ran_gumbel2_pdf` | `gumbel_type_2_density` |
| `gsl_ran_logistic` | `logistic` |
| `gsl_ran_logistic_pdf` | `logistic_density` |
| `gsl_ran_lognormal` | `lognormal` |
| `gsl_ran_lognormal_pdf` | `lognormal_density` |
| `gsl_ran_logarithmic` | `logarithmic_series` |
| `gsl_ran_logarithmic_pdf` | `logarithmic_series_pdf` |
| `gsl_ran_multinomial` | `multinomial` |
| `gsl_ran_multinomial_pdf` | `multinomial_pdf` |
| `gsl_ran_multinomial_lnpdf` | `multinomial_log_pdf` |
| `gsl_ran_negative_binomial` | `negative_binomial` |
| `gsl_ran_negative_binomial_pdf` | `negative_binomial_pdf` |
| `gsl_ran_pascal` | `pascal` |
| `gsl_ran_pascal_pdf` | `pascal_pdf` |
| `gsl_ran_pareto` | `pareto` |
| `gsl_ran_pareto_pdf` | `pareto_density` |
| `gsl_ran_poisson` | `poisson` |
| `gsl_ran_poisson_array` | `fill_poisson` |
| `gsl_ran_poisson_pdf` | `poisson_pdf` |
| `gsl_ran_rayleigh` | `rayleigh` |
| `gsl_ran_rayleigh_pdf` | `rayleigh_density` |
| `gsl_ran_rayleigh_tail` | `rayleigh_tail` |
| `gsl_ran_rayleigh_tail_pdf` | `rayleigh_tail_density` |
| `gsl_ran_tdist` | `student_t` |
| `gsl_ran_tdist_pdf` | `student_t_density` |
| `gsl_ran_laplace` | `laplace` |
| `gsl_ran_laplace_pdf` | `laplace_density` |
| `gsl_ran_levy` | `stable_symmetric` |
| `gsl_ran_levy_skew` | `stable_skewed` |
| `gsl_ran_weibull` | `weibull` |
| `gsl_ran_weibull_pdf` | `weibull_density` |
| `gsl_ran_dir_2d` | `direction_2d` |
| `gsl_ran_dir_2d_trig_method` | `direction_2d_trigonometric` |
| `gsl_ran_dir_3d` | `direction_3d` |
| `gsl_ran_dir_nd` | `direction_nd` |
| `gsl_ran_shuffle` | `shuffle` |
| `gsl_ran_choose` | `choose` |
| `gsl_ran_sample` | `sample_with_replacement` |
| `gsl_ran_discrete_preproc` | `DiscreteDistribution::new` |
| `gsl_ran_discrete_free` | Automatic cleanup when DiscreteDistribution is dropped. |
| `gsl_ran_discrete` | `DiscreteDistribution::sample` |
| `gsl_ran_discrete_pdf` | `DiscreteDistribution::pdf` |

Slices replace pointer/count pairs; returned vectors and tuples replace output pointers, and scratch storage is managed internally. `multinomial(rng, trials, weights)` accepts unnormalized nonnegative weights. `dirichlet` takes concentration parameters; its density takes concentrations followed by simplex coordinates. `geometric` counts trials through the first success; negative binomial and Pascal count failures before the requested successes. Hypergeometric draws are clamped to the population as in GSL.

Multivariate matrices are flattened row-major slices. Normal sampling and density take a mean vector and the lower Cholesky factor of its covariance; the dimension comes from the mean. Upper-triangular factor entries are ignored. The mean and sample-covariance estimators take observations followed by their column count; covariance divides by `rows - 1`. `wishart(rng, degrees, scale_cholesky, dimensions)` returns a square matrix. Wishart density functions take `(sample_cholesky, degrees, scale_cholesky, dimensions)`: the observation factor supplies its determinant and the trace term, so GSL's redundant observation matrix and work matrix are unnecessary.

`choose` preserves source order and samples without replacement; `sample_with_replacement` permits repeated items. Both return cloned items. `shuffle` mutates a borrowed slice. `DiscreteDistribution::new` preprocesses weights into a Walker alias table, with automatic cleanup. Its `pdf(category)` method returns the normalized probability for a zero-based category index, or zero outside the table.

## GSL references

The GSL 2.8 definitions are documented in [Statistics](https://www.gnu.org/software/gsl/doc/html/statistics.html), [Running Statistics](https://www.gnu.org/software/gsl/doc/html/rstat.html), [Special Functions](https://www.gnu.org/software/gsl/doc/html/specfunc.html), [Numerical Integration](https://www.gnu.org/software/gsl/doc/html/integration.html), and [Random Number Distributions](https://www.gnu.org/software/gsl/doc/html/randist.html).

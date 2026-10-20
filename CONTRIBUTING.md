# Contributing

Contributions are welcome! Please follow the guidelines below.

## References and testing

The public [GSL 2.8 definitions and documented behavior](https://www.gnu.org/software/gsl/doc/html/) are the behavior reference. Implementations and tests are independently written from mathematical definitions and established numerical methods, without copying or translating GSL implementation code or its coefficient tables. Short mathematical provenance notes live beside the Rust source.

Methods use DLMF formulas: the error-function Laplace fraction, gamma recurrence/reflection and Bernoulli Stirling terms, combined beta log-ratios, incomplete-gamma series/fractions, and incomplete-beta symmetry/series/fractions. Bessel base orders use small-argument series, Taylor propagation of the ordinary Bessel equation, a private K integral kernel, and asymptotic expansions with generated product coefficients. Integer-order recurrences use stable directions, rescaling, and normalization. There are no fitted coefficient tables or added dependencies.

Unit, integration, and documentation tests run through Cargo with expected values stored directly in Rust, justified floating-point tolerances, and explicit edge cases. Tests require no external scripts, numerical fixtures, GSL, or C/C++ numerical compiler.

Reference cases specify both absolute and relative tolerances. Finite comparisons validate the tolerances and use `abs(actual - expected) <= max(abs_tol, rel_tol * abs(expected))`; NaN, infinity, and signed-zero behavior use direct assertions, including zero sign bits. Deterministic identity checks supplement the reference values, with explicit coverage of subnormals and adjacent representable branch-boundary inputs.

Numerical anchors come from published mpmath/SciPy documentation examples and exact mathematical results or explicitly bounded formula truncations, cited beside each test. Neither library is executed or required. Most full-precision anchors allow roughly 1e-14–1e-13 relative roundoff; rounded published tables and asymptotic checks use bounds appropriate to their precision. Tiny-tail comparisons exclude zero. Private kernels exercise real iteration exhaustion with reduced budgets. Correct rounding, GSL execution comparisons, and cross-platform bitwise agreement are not claimed.

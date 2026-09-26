# Changelog

## ppgt 0.1.0 (2026-04-22)

### Initial stable release

This version represents the first stable, internally verified release of
**ppgt** following a full end-to-end validation workflow (Sessions A–H).

### Core functionality

- Implemented and verified:
  - Polynomial Pooling (PP)
  - PBEST design
  - HYPER design (Hong et al. 2022)
  - Separable and disjunct matrix constructions
  - Finite-field GF(q) arithmetic utilities

### Major corrections

- **HYPER-EC re-scoped**
  - [`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md)
    and
    [`hyper_ec_decode()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_decode.md)
    are now documented as **experimental parity-augmented extensions**
  - Removed incorrect claims of deterministic error-correction
    guarantees
  - Documentation now clearly distinguishes from Hong et al. (2022)
- **PP decoding improved**
  - Exhaustive COMP fallback added before GPSR
  - Eliminates false negatives in tested scenarios
- **Design validation fixed**
  - [`pp_design()`](https://myaseen208.github.io/ppgt/reference/pp_design.md)
    now enforces valid compression constraints
  - Correct filtering of infeasible configurations

### Finite-field implementation

- Verified correctness of:
  - [`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
    [`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
    [`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
    [`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md)
- Supports:
  - Prime fields GF(p)
  - Binary extensions GF(2^n)
  - Selected GF(p^n) extensions

### Documentation

- Full roxygen2 coverage for all exported functions
- Mathematical notation aligned with `Notations.tex`
- Updated examples to reflect correct design usage
- Vignette suite split into focused topic articles for PP, P-BEST,
  HYPER, HyperEC, GF(q), separable/disjunct designs, comparisons, and
  diagnostics
- Release-facing vignette formatting, printed example outputs, and
  pkgdown site metadata refined for the 0.1.0 documentation build

### Testing

- Full testthat Edition 3 suite
- High coverage across all core modules
- 100% coverage achieved for critical components:
  - PP designs
  - decoding
  - separable/disjunct constructions

### Package quality

- [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  clean:
  - 0 errors
  - 0 warnings
  - 0 notes

------------------------------------------------------------------------

## ppgt 0.0.0

- Development versions (pre-validation)

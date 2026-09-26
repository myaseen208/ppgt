

# ppgt 0.1.2 (2026-09-26)

## New vignette

- Added `ppgt_pp_188x282`, “Building a 188 x 282 pooling matrix”: a
  worked Polynomial Pools example
  (`pp_matrix(q = 47, d = 3, nl = 4, N = 282)`) covering construction,
  verification (dimensions, row/column sums, the pairwise-overlap
  bound), efficiency, and a `pp_decode()` example.

## Kirkman triple system (KTS) helpers: new deterministic constructions

- `kirkman_matrix()` now uses the affine geometry `AG(n,3)` construction
  for every `v` that is an exact power of 3 (`v = 3, 9, 27, 81, ...`),
  not just `v = 3` and `v = 9`. This is a direct, deterministic
  construction (no randomized search), so `v = 27` and `v = 81` now join
  `v = 3`, `9`, and `15` as fast, guaranteed constructions; `v = 9` also
  moved from the general randomized search onto this same deterministic
  path.
- `v = 21` and `v = 33` were investigated for a dedicated deterministic
  construction (in the same spirit as `v = 15`’s PG(3,2) spread packing,
  and a `PG(5,2)`-spread packing was attempted for `v = 63`). None was
  completed within a reasonable computational budget in this round of
  work, and empirically the existing general randomized search is not
  reliable for `v = 21` either (it did not find a system within its
  budget in testing). This is now documented explicitly in
  `?kirkman_matrix`: guaranteed-fast `v` are powers of 3 and `v = 15`;
  all other admissible `v` use the general search, which is not
  guaranteed to succeed.
- Added tests for every fast-path `v` (`3, 9, 15, 27, 81`) and for the
  internal dispatch helper; the general search itself is now also
  covered directly (via a small, fast `v`) rather than only through the
  slow, skip-gated end-to-end tests for `v = 21`/`v = 33`.

## Bug fix: `compare_all_designs_honest()` reported hard-coded, sometimes-wrong efficiency comparisons

`compare_all_designs_honest()`’s printed “PACKAGE-SPECIFIC SUMMARY”
lines and its `Rank_Efficiency` columns were hard-coded rather than
computed from the function’s own tables, and had drifted out of sync
with them – e.g. it reported HYPER and Tapestry as “tied” at `N = 256`
when their actual computed `m/N` differ (`0.172` vs `0.176`), and its
`N = 384` ranking placed Tapestry above HYPER-EC even though HYPER-EC’s
computed ratio is lower (more efficient). `Rank_Efficiency` is now
computed from the realized `m_over_N` column, and the printed narrative
and returned `$summary` are derived from the same computation. Also
switched this function from `cat()`/`print()` to `message()`, matching
the rest of the package.

## Other fixes

- The `URL` field’s `https://myaseen208.github.io/ppgt` link (and the
  same link in `README`) redirected (301) to a custom domain that itself
  404s, which both triggered an `R CMD check` NOTE and pointed users at
  a dead page; the link has been removed rather than replaced with its
  (broken) redirect target. The GitHub repository URL remains.
- `.Rbuildignore` now also excludes `vignettes/.quarto` and
  `vignettes/.Rhistory`, fixing an `R CMD check` NOTE about hidden
  files.

## Contributors

Added four contributors (`ctb`) to `Authors@R`: Christopher McMahan,
Christopher Bilder, Joshua Tebbs, and Pranta Das. Muhammad Yaseen
remains the sole author, maintainer, and copyright holder.

# ppgt 0.1.1 (2026-09-26)

## Bug fix: Kirkman triple system (KTS) helpers did not satisfy the KTS defining properties

A reviewer reported that `kirkman_matrix()` did not produce genuine
Kirkman triple systems for some `v`. Investigation confirmed and
precisely quantified the problem:

- `kirkman_matrix(15)` used a hand-written block set whose last three
  “parallel classes” were not valid partitions of the 15 points: 19
  pairs of points were covered by more than one block and 24 pairs were
  never covered at all (out of the 105 pairs a KTS(15) must cover
  exactly once). Verified against an independently constructed and
  checked reference KTS(15) (a packing of the 35 lines of PG(3,2) into 7
  spreads).
- `kirkman_matrix(21)` used an internal helper that generated only 7 of
  the 70 blocks a KTS(21) requires – an incomplete construction, not a
  numerical error, that silently returned a badly undersized system.
- `kirkman_matrix(3)` was rejected outright, even though `v = 3` is the
  trivial (one block, one parallel class) Kirkman triple system.
- Invalid `v` (not `== 3 mod 6`) returned `NULL` with a `warning()`
  instead of failing, so a caller that didn’t check for `NULL` could
  silently propagate a missing design.

### What changed

- `kirkman_matrix()` was rewritten around a general, uniform
  construction rather than per-`v` lookup tables: a dedicated, verified
  construction for `v = 15` (packing the lines of PG(3,2) into spreads –
  the classical solution to Kirkman’s schoolgirl problem), and a general
  randomized candidate-pool exact-cover search for every other
  admissible `v` (including `v = 3` and `v = 9`, which construct in well
  under a second). Every result is independently re-verified before
  being returned; the function errors rather than ever returning an
  unverified or invalid system.
- Invalid `v` (not an integer `>= 3`, or not `== 3 mod 6`) now raises an
  informative `stop()` instead of a `warning()` + `NULL`.
- For `v` where the general search does not find a system within its
  computational budget (currently only guaranteed fast for `v` in
  `{3, 9, 15}`), `kirkman_matrix()` raises an informative error
  explaining this is a search-budget limitation, not evidence that no
  KTS exists for that `v`.
- New exported validator `is_kts()` independently checks all four
  defining KTS properties (block size 3, correct block count, every pair
  of points covered exactly once, and resolvability into parallel
  classes that each partition the points) directly against a candidate
  incidence matrix. Used throughout the package’s own tests, and
  available for users to verify any KTS-shaped matrix themselves.
- `kirkman_matrix()`’s returned `design` attribute now includes
  `parallel_class`, the resolution-class label of each block, so
  resolvability can be checked (by `is_kts()` or by users) directly from
  the returned object.

### Who was affected

Only `kirkman_matrix()` itself, and only for `v = 3`, `v = 21`, and
invalid `v`; `kirkman_matrix(9)` was already correct and is unaffected.
No other exported function calls `kirkman_matrix()` with `v` other than
`9` (`PoolMatrix(family = "kirkman")` and its tests all use `v = 9`), so
no other function, vignette, or previously published result was affected
by this bug.

## Other fixes: RNG state leaked to the caller

Four functions called `set.seed()` internally without saving and
restoring the caller’s prior random-number generator state, so calling
them could silently perturb any subsequent random draws elsewhere in a
user’s session or script:

- `SeparableMatrix()` and `DisjunctMatrix()` (`method = "random"` path)
- `simulate_group_testing()` (when `seed` is supplied)
- `pbest_clinical_matrix()`

All four now save the caller’s `.Random.seed` before seeding and restore
it on exit, matching the pattern already used correctly by
`regular_pooling_matrix()`. Output for a given `seed` is unchanged; only
the side effect on the caller’s RNG stream is fixed.

## Documentation

- `README.qmd` and `NEWS.qmd` were filled in for the first time
  (previously placeholder files).

------------------------------------------------------------------------

# ppgt 0.1.0 (2026-04-22)

## Initial stable release

This version represents the first stable, internally verified release of
**ppgt** following a full end-to-end validation workflow (Sessions A-H).

## Core functionality

- Implemented and verified:
  - Polynomial Pooling (PP)
  - PBEST design
  - HYPER design (Hong et al. 2022)
  - Separable and disjunct matrix constructions
  - Finite-field GF(q) arithmetic utilities

## Major corrections

- **HYPER-EC re-scoped**
  - `hyper_ec_matrix()` and `hyper_ec_decode()` are now documented as
    **experimental parity-augmented extensions**
  - Removed incorrect claims of deterministic error-correction
    guarantees
  - Documentation now clearly distinguishes from Hong et al. (2022)
- **PP decoding improved**
  - Exhaustive COMP fallback added before GPSR
  - Eliminates false negatives in tested scenarios
- **Design validation fixed**
  - `pp_design()` now enforces valid compression constraints
  - Correct filtering of infeasible configurations

## Finite-field implementation

- Verified correctness of:
  - `gf()`, `gf_add()`, `gf_mult()`, `gf_pow()`
- Supports:
  - Prime fields GF(p)
  - Binary extensions GF(2^n)
  - Selected GF(p^n) extensions

## Documentation

- Full roxygen2 coverage for all exported functions
- Mathematical notation aligned with `Notations.tex`
- Updated examples to reflect correct design usage
- Vignette suite split into focused topic articles for PP, P-BEST,
  HYPER, HyperEC, GF(q), separable/disjunct designs, comparisons, and
  diagnostics
- Release-facing vignette formatting, printed example outputs, and
  pkgdown site metadata refined for the 0.1.0 documentation build

## Testing

- Full testthat Edition 3 suite
- High coverage across all core modules
- 100% coverage achieved for critical components:
  - PP designs
  - decoding
  - separable/disjunct constructions

## Package quality

- `devtools::check()` clean:
  - 0 errors
  - 0 warnings
  - 0 notes

------------------------------------------------------------------------

# ppgt 0.0.0

- Development versions (pre-validation)

# ppgt: Theory, Notation, and Package Overview

## 1 Abstract

The `ppgt` package provides pooled group testing matrix constructors,
decoders, finite-field utilities, and comparison helpers. This
introductory article summarizes the package notation, organizes the
design families covered in the package, introduces the incidence-matrix
diagnostics used throughout the vignette suite, and closes with a short
end-to-end example.

## 2 Package Introduction

In the package, a pooling design is represented by a binary incidence
matrix \\\mathbf{M} \in \\0,1\\^{J \times N}\\ with \\J\\ pools and
\\N\\ individuals. The major exported families are:

- Polynomial Pools and PP-based P-BEST matrices
- HYPER-style matrices
- the package-specific experimental HyperEC extension
- finite-field arithmetic utilities for \\GF(q)\\
- separable and disjunct constructions
- package comparison helpers

The vignettes are split by topic so that the package overview remains
concise and the longer mathematical exposition stays with the relevant
design family.

Given latent individual statuses \\\widetilde{\mathbf{y}} =
(\widetilde{y}\_1, \ldots, \widetilde{y}\_N)^\top\\, the package uses
the pooled binary status map \\ \widetilde{z}\_j =
\mathbb{I}\left(\sum\_{i = 1}^N M\_{ji}\widetilde{y}\_i \> 0\right),
\qquad j \in \\1, \ldots, J\\, \\ with latent pool-status vector
\\\widetilde{\mathbf{z}} = (\widetilde{z}\_1, \ldots,
\widetilde{z}\_J)^\top\\ and observed pool outcomes \\\mathbf{z} = (z_1,
\ldots, z_J)^\top\\.

## 3 Package Layers

The package architecture is intentionally layered.

1.  Matrix construction:
    [`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md),
    [`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
    [`SeparableMatrix()`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md),
    [`DisjunctMatrix()`](https://myaseen208.github.io/ppgt/reference/DisjunctMatrix.md),
    the additional design-family constructors, and
    [`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md)
    return a pooling design matrix \\\mathbf{M}\\ or a
    design-construction object.
2.  Decoding: functions such as
    [`pp_comp()`](https://myaseen208.github.io/ppgt/reference/pp_comp.md),
    [`pp_gpsr()`](https://myaseen208.github.io/ppgt/reference/pp_gpsr.md),
    [`pp_decode()`](https://myaseen208.github.io/ppgt/reference/pp_decode.md),
    and
    [`hyper_ec_decode()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_decode.md)
    interpret observed outcomes relative to a fixed design.
3.  Simulation:
    [`simulate_group_testing()`](https://myaseen208.github.io/ppgt/reference/simulate_group_testing.md)
    maps a design and \\\widetilde{\mathbf{y}}\\ to latent and observed
    pool outcomes.
4.  Workflow execution:
    [`run_testing_workflow()`](https://myaseen208.github.io/ppgt/reference/run_testing_workflow.md)
    coordinates the currently implemented workflow skeleton above the
    constructor layer. In the present release, that means either stage-1
    only (`workflow = "non_adaptive"`, `stage2 = "none"`) or stage-1
    followed by individual retesting of the stage-1 non-negative set
    (`workflow = "adaptive_retest"`, `stage2 = "individual_retest"`).
5.  Reporting:
    [`protocol_summary()`](https://myaseen208.github.io/ppgt/reference/protocol_summary.md)
    summarizes a workflow object without constructing a matrix or
    re-running a decoder.

These layers are separate by design. In particular,
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md)
remains a constructor-only dispatcher that returns \\\mathbf{M}\\ and
does not choose adaptive versus non-adaptive workflows, retesting
scenarios, stage-2 protocol logic, or decoder behavior.

## 4 Notation Summary

This package follows the notation defined in `Notations.tex`.

- \\N\\: total number of individuals
- \\i \in \\1, \ldots, N\\\\: individual index
- \\\widetilde{\mathbf{y}} = (\widetilde{y}\_1, \ldots,
  \widetilde{y}\_N)^\top\\: vector of latent infection indicators
- \\J\\: total number of pools
- \\j \in \\1, \ldots, J\\\\: pool index
- \\\mathcal{P}\_j\\: set of individuals in pool \\j\\
- \\n_j = \|\mathcal{P}\_j\|\\: size of pool \\j\\
- \\\mathbf{M} = (M\_{ji}) \in \\0,1\\^{J \times N}\\: pooling design
  matrix
- \\\mathcal{J}\_i = \\j : M\_{ji} = 1\\\\: set of pools containing
  individual \\i\\
- \\w_i = \|\mathcal{J}\_i\|\\: number of pools containing individual
  \\i\\
- \\\widetilde{z}\_j\\: latent status of pool \\j\\
- \\z_j\\: observed test outcome of pool \\j\\
- \\\widetilde{\mathbf{z}} = (\widetilde{z}\_1, \ldots,
  \widetilde{z}\_J)^\top\\: vector of latent pool statuses
- \\\mathbf{z} = (z_1, \ldots, z_J)^\top\\: vector of observed pool
  results

For noiseless binary testing, the latent pool status is
\\\widetilde{z}\_j = \mathbb{I}\left(\sum\_{i = 1}^N M\_{ji}
\widetilde{y}\_i \> 0\right)\\. Collecting the pool-level statuses gives
\\\widetilde{\mathbf{z}} = (\widetilde{z}\_1, \ldots,
\widetilde{z}\_J)^\top\\.

### 4.1 Package Notation Registry

The function
[`ppgt_notations()`](https://myaseen208.github.io/ppgt/reference/ppgt_notations.md)
exposes the package notation registry as a plain data frame. It records
the package-wide canonical notation, family-local notation,
decoder-local notation, display aliases, and deprecated aliases in one
place.

``` r
notations <- ppgt_notations()
notations
                         symbol                     render
  1                     \\alpha                    \\alpha
  2                      \\beta                     \\beta
  3          \\hat{\\mathbf{x}}         \\hat{\\mathbf{x}}
  4                 \\mathbf{A}                \\mathbf{A}
  5               \\mathbf{g}_k              \\mathbf{g}_k
  6                 \\mathbf{x}                \\mathbf{x}
  7                 \\mathbf{y}                \\mathbf{y}
  8                     \\sigma                    \\sigma
  9                       \\tau                      \\tau
  10       F_\\tau(\\mathbf{z})       F_\\tau(\\mathbf{z})
  11                compression                compression
  12                 Lambda_max                 Lambda_max
  13            lines_per_point            lines_per_point
  14                          M                          M
  15                    n_pools                    n_pools
  16                  n_samples                  n_samples
  17            points_per_line            points_per_line
  18                  pool_size                  pool_size
  19           pools_per_sample           pools_per_sample
  20            sample_coverage            sample_coverage
  21                    \\alpha                    \\alpha
  22                     \\beta                     \\beta
  23                    \\gamma                    \\gamma
  24                   \\lambda                   \\lambda
  25                          b                          b
  26                          c                          c
  27                          g                          g
  28                          k                          k
  29                          r                          r
  30                          v                          v
  31                          m                          m
  32                    reorder                    reorder
  33                          d                          d
  34                  k_{\\max}                  k_{\\max}
  35                        n_l                        n_l
  36                          q                          q
  37 \\mathbf{u}_{\\mathcal{S}} \\mathbf{u}_{\\mathcal{S}}
  38               \\mathcal{S}               \\mathcal{S}
  39                \\mathbf{y}                \\mathbf{y}
  40   \\widetilde{\\mathbf{Y}}   \\widetilde{\\mathbf{Y}}
  41                        n_i                        n_i
  42                \\mathbf{M}                \\mathbf{M}
  43        \\mathbf{m}_j^\\top        \\mathbf{m}_j^\\top
  44              \\mathbf{s}_e              \\mathbf{s}_e
  45              \\mathbf{s}_p              \\mathbf{s}_p
  46                \\mathbf{z}                \\mathbf{z}
  47             \\mathcal{J}_i             \\mathcal{J}_i
  48             \\mathcal{P}_j             \\mathcal{P}_j
  49   \\widetilde{\\mathbf{y}}   \\widetilde{\\mathbf{y}}
  50   \\widetilde{\\mathbf{z}}   \\widetilde{\\mathbf{z}}
  51           \\widetilde{y}_i           \\widetilde{y}_i
  52           \\widetilde{z}_j           \\widetilde{z}_j
  53                          i                          i
  54                          j                          j
  55                          J                          J
  56                     M_{ji}                     M_{ji}
  57                          N                          N
  58                        n_j                        n_j
  59                        p_i                        p_i
  60                    S_{e_j}                    S_{e_j}
  61                    S_{p_j}                    S_{p_j}
  62                        w_i                        w_i
  63                        z_j                        z_j
                                                                                meaning
  1                                                            GPSR step-size parameter
  2                                                   GPSR line-search shrinkage factor
  3                                         decoder-local estimate of the sparse signal
  4                                              decoder-local sensing or design matrix
  5                                        decoder-local gradient vector at iteration k
  6                                                 decoder-local unknown sparse signal
  7                                                       decoder-local response vector
  8                                                          GPSR line-search parameter
  9                                                       GPSR regularization parameter
  10                                                            GPSR objective function
  11                                                  display label for the ratio N / J
  12                                        display label for the maximum column weight
  13                                   display label for projective-plane column weight
  14                                       legacy display alias for the number of pools
  15                                              display label for the number of pools
  16                                        display label for the number of individuals
  17                                       display label for projective-plane pool size
  18                                                        display label for pool size
  19                                                    display label for column weight
  20                          display label for column weight in PP verification output
  21                   projective-coordinate parameter or family-local tuning parameter
  22                   projective-coordinate parameter or family-local tuning parameter
  23                                                    projective-coordinate parameter
  24                                                      BIBD pair-incidence parameter
  25                                                         number of blocks in a BIBD
  26                                number of columns in a two-dimensional array design
  27                                                                  Dorfman pool size
  28                                  family-local block size or active-count parameter
  29 family-local row count, replication number, or clinical pools-per-sample parameter
  30                                                           family-local point count
  31                                   HYPER row-count parameter for the stage-1 design
  32                                                                 HYPER reorder flag
  33                                   family-local dimension or target order parameter
  34                                                       PP guaranteed-positive bound
  35                                                                PP number of layers
  36                             family-local field order or pools-per-sample parameter
  37                                     Boolean union pattern induced by an active set
  38                                       active set in separable/disjunct definitions
  39                             forbidden package-wide alias for observed pool results
  40                                  deprecated uppercase-bold latent infection vector
  41                                    deprecated package-wide alias for column weight
  42                                                              pooling design matrix
  43                                                        row j of the pooling matrix
  44                                                 vector of pool-level sensitivities
  45                                                 vector of pool-level specificities
  46                                                    vector of observed pool results
  47                                               set of pools containing individual i
  48                                              set of individuals assigned to pool j
  49                                              vector of latent infection indicators
  50                                                     vector of latent pool statuses
  51                                        latent infection indicator for individual i
  52                                                            latent status of pool j
  53                                                                   individual index
  54                                                                         pool index
  55                                                              total number of pools
  56                                                        entry of the pooling matrix
  57                                                        total number of individuals
  58                                                                     size of pool j
  59                                             infection probability for individual i
  60                                                  pool-level sensitivity for pool j
  61                                                  pool-level specificity for pool j
  62                                                     column weight for individual i
  63                                                         observed outcome of pool j
             category     status                       scope             family
  1     decoder_local      local          decoder-local only            decoder
  2     decoder_local      local          decoder-local only            decoder
  3     decoder_local      local          decoder-local only            decoder
  4     decoder_local      local          decoder-local only            decoder
  5     decoder_local      local          decoder-local only            decoder
  6     decoder_local      local          decoder-local only            decoder
  7     decoder_local      local          decoder-local only            decoder
  8     decoder_local      local          decoder-local only            decoder
  9     decoder_local      local          decoder-local only            decoder
  10    decoder_local      local          decoder-local only            decoder
  11    display_alias      alias summary and display outputs          summaries
  12    display_alias      alias summary and display outputs          summaries
  13    display_alias      alias summary and display outputs          summaries
  14    display_alias      alias summary and display outputs          summaries
  15    display_alias      alias summary and display outputs          summaries
  16    display_alias      alias summary and display outputs          summaries
  17    display_alias      alias summary and display outputs          summaries
  18    display_alias      alias summary and display outputs          summaries
  19    display_alias      alias summary and display outputs          summaries
  20    display_alias      alias summary and display outputs          summaries
  21     family_local      local           family-local only  classical_designs
  22     family_local      local           family-local only  classical_designs
  23     family_local      local           family-local only  classical_designs
  24     family_local      local           family-local only  classical_designs
  25     family_local      local           family-local only  classical_designs
  26     family_local      local           family-local only  classical_designs
  27     family_local      local           family-local only  classical_designs
  28     family_local      local           family-local only  classical_designs
  29     family_local      local           family-local only  classical_designs
  30     family_local      local           family-local only  classical_designs
  31     family_local      local           family-local only              hyper
  32     family_local      local           family-local only              hyper
  33     family_local      local           family-local only           pp/pbest
  34     family_local      local           family-local only           pp/pbest
  35     family_local      local           family-local only           pp/pbest
  36     family_local      local           family-local only           pp/pbest
  37     family_local      local           family-local only separable/disjunct
  38     family_local      local           family-local only separable/disjunct
  39  forbidden_alias deprecated   legacy or forbidden usage             legacy
  40  forbidden_alias deprecated   legacy or forbidden usage             legacy
  41  forbidden_alias deprecated   legacy or forbidden usage             legacy
  42 global_canonical  canonical                package-wide                all
  43 global_canonical  canonical                package-wide                all
  44 global_canonical  canonical                package-wide                all
  45 global_canonical  canonical                package-wide                all
  46 global_canonical  canonical                package-wide                all
  47 global_canonical  canonical                package-wide                all
  48 global_canonical  canonical                package-wide                all
  49 global_canonical  canonical                package-wide                all
  50 global_canonical  canonical                package-wide                all
  51 global_canonical  canonical                package-wide                all
  52 global_canonical  canonical                package-wide                all
  53 global_canonical  canonical                package-wide                all
  54 global_canonical  canonical                package-wide                all
  55 global_canonical  canonical                package-wide                all
  56 global_canonical  canonical                package-wide                all
  57 global_canonical  canonical                package-wide                all
  58 global_canonical  canonical                package-wide                all
  59 global_canonical  canonical                package-wide                all
  60 global_canonical  canonical                package-wide                all
  61 global_canonical  canonical                package-wide                all
  62 global_canonical  canonical                package-wide                all
  63 global_canonical  canonical                package-wide                all
                     alias_of                                     source
  1                      <NA>                             pp_decode docs
  2                      <NA>                             pp_decode docs
  3                      <NA>                             pp_decode docs
  4                      <NA>      pp_decode docs; tapestry_hyperec docs
  5                      <NA>                             pp_decode docs
  6                      <NA>      pp_decode docs; tapestry_hyperec docs
  7                      <NA>      pp_decode docs; tapestry_hyperec docs
  8                      <NA>                             pp_decode docs
  9                      <NA>                             pp_decode docs
  10                     <NA>                             pp_decode docs
  11                    N / J             compare_* docs; pp_design docs
  12         \\lambda_{\\max}                             compare_* docs
  13                      w_i                                    pg docs
  14                        J                notation audit; legacy text
  15                        J   pp_verify docs; separable/disjunct attrs
  16                        N   pp_verify docs; separable/disjunct attrs
  17                      n_j                                    pg docs
  18                      n_j             pp_matrix docs; compare_* docs
  19                      w_i list_designs docs; additional design attrs
  20                      w_i                             pp_verify docs
  21                     <NA>                      classical-design docs
  22                     <NA>                      classical-design docs
  23                     <NA>                      classical-design docs
  24                     <NA>                      classical-design docs
  25                     <NA>                      classical-design docs
  26                     <NA>                      classical-design docs
  27                     <NA>                      classical-design docs
  28                     <NA>                      classical-design docs
  29                     <NA>                      classical-design docs
  30                     <NA>                      classical-design docs
  31                     <NA>                                 hyper docs
  32                     <NA>                                 hyper docs
  33                     <NA>           pp_matrix docs; ppgt_pp vignette
  34                     <NA>           pp_matrix docs; ppgt_pp vignette
  35                     <NA>           pp_matrix docs; ppgt_pp vignette
  36                     <NA>           pp_matrix docs; ppgt_pp vignette
  37                     <NA>                    separable/disjunct docs
  38                     <NA>                    separable/disjunct docs
  39              \\mathbf{z}                notation audit; legacy text
  40 \\widetilde{\\mathbf{y}}                notation audit; legacy text
  41                      w_i                notation audit; legacy text
  42                     <NA>                Notations.tex; package docs
  43                     <NA>                              Notations.tex
  44                     <NA>                              Notations.tex
  45                     <NA>                              Notations.tex
  46                     <NA>                Notations.tex; package docs
  47                     <NA>                Notations.tex; package docs
  48                     <NA>                Notations.tex; package docs
  49                     <NA>                Notations.tex; package docs
  50                     <NA>                Notations.tex; package docs
  51                     <NA>                Notations.tex; package docs
  52                     <NA>                Notations.tex; package docs
  53                     <NA>                Notations.tex; package docs
  54                     <NA>                Notations.tex; package docs
  55                     <NA>                Notations.tex; package docs
  56                     <NA>                Notations.tex; package docs
  57                     <NA>                Notations.tex; package docs
  58                     <NA>                Notations.tex; package docs
  59                     <NA>                              Notations.tex
  60                     <NA>                              Notations.tex
  61                     <NA>                              Notations.tex
  62                     <NA>                Notations.tex; package docs
  63                     <NA>                Notations.tex; package docs
                                                                                                   notes
  1                                                 Use only inside GPSR or related decoder derivations.
  2                                                 Use only inside GPSR or related decoder derivations.
  3                                                   Use only inside decoder-local estimation formulas.
  4                              Use only in decoder-local least-squares or compressed-sensing notation.
  5                                       Do not confuse with measurement-noise notation used elsewhere.
  6  Do not replace package-wide latent-status notation with this symbol outside decoder-local sections.
  7           Use only inside decoder-local sections; package-wide observed outcomes remain \\mathbf{z}.
  8                                                 Use only inside GPSR or related decoder derivations.
  9                                                 Use only inside GPSR or related decoder derivations.
  10                                       Valid only inside GPSR and related decoder-local derivations.
  11                          Alias used in design summaries; the canonical mathematical ratio is N / J.
  12             Alias used in comparison outputs; the underlying quantity is the maximum column weight.
  13                                         Alias used in projective-plane summaries for column weight.
  14    Legacy output label only; do not use plain M as the package-wide symbol for the number of pools.
  15                              Alias used in return metadata; the canonical mathematical symbol is J.
  16                              Alias used in return metadata; the canonical mathematical symbol is N.
  17                                             Alias used in projective-plane summaries for pool size.
  18                                                      Alias used widely in design summaries for n_j.
  19                                                      Alias used widely in design summaries for w_i.
  20                                                     Alias used in PP verification metadata for w_i.
  21                                      Valid only in family-local sections such as projective planes.
  22                                      Valid only in family-local sections such as projective planes.
  23                                                            Valid only in projective-plane sections.
  24                                                                        Valid only in BIBD sections.
  25                                                                        Valid only in BIBD sections.
  26                                                       Valid only in two-dimensional array sections.
  27                                                                       Use only for Dorfman designs.
  28                                                          Interpret only within the selected family.
  29                        Interpret only within the selected family; it is overloaded across families.
  30                                      Use only in classical-design sections such as BIBD or Kirkman.
  31                                                       Use only in HYPER matrix-generation sections.
  32                                                       Use only in HYPER matrix-generation sections.
  33                                                          Interpret only within the selected family.
  34                                                                   Use only in PP / P-BEST sections.
  35                                     Use only in PP / P-BEST sections; API argument name remains nl.
  36                        Interpret only within the selected family; it is overloaded across families.
  37                                                         Use only in separable/disjunct definitions.
  38                                                         Use only in separable/disjunct definitions.
  39         Outside decoder-local derivations, observed pool outcomes should be written as \\mathbf{z}.
  40                                 Do not use uppercase-bold latent infection vectors in package docs.
  41                                                    Use w_i for package-wide column weight notation.
  42                                                                 Primary package-wide matrix symbol.
  43                               Defined in Notations.tex; used rarely outside the notation reference.
  44                                        Defined in Notations.tex for poolwise simulation parameters.
  45                                        Defined in Notations.tex for poolwise simulation parameters.
  46                                           Primary package-wide notation for observed pool outcomes.
  47                                                 Use to define column weights and overlap structure.
  48                                           Use to define pool sizes and latent pool-status formulas.
  49                                        Primary package-wide notation for latent infection statuses.
  50                                             Primary package-wide notation for latent pool statuses.
  51                                          Canonical scalar component of the latent infection vector.
  52                                        Canonical scalar component of the latent pool-status vector.
  53                                                                    Canonical index for individuals.
  54                                                                          Canonical index for pools.
  55                                                         Primary package-wide scalar for pool count.
  56                                            Use with \\mathbf{M}, J, and N in package-wide notation.
  57                                                       Primary package-wide scalar for sample count.
  58                                                      Canonical package-wide notation for pool size.
  59            Defined in Notations.tex but used less often than population-level prevalence arguments.
  60                                       Defined in Notations.tex; used mainly in simulation notation.
  61                                       Defined in Notations.tex; used mainly in simulation notation.
  62                    Canonical package-wide notation for the number of pools containing individual i.
  63                  Use lowercase z_j only in plain-text descriptions when math mode is not available.

notations[
  notations$category == "global_canonical",
  c("symbol", "meaning", "family", "status")
]
                       symbol                                     meaning family
  42              \\mathbf{M}                       pooling design matrix    all
  43      \\mathbf{m}_j^\\top                 row j of the pooling matrix    all
  44            \\mathbf{s}_e          vector of pool-level sensitivities    all
  45            \\mathbf{s}_p          vector of pool-level specificities    all
  46              \\mathbf{z}             vector of observed pool results    all
  47           \\mathcal{J}_i        set of pools containing individual i    all
  48           \\mathcal{P}_j       set of individuals assigned to pool j    all
  49 \\widetilde{\\mathbf{y}}       vector of latent infection indicators    all
  50 \\widetilde{\\mathbf{z}}              vector of latent pool statuses    all
  51         \\widetilde{y}_i latent infection indicator for individual i    all
  52         \\widetilde{z}_j                     latent status of pool j    all
  53                        i                            individual index    all
  54                        j                                  pool index    all
  55                        J                       total number of pools    all
  56                   M_{ji}                 entry of the pooling matrix    all
  57                        N                 total number of individuals    all
  58                      n_j                              size of pool j    all
  59                      p_i      infection probability for individual i    all
  60                  S_{e_j}           pool-level sensitivity for pool j    all
  61                  S_{p_j}           pool-level specificity for pool j    all
  62                      w_i              column weight for individual i    all
  63                      z_j                  observed outcome of pool j    all
        status
  42 canonical
  43 canonical
  44 canonical
  45 canonical
  46 canonical
  47 canonical
  48 canonical
  49 canonical
  50 canonical
  51 canonical
  52 canonical
  53 canonical
  54 canonical
  55 canonical
  56 canonical
  57 canonical
  58 canonical
  59 canonical
  60 canonical
  61 canonical
  62 canonical
  63 canonical
```

The registry is intended to be the programmatic reference point for
notation used in the package documentation. Family-specific symbols such
as PP parameters or decoder-local symbols remain in the registry, but
they are explicitly marked as local rather than package-wide canonical
notation.

## 5 Package Organization

The vignette suite is organized as follows.

1.  `ppgt_theory`: package introduction, notation, design-family tour,
    and a short workflow example
2.  `ppgt_pp`: Polynomial Pools construction, truncation, and decoding
3.  `ppgt_pbest`: PP-based P-BEST, later clinical P-BEST discussion, and
    package-specific boundaries
4.  `ppgt_hyper`: HYPER matrix generation
5.  `ppgt_hyperec`: package-specific experimental HyperEC extension
6.  `ppgt_galois_field`: finite-field arithmetic
7.  `ppgt_separable_disjunct`: separable and disjunct matrix
    constructors
8.  `ppgt_comparisons`: package comparison helpers
9.  `ppgt_diagnostics`: incidence-matrix diagnostics for overlap,
    balance, and interpretation

## 6 Brief Tour Of The Design Families

The package spans several distinct construction families.

- Polynomial Pools: finite-field constructions indexed by \\q\\, \\d\\,
  and \\n_l\\, with sparse binary matrices produced by
  [`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md).
- PP-based P-BEST: the standard package configuration
  `pp_matrix(q = 8, d = 3, nl = 6, N = 384)` corresponds to the 48-pool
  PP-based family used in the bundled P-BEST references.
- HYPER:
  [`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md)
  is an upstream-compatible HYPER matrix generator for \\q = 1\\, \\2\\,
  or \\3\\.
- HyperEC:
  [`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md)
  and
  [`hyper_ec_decode()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_decode.md)
  are package-specific experimental extensions layered on top of a base
  HYPER design.
- Separable and disjunct designs:
  [`SeparableMatrix()`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md)
  and
  [`DisjunctMatrix()`](https://myaseen208.github.io/ppgt/reference/DisjunctMatrix.md)
  provide random and Reed-Solomon-style paths, together with
  verification metadata.
- Comparisons:
  [`compare_all_designs()`](https://myaseen208.github.io/ppgt/reference/compare_all_designs.md)
  and
  [`compare_all_designs_honest()`](https://myaseen208.github.io/ppgt/reference/compare_all_designs_honest.md)
  are package summaries of constructors, not literature results.

At a high level, the package covers both algebraic and combinatorial
matrix families. Polynomial Pools and the PP-based P-BEST family are
driven by finite-field arithmetic over \\\mathrm{GF}(q)\\. HYPER and
HyperEC emphasize balanced low-weight incidence patterns. The separable
and disjunct constructors target exact distinguishability properties for
bounded active sets. The additional design families expose classical
screening layouts such as Dorfman, arrays, BIBDs, hypercubes, Kirkman
systems, and projective planes.

## 7 Paper-Aligned, Compatible, And Package-Specific Scope

The package documentation distinguishes:

- paper-aligned construction claims, such as the PP construction and the
  standard PP-based P-BEST matrix family
- upstream-compatible implementation claims, such as the
  [`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md)
  generator
- package-specific or experimental behavior, such as HyperEC and the
  summary labels used by the comparison helpers

## 8 Workflow Boundary

Adaptive workflows, non-adaptive workflows, retesting scenarios, and
protocol-level decision rules belong above matrix construction.

- Use
  [`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md)
  or a family-specific constructor to obtain \\\mathbf{M}\\.
- Use a decoder to interpret \\\mathbf{z}\\ relative to \\\mathbf{M}\\.
- Use
  [`simulate_group_testing()`](https://myaseen208.github.io/ppgt/reference/simulate_group_testing.md)
  only when the task is to simulate \\\widetilde{\mathbf{z}}\\ and
  \\\mathbf{z}\\ from a fixed design.
- Use
  [`run_testing_workflow()`](https://myaseen208.github.io/ppgt/reference/run_testing_workflow.md)
  only when the task is to coordinate the currently implemented stage-1
  and optional stage-2 protocol logic.
- Use
  [`protocol_summary()`](https://myaseen208.github.io/ppgt/reference/protocol_summary.md)
  only when the task is to summarize an existing workflow object.

This keeps matrix construction honest and keeps workflow behavior
explicit.

## 9 Constructor-Layer Example With `PoolMatrix()`

The constructor wrapper
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md)
belongs only to the matrix-construction layer. A simple deterministic PP
example is enough to show that boundary cleanly.

``` r
M_pool <- PoolMatrix(family = "pp", q = 4, d = 3, nl = 5, N = 64)
M_pool
  20 x 64 sparse Matrix of class "dgCMatrix"
                                                                                 
   [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
   [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
   [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
   [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
   [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
  [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
  [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
  [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
  [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
  [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
  [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
  [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
  [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
  [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
                                                             
   [1,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [2,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [3,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [4,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [5,] . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . .
   [6,] . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . .
   [7,] 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 .
   [8,] . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1
   [9,] . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1
  [10,] . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 .
  [11,] 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . .
  [12,] . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . .
  [13,] . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1
  [14,] 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 .
  [15,] . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . .
  [16,] . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
  [17,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [18,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [19,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
  [20,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1

Y_tilde_pool <- rep(0L, 64)
Y_tilde_pool[c(2, 13)] <- 1L
Y_tilde_pool
   [1] 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

sim_pool <- simulate_group_testing(M_pool, Y_tilde_pool)
sim_pool
  $matrix
  20 x 64 sparse Matrix of class "dgCMatrix"
                                                                                 
   [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
   [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
   [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
   [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
   [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
  [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
  [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
  [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
  [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
  [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
  [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
  [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
  [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
  [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
                                                             
   [1,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [2,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [3,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [4,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [5,] . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . .
   [6,] . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . .
   [7,] 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 .
   [8,] . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1
   [9,] . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1
  [10,] . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 .
  [11,] 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . .
  [12,] . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . .
  [13,] . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1
  [14,] 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 .
  [15,] . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . .
  [16,] . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
  [17,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [18,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [19,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
  [20,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
  
  $Y_tilde
   [1] 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  
  $z_tilde
   [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
  
  $z
   [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
  
  $s_e
  NULL
  
  $s_p
  NULL
  
  $mode
  [1] "noiseless"
sim_pool$z_tilde
   [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
sim_pool$z
   [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0

wf_pool <- run_testing_workflow(
  design = M_pool,
  z = sim_pool$z,
  workflow = "non_adaptive",
  decoder = "pp_decode",
  decoder_args = list(verbose = FALSE),
  stage2 = "none"
)
wf_pool
  $workflow
  [1] "non_adaptive"
  
  $matrix
  20 x 64 sparse Matrix of class "dgCMatrix"
                                                                                 
   [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
   [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
   [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
   [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
   [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
  [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
  [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
  [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
  [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
  [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
  [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
  [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
  [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
  [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
                                                             
   [1,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [2,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [3,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [4,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [5,] . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . .
   [6,] . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . .
   [7,] 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 .
   [8,] . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1
   [9,] . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1
  [10,] . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 .
  [11,] 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . .
  [12,] . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . .
  [13,] . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1
  [14,] 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 .
  [15,] . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . .
  [16,] . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
  [17,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [18,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [19,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
  [20,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
  
  $z
   [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
  
  $stage1
  $stage1$decoder
  [1] "pp_decode"
  
  $stage1$raw
  $stage1$raw$positives
  [1]  2 13
  
  $stage1$raw$n_positives
  [1] 2
  
  $stage1$raw$suspected
  integer(0)
  
  $stage1$raw$n_suspected
  [1] 0
  
  $stage1$raw$candidates
  [1]  2 13
  
  $stage1$raw$n_candidates
  [1] 2
  
  $stage1$raw$error
  [1] 0
  
  $stage1$raw$method
  [1] "COMP (exact)"
  
  
  $stage1$positive
  [1]  2 13
  
  $stage1$suspected
  integer(0)
  
  $stage1$negative
  integer(0)
  
  $stage1$candidates
  [1]  2 13
  
  
  $stage2
  $stage2$method
  [1] "none"
  
  $stage2$planned
  integer(0)
  
  $stage2$results
   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA
  
  $stage2$completed
  [1] TRUE
  
  
  $final_calls
  $final_calls$positive
  [1]  2 13
  
  $final_calls$negative
  integer(0)
  
  $final_calls$unresolved
   [1]  1  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27
  [26] 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52
  [51] 53 54 55 56 57 58 59 60 61 62 63 64
  
  
  attr(,"class")
  [1] "ppgt_workflow"
wf_pool$stage1
  $decoder
  [1] "pp_decode"
  
  $raw
  $raw$positives
  [1]  2 13
  
  $raw$n_positives
  [1] 2
  
  $raw$suspected
  integer(0)
  
  $raw$n_suspected
  [1] 0
  
  $raw$candidates
  [1]  2 13
  
  $raw$n_candidates
  [1] 2
  
  $raw$error
  [1] 0
  
  $raw$method
  [1] "COMP (exact)"
  
  
  $positive
  [1]  2 13
  
  $suspected
  integer(0)
  
  $negative
  integer(0)
  
  $candidates
  [1]  2 13
wf_pool$final_calls
  $positive
  [1]  2 13
  
  $negative
  integer(0)
  
  $unresolved
   [1]  1  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27
  [26] 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52
  [51] 53 54 55 56 57 58 59 60 61 62 63 64

summary_pool <- protocol_summary(wf_pool)
summary_pool
  $workflow
  [1] "non_adaptive"
  
  $n_pools
  [1] 20
  
  $n_individuals
  [1] 64
  
  $n_stage1_positive
  [1] 2
  
  $n_stage1_suspected
  [1] 0
  
  $n_stage1_negative
  [1] 0
  
  $n_stage2_planned
  [1] 0
  
  $n_final_positive
  [1] 2
  
  $n_final_negative
  [1] 0
  
  $n_final_unresolved
  [1] 62
```

In this example:

- [`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md)
  returns only the stage-1 design matrix \\\mathbf{M}\\
- [`simulate_group_testing()`](https://myaseen208.github.io/ppgt/reference/simulate_group_testing.md)
  generates \\\widetilde{\mathbf{z}}\\ and \\\mathbf{z}\\ from the fixed
  design
- [`run_testing_workflow()`](https://myaseen208.github.io/ppgt/reference/run_testing_workflow.md)
  applies the currently implemented non-adaptive workflow skeleton above
  that design
- [`protocol_summary()`](https://myaseen208.github.io/ppgt/reference/protocol_summary.md)
  reports the stored workflow result

The constructor itself still does not choose a workflow, perform
retesting, or run a decoder.

## 10 Incidence-Matrix Diagnostics

Two matrix products are especially useful for reading a pooling design.

The pool-overlap matrix is \\\mathbf{M}\mathbf{M}^\top =
\operatorname{tcrossprod}(\mathbf{M})\\.

Its \\(j, j')\\ entry counts the number of individuals shared by pools
\\\mathcal{P}\_j\\ and \\\mathcal{P}\_{j'}\\. The diagonal entries are
the pool sizes \\n_j\\.

The sample co-occurrence matrix is \\\mathbf{M}^\top \mathbf{M} =
\operatorname{crossprod}(\mathbf{M})\\.

Its \\(i, i')\\ entry counts the number of pools shared by individuals
\\i\\ and \\i'\\. The diagonal entries are the sample coverages \\w_i\\.

These two diagnostics help interpret:

- balance, through the diagonal structure
- overlap patterns, through the off-diagonal structure
- how strongly different pools or individuals are coupled
- why low shared coverage supports \\d\\-separable and \\d\\-disjunct
  interpretations

In particular, the diagonal entries recover the row and column weights
\\ (\mathbf{M}\mathbf{M}^\top)\_{jj} = n_j = \|\mathcal{P}\_j\|
\qquad\text{and}\qquad (\mathbf{M}^\top\mathbf{M})\_{ii} = w_i =
\|\mathcal{J}\_i\|, \\ while the off-diagonal entries measure overlap:
\\ (\mathbf{M}\mathbf{M}^\top)\_{jj'} = \|\mathcal{P}\_j \cap
\mathcal{P}\_{j'}\|, \qquad (\mathbf{M}^\top\mathbf{M})\_{ii'} =
\|\mathcal{J}\_i \cap \mathcal{J}\_{i'}\|. \\ These identities are
useful across the entire package, regardless of whether the design
arises from PP, HYPER, a combinatorial block design, or a verified
separable/disjunct construction.

The full diagnostics discussion is developed in `ppgt_diagnostics`.

## 11 Short End-To-End Example

``` r
M <- pp_matrix(q = 4, d = 3, nl = 5)
M
  20 x 64 sparse Matrix of class "dgCMatrix"
                                                                                 
   [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
   [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
   [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
   [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
   [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
  [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
  [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
  [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
  [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
  [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
  [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
  [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
  [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
  [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
                                                             
   [1,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [2,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [3,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [4,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [5,] . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . .
   [6,] . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . .
   [7,] 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 .
   [8,] . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1
   [9,] . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1
  [10,] . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 .
  [11,] 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . .
  [12,] . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . .
  [13,] . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1
  [14,] 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 .
  [15,] . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . .
  [16,] . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
  [17,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [18,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [19,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
  [20,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1

Y_tilde <- rep(0L, 64)
Y_tilde[c(2, 13)] <- 1L
Y_tilde
   [1] 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

z <- as.integer(as.vector(M %*% Y_tilde) > 0)
z
   [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0

decode <- pp_decode(M, z, verbose = FALSE)
decode
  $positives
  [1]  2 13
  
  $n_positives
  [1] 2
  
  $suspected
  integer(0)
  
  $n_suspected
  [1] 0
  
  $candidates
  [1]  2 13
  
  $n_candidates
  [1] 2
  
  $error
  [1] 0
  
  $method
  [1] "COMP (exact)"

pool_overlap <- Matrix::tcrossprod(M)
pool_overlap
  20 x 20 sparse Matrix of class "dsCMatrix"
                                                                   
   [1,] 16  .  .  .  4  4  4  4  4  4  4  4  4  4  4  4 16  .  .  .
   [2,]  . 16  .  .  4  4  4  4  4  4  4  4  4  4  4  4  . 16  .  .
   [3,]  .  . 16  .  4  4  4  4  4  4  4  4  4  4  4  4  .  . 16  .
   [4,]  .  .  . 16  4  4  4  4  4  4  4  4  4  4  4  4  .  .  . 16
   [5,]  4  4  4  4 16  .  .  .  4  4  4  4  4  4  4  4  4  4  4  4
   [6,]  4  4  4  4  . 16  .  .  4  4  4  4  4  4  4  4  4  4  4  4
   [7,]  4  4  4  4  .  . 16  .  4  4  4  4  4  4  4  4  4  4  4  4
   [8,]  4  4  4  4  .  .  . 16  4  4  4  4  4  4  4  4  4  4  4  4
   [9,]  4  4  4  4  4  4  4  4 16  .  .  .  4  4  4  4  4  4  4  4
  [10,]  4  4  4  4  4  4  4  4  . 16  .  .  4  4  4  4  4  4  4  4
  [11,]  4  4  4  4  4  4  4  4  .  . 16  .  4  4  4  4  4  4  4  4
  [12,]  4  4  4  4  4  4  4  4  .  .  . 16  4  4  4  4  4  4  4  4
  [13,]  4  4  4  4  4  4  4  4  4  4  4  4 16  .  .  .  4  4  4  4
  [14,]  4  4  4  4  4  4  4  4  4  4  4  4  . 16  .  .  4  4  4  4
  [15,]  4  4  4  4  4  4  4  4  4  4  4  4  .  . 16  .  4  4  4  4
  [16,]  4  4  4  4  4  4  4  4  4  4  4  4  .  .  . 16  4  4  4  4
  [17,] 16  .  .  .  4  4  4  4  4  4  4  4  4  4  4  4 16  .  .  .
  [18,]  . 16  .  .  4  4  4  4  4  4  4  4  4  4  4  4  . 16  .  .
  [19,]  .  . 16  .  4  4  4  4  4  4  4  4  4  4  4  4  .  . 16  .
  [20,]  .  .  . 16  4  4  4  4  4  4  4  4  4  4  4  4  .  .  . 16

sample_overlap <- Matrix::crossprod(M)
sample_overlap[1:12, 1:12]
  12 x 12 sparse Matrix of class "dsCMatrix"
                               
   [1,] 5 . . . 2 1 1 1 2 1 1 1
   [2,] . 5 . . 1 2 1 1 1 2 1 1
   [3,] . . 5 . 1 1 2 1 1 1 2 1
   [4,] . . . 5 1 1 1 2 1 1 1 2
   [5,] 2 1 1 1 5 . . . 2 1 1 1
   [6,] 1 2 1 1 . 5 . . 1 2 1 1
   [7,] 1 1 2 1 . . 5 . 1 1 2 1
   [8,] 1 1 1 2 . . . 5 1 1 1 2
   [9,] 2 1 1 1 2 1 1 1 5 . . .
  [10,] 1 2 1 1 1 2 1 1 . 5 . .
  [11,] 1 1 2 1 1 1 2 1 . . 5 .
  [12,] 1 1 1 2 1 1 1 2 . . . 5
```

This short example shows the full package pattern:

1.  build a matrix \\\mathbf{M}\\
2.  simulate latent statuses \\\widetilde{\mathbf{y}}\\
3.  compute observed pool results \\\mathbf{z}\\
4.  decode
5.  inspect overlap structure through `tcrossprod(M)`

## 12 References

- Tan, Y. H. I. (2020). Pooling matrix designs for group testing. *SIAM
  Undergraduate Research Online*, 13, 1-21.
- Shental, N., Levy, S., Wuvshet, V., et al. (2020). Efficient
  high-throughput SARS-CoV-2 testing to detect asymptomatic carriers.
  *Science Advances*, 6(37), eabc5961.
- Zismanov, V., Yelin, I., Klochendler, A., et al. (2024). High capacity
  clinical SARS-CoV-2 molecular testing using combinatorial pooling.
  *Communications Medicine*, 4, Article 121.
- Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing
  for SARS-CoV-2: A flexible, easy-to-implement, and highly efficient
  method. *Nature Communications*, 13, Article 3626.

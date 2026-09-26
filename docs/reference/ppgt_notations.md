# Inspect the Package Notation Registry

Return the package-wide notation registry used by ppgt documentation.
The registry records canonical symbols, family-local notation,
decoder-local notation, display aliases used in summaries, and
deprecated or forbidden aliases that should not be used as package-wide
notation.

## Usage

``` r
ppgt_notations()
```

## Value

A base `data.frame` with one row per notation entry and the columns:

- `symbol`:

  Human-readable symbol label.

- `render`:

  Render-safe LaTeX-ready representation used in the package
  documentation.

- `meaning`:

  Meaning of the symbol in context.

- `category`:

  One of `"global_canonical"`, `"family_local"`, `"decoder_local"`,
  `"display_alias"`, or `"forbidden_alias"`.

- `status`:

  Notation status, such as `"canonical"`, `"local"`, `"alias"`, or
  `"deprecated"`.

- `scope`:

  Where the notation should be used.

- `family`:

  Associated design family, layer, or package-wide scope.

- `alias_of`:

  Canonical target symbol for aliases, or `NA_character_` when not
  applicable.

- `source`:

  Primary package source for the entry.

- `notes`:

  Short notes describing usage boundaries or deprecation status.

## Details

The package follows the notation defined in `Notations.tex`. The core
package-wide notation is: \$\$ \mathbf{M} \in \\0,1\\^{J \times N},
\qquad \widetilde{\mathbf{y}} = (\widetilde{y}\_1, \ldots,
\widetilde{y}\_N)^\top, \qquad \widetilde{\mathbf{z}} =
(\widetilde{z}\_1, \ldots, \widetilde{z}\_J)^\top, \qquad \mathbf{z} =
(z_1, \ldots, z_J)^\top. \$\$ The corresponding pool-status map is \$\$
\widetilde{z}\_j = \mathbb{I}\left(\sum\_{i = 1}^N
M\_{ji}\widetilde{y}\_i \> 0\right), \qquad j \in \\1, \ldots, J\\. \$\$

The registry makes three distinctions explicit:

- package-wide canonical notation used throughout the main
  documentation;

- family-local or decoder-local symbols that are valid only inside those
  sections;

- aliases and deprecated forms that should not replace the canonical
  package notation.

In particular, design-summary fields such as `pool_size`,
`pools_per_sample`, `compression`, `n_pools`, and `n_samples` are
recorded as display aliases rather than as primary mathematical symbols.
Likewise, deprecated forms such as \\\widetilde{\mathbf{Y}}\\, scalar
\\M\\ for the number of pools, or \\n_i\\ as a package-wide
column-weight symbol are tracked explicitly when they would otherwise
introduce ambiguity.

## References

Tan, Y. H. I. (2020). Pooling matrix designs for group testing. *SIAM
Undergraduate Research Online*, 13, 1-21.

Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing
for SARS-CoV-2: A flexible, easy-to-implement, and highly efficient
method. *Nature Communications*, 13, Article 3626.

Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
information theory perspective. *Foundations and Trends in
Communications and Information Theory*, 15(3-4), 196-392.

## See also

[`ppgt-package`](https://myaseen208.github.io/ppgt/reference/ppgt-package.md),
[`PoolMatrix`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`pp_matrix`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md),
[`pp_decode`](https://myaseen208.github.io/ppgt/reference/pp_decode.md),
[`hyper_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
[`SeparableMatrix`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md),
[`DisjunctMatrix`](https://myaseen208.github.io/ppgt/reference/DisjunctMatrix.md)

## Examples

``` r
notations <- ppgt_notations()
notations
#>                        symbol                     render
#> 1        F_\\tau(\\mathbf{z})       F_\\tau(\\mathbf{z})
#> 2                     \\alpha                    \\alpha
#> 3                      \\beta                     \\beta
#> 4          \\hat{\\mathbf{x}}         \\hat{\\mathbf{x}}
#> 5                 \\mathbf{A}                \\mathbf{A}
#> 6               \\mathbf{g}_k              \\mathbf{g}_k
#> 7                 \\mathbf{x}                \\mathbf{x}
#> 8                 \\mathbf{y}                \\mathbf{y}
#> 9                     \\sigma                    \\sigma
#> 10                      \\tau                      \\tau
#> 11                 Lambda_max                 Lambda_max
#> 12                          M                          M
#> 13                compression                compression
#> 14            lines_per_point            lines_per_point
#> 15                    n_pools                    n_pools
#> 16                  n_samples                  n_samples
#> 17            points_per_line            points_per_line
#> 18                  pool_size                  pool_size
#> 19           pools_per_sample           pools_per_sample
#> 20            sample_coverage            sample_coverage
#> 21                    \\alpha                    \\alpha
#> 22                     \\beta                     \\beta
#> 23                    \\gamma                    \\gamma
#> 24                   \\lambda                   \\lambda
#> 25                          b                          b
#> 26                          c                          c
#> 27                          g                          g
#> 28                          k                          k
#> 29                          r                          r
#> 30                          v                          v
#> 31                          m                          m
#> 32                    reorder                    reorder
#> 33                          d                          d
#> 34                  k_{\\max}                  k_{\\max}
#> 35                        n_l                        n_l
#> 36                          q                          q
#> 37 \\mathbf{u}_{\\mathcal{S}} \\mathbf{u}_{\\mathcal{S}}
#> 38               \\mathcal{S}               \\mathcal{S}
#> 39                \\mathbf{y}                \\mathbf{y}
#> 40   \\widetilde{\\mathbf{Y}}   \\widetilde{\\mathbf{Y}}
#> 41                        n_i                        n_i
#> 42                          J                          J
#> 43                     M_{ji}                     M_{ji}
#> 44                          N                          N
#> 45                    S_{e_j}                    S_{e_j}
#> 46                    S_{p_j}                    S_{p_j}
#> 47                \\mathbf{M}                \\mathbf{M}
#> 48        \\mathbf{m}_j^\\top        \\mathbf{m}_j^\\top
#> 49              \\mathbf{s}_e              \\mathbf{s}_e
#> 50              \\mathbf{s}_p              \\mathbf{s}_p
#> 51                \\mathbf{z}                \\mathbf{z}
#> 52             \\mathcal{J}_i             \\mathcal{J}_i
#> 53             \\mathcal{P}_j             \\mathcal{P}_j
#> 54   \\widetilde{\\mathbf{y}}   \\widetilde{\\mathbf{y}}
#> 55   \\widetilde{\\mathbf{z}}   \\widetilde{\\mathbf{z}}
#> 56           \\widetilde{y}_i           \\widetilde{y}_i
#> 57           \\widetilde{z}_j           \\widetilde{z}_j
#> 58                          i                          i
#> 59                          j                          j
#> 60                        n_j                        n_j
#> 61                        p_i                        p_i
#> 62                        w_i                        w_i
#> 63                        z_j                        z_j
#>                                                                               meaning
#> 1                                                             GPSR objective function
#> 2                                                            GPSR step-size parameter
#> 3                                                   GPSR line-search shrinkage factor
#> 4                                         decoder-local estimate of the sparse signal
#> 5                                              decoder-local sensing or design matrix
#> 6                                        decoder-local gradient vector at iteration k
#> 7                                                 decoder-local unknown sparse signal
#> 8                                                       decoder-local response vector
#> 9                                                          GPSR line-search parameter
#> 10                                                      GPSR regularization parameter
#> 11                                        display label for the maximum column weight
#> 12                                       legacy display alias for the number of pools
#> 13                                                  display label for the ratio N / J
#> 14                                   display label for projective-plane column weight
#> 15                                              display label for the number of pools
#> 16                                        display label for the number of individuals
#> 17                                       display label for projective-plane pool size
#> 18                                                        display label for pool size
#> 19                                                    display label for column weight
#> 20                          display label for column weight in PP verification output
#> 21                   projective-coordinate parameter or family-local tuning parameter
#> 22                   projective-coordinate parameter or family-local tuning parameter
#> 23                                                    projective-coordinate parameter
#> 24                                                      BIBD pair-incidence parameter
#> 25                                                         number of blocks in a BIBD
#> 26                                number of columns in a two-dimensional array design
#> 27                                                                  Dorfman pool size
#> 28                                  family-local block size or active-count parameter
#> 29 family-local row count, replication number, or clinical pools-per-sample parameter
#> 30                                                           family-local point count
#> 31                                   HYPER row-count parameter for the stage-1 design
#> 32                                                                 HYPER reorder flag
#> 33                                   family-local dimension or target order parameter
#> 34                                                       PP guaranteed-positive bound
#> 35                                                                PP number of layers
#> 36                             family-local field order or pools-per-sample parameter
#> 37                                     Boolean union pattern induced by an active set
#> 38                                       active set in separable/disjunct definitions
#> 39                             forbidden package-wide alias for observed pool results
#> 40                                  deprecated uppercase-bold latent infection vector
#> 41                                    deprecated package-wide alias for column weight
#> 42                                                              total number of pools
#> 43                                                        entry of the pooling matrix
#> 44                                                        total number of individuals
#> 45                                                  pool-level sensitivity for pool j
#> 46                                                  pool-level specificity for pool j
#> 47                                                              pooling design matrix
#> 48                                                        row j of the pooling matrix
#> 49                                                 vector of pool-level sensitivities
#> 50                                                 vector of pool-level specificities
#> 51                                                    vector of observed pool results
#> 52                                               set of pools containing individual i
#> 53                                              set of individuals assigned to pool j
#> 54                                              vector of latent infection indicators
#> 55                                                     vector of latent pool statuses
#> 56                                        latent infection indicator for individual i
#> 57                                                            latent status of pool j
#> 58                                                                   individual index
#> 59                                                                         pool index
#> 60                                                                     size of pool j
#> 61                                             infection probability for individual i
#> 62                                                     column weight for individual i
#> 63                                                         observed outcome of pool j
#>            category     status                       scope             family
#> 1     decoder_local      local          decoder-local only            decoder
#> 2     decoder_local      local          decoder-local only            decoder
#> 3     decoder_local      local          decoder-local only            decoder
#> 4     decoder_local      local          decoder-local only            decoder
#> 5     decoder_local      local          decoder-local only            decoder
#> 6     decoder_local      local          decoder-local only            decoder
#> 7     decoder_local      local          decoder-local only            decoder
#> 8     decoder_local      local          decoder-local only            decoder
#> 9     decoder_local      local          decoder-local only            decoder
#> 10    decoder_local      local          decoder-local only            decoder
#> 11    display_alias      alias summary and display outputs          summaries
#> 12    display_alias      alias summary and display outputs          summaries
#> 13    display_alias      alias summary and display outputs          summaries
#> 14    display_alias      alias summary and display outputs          summaries
#> 15    display_alias      alias summary and display outputs          summaries
#> 16    display_alias      alias summary and display outputs          summaries
#> 17    display_alias      alias summary and display outputs          summaries
#> 18    display_alias      alias summary and display outputs          summaries
#> 19    display_alias      alias summary and display outputs          summaries
#> 20    display_alias      alias summary and display outputs          summaries
#> 21     family_local      local           family-local only  classical_designs
#> 22     family_local      local           family-local only  classical_designs
#> 23     family_local      local           family-local only  classical_designs
#> 24     family_local      local           family-local only  classical_designs
#> 25     family_local      local           family-local only  classical_designs
#> 26     family_local      local           family-local only  classical_designs
#> 27     family_local      local           family-local only  classical_designs
#> 28     family_local      local           family-local only  classical_designs
#> 29     family_local      local           family-local only  classical_designs
#> 30     family_local      local           family-local only  classical_designs
#> 31     family_local      local           family-local only              hyper
#> 32     family_local      local           family-local only              hyper
#> 33     family_local      local           family-local only           pp/pbest
#> 34     family_local      local           family-local only           pp/pbest
#> 35     family_local      local           family-local only           pp/pbest
#> 36     family_local      local           family-local only           pp/pbest
#> 37     family_local      local           family-local only separable/disjunct
#> 38     family_local      local           family-local only separable/disjunct
#> 39  forbidden_alias deprecated   legacy or forbidden usage             legacy
#> 40  forbidden_alias deprecated   legacy or forbidden usage             legacy
#> 41  forbidden_alias deprecated   legacy or forbidden usage             legacy
#> 42 global_canonical  canonical                package-wide                all
#> 43 global_canonical  canonical                package-wide                all
#> 44 global_canonical  canonical                package-wide                all
#> 45 global_canonical  canonical                package-wide                all
#> 46 global_canonical  canonical                package-wide                all
#> 47 global_canonical  canonical                package-wide                all
#> 48 global_canonical  canonical                package-wide                all
#> 49 global_canonical  canonical                package-wide                all
#> 50 global_canonical  canonical                package-wide                all
#> 51 global_canonical  canonical                package-wide                all
#> 52 global_canonical  canonical                package-wide                all
#> 53 global_canonical  canonical                package-wide                all
#> 54 global_canonical  canonical                package-wide                all
#> 55 global_canonical  canonical                package-wide                all
#> 56 global_canonical  canonical                package-wide                all
#> 57 global_canonical  canonical                package-wide                all
#> 58 global_canonical  canonical                package-wide                all
#> 59 global_canonical  canonical                package-wide                all
#> 60 global_canonical  canonical                package-wide                all
#> 61 global_canonical  canonical                package-wide                all
#> 62 global_canonical  canonical                package-wide                all
#> 63 global_canonical  canonical                package-wide                all
#>                    alias_of                                     source
#> 1                      <NA>                             pp_decode docs
#> 2                      <NA>                             pp_decode docs
#> 3                      <NA>                             pp_decode docs
#> 4                      <NA>                             pp_decode docs
#> 5                      <NA>      pp_decode docs; tapestry_hyperec docs
#> 6                      <NA>                             pp_decode docs
#> 7                      <NA>      pp_decode docs; tapestry_hyperec docs
#> 8                      <NA>      pp_decode docs; tapestry_hyperec docs
#> 9                      <NA>                             pp_decode docs
#> 10                     <NA>                             pp_decode docs
#> 11         \\lambda_{\\max}                             compare_* docs
#> 12                        J                notation audit; legacy text
#> 13                    N / J             compare_* docs; pp_design docs
#> 14                      w_i                                    pg docs
#> 15                        J   pp_verify docs; separable/disjunct attrs
#> 16                        N   pp_verify docs; separable/disjunct attrs
#> 17                      n_j                                    pg docs
#> 18                      n_j             pp_matrix docs; compare_* docs
#> 19                      w_i list_designs docs; additional design attrs
#> 20                      w_i                             pp_verify docs
#> 21                     <NA>                      classical-design docs
#> 22                     <NA>                      classical-design docs
#> 23                     <NA>                      classical-design docs
#> 24                     <NA>                      classical-design docs
#> 25                     <NA>                      classical-design docs
#> 26                     <NA>                      classical-design docs
#> 27                     <NA>                      classical-design docs
#> 28                     <NA>                      classical-design docs
#> 29                     <NA>                      classical-design docs
#> 30                     <NA>                      classical-design docs
#> 31                     <NA>                                 hyper docs
#> 32                     <NA>                                 hyper docs
#> 33                     <NA>           pp_matrix docs; ppgt_pp vignette
#> 34                     <NA>           pp_matrix docs; ppgt_pp vignette
#> 35                     <NA>           pp_matrix docs; ppgt_pp vignette
#> 36                     <NA>           pp_matrix docs; ppgt_pp vignette
#> 37                     <NA>                    separable/disjunct docs
#> 38                     <NA>                    separable/disjunct docs
#> 39              \\mathbf{z}                notation audit; legacy text
#> 40 \\widetilde{\\mathbf{y}}                notation audit; legacy text
#> 41                      w_i                notation audit; legacy text
#> 42                     <NA>                Notations.tex; package docs
#> 43                     <NA>                Notations.tex; package docs
#> 44                     <NA>                Notations.tex; package docs
#> 45                     <NA>                              Notations.tex
#> 46                     <NA>                              Notations.tex
#> 47                     <NA>                Notations.tex; package docs
#> 48                     <NA>                              Notations.tex
#> 49                     <NA>                              Notations.tex
#> 50                     <NA>                              Notations.tex
#> 51                     <NA>                Notations.tex; package docs
#> 52                     <NA>                Notations.tex; package docs
#> 53                     <NA>                Notations.tex; package docs
#> 54                     <NA>                Notations.tex; package docs
#> 55                     <NA>                Notations.tex; package docs
#> 56                     <NA>                Notations.tex; package docs
#> 57                     <NA>                Notations.tex; package docs
#> 58                     <NA>                Notations.tex; package docs
#> 59                     <NA>                Notations.tex; package docs
#> 60                     <NA>                Notations.tex; package docs
#> 61                     <NA>                              Notations.tex
#> 62                     <NA>                Notations.tex; package docs
#> 63                     <NA>                Notations.tex; package docs
#>                                                                                                  notes
#> 1                                        Valid only inside GPSR and related decoder-local derivations.
#> 2                                                 Use only inside GPSR or related decoder derivations.
#> 3                                                 Use only inside GPSR or related decoder derivations.
#> 4                                                   Use only inside decoder-local estimation formulas.
#> 5                              Use only in decoder-local least-squares or compressed-sensing notation.
#> 6                                       Do not confuse with measurement-noise notation used elsewhere.
#> 7  Do not replace package-wide latent-status notation with this symbol outside decoder-local sections.
#> 8           Use only inside decoder-local sections; package-wide observed outcomes remain \\mathbf{z}.
#> 9                                                 Use only inside GPSR or related decoder derivations.
#> 10                                                Use only inside GPSR or related decoder derivations.
#> 11             Alias used in comparison outputs; the underlying quantity is the maximum column weight.
#> 12    Legacy output label only; do not use plain M as the package-wide symbol for the number of pools.
#> 13                          Alias used in design summaries; the canonical mathematical ratio is N / J.
#> 14                                         Alias used in projective-plane summaries for column weight.
#> 15                              Alias used in return metadata; the canonical mathematical symbol is J.
#> 16                              Alias used in return metadata; the canonical mathematical symbol is N.
#> 17                                             Alias used in projective-plane summaries for pool size.
#> 18                                                      Alias used widely in design summaries for n_j.
#> 19                                                      Alias used widely in design summaries for w_i.
#> 20                                                     Alias used in PP verification metadata for w_i.
#> 21                                      Valid only in family-local sections such as projective planes.
#> 22                                      Valid only in family-local sections such as projective planes.
#> 23                                                            Valid only in projective-plane sections.
#> 24                                                                        Valid only in BIBD sections.
#> 25                                                                        Valid only in BIBD sections.
#> 26                                                       Valid only in two-dimensional array sections.
#> 27                                                                       Use only for Dorfman designs.
#> 28                                                          Interpret only within the selected family.
#> 29                        Interpret only within the selected family; it is overloaded across families.
#> 30                                      Use only in classical-design sections such as BIBD or Kirkman.
#> 31                                                       Use only in HYPER matrix-generation sections.
#> 32                                                       Use only in HYPER matrix-generation sections.
#> 33                                                          Interpret only within the selected family.
#> 34                                                                   Use only in PP / P-BEST sections.
#> 35                                     Use only in PP / P-BEST sections; API argument name remains nl.
#> 36                        Interpret only within the selected family; it is overloaded across families.
#> 37                                                         Use only in separable/disjunct definitions.
#> 38                                                         Use only in separable/disjunct definitions.
#> 39         Outside decoder-local derivations, observed pool outcomes should be written as \\mathbf{z}.
#> 40                                 Do not use uppercase-bold latent infection vectors in package docs.
#> 41                                                    Use w_i for package-wide column weight notation.
#> 42                                                         Primary package-wide scalar for pool count.
#> 43                                            Use with \\mathbf{M}, J, and N in package-wide notation.
#> 44                                                       Primary package-wide scalar for sample count.
#> 45                                       Defined in Notations.tex; used mainly in simulation notation.
#> 46                                       Defined in Notations.tex; used mainly in simulation notation.
#> 47                                                                 Primary package-wide matrix symbol.
#> 48                               Defined in Notations.tex; used rarely outside the notation reference.
#> 49                                        Defined in Notations.tex for poolwise simulation parameters.
#> 50                                        Defined in Notations.tex for poolwise simulation parameters.
#> 51                                           Primary package-wide notation for observed pool outcomes.
#> 52                                                 Use to define column weights and overlap structure.
#> 53                                           Use to define pool sizes and latent pool-status formulas.
#> 54                                        Primary package-wide notation for latent infection statuses.
#> 55                                             Primary package-wide notation for latent pool statuses.
#> 56                                          Canonical scalar component of the latent infection vector.
#> 57                                        Canonical scalar component of the latent pool-status vector.
#> 58                                                                    Canonical index for individuals.
#> 59                                                                          Canonical index for pools.
#> 60                                                      Canonical package-wide notation for pool size.
#> 61            Defined in Notations.tex but used less often than population-level prevalence arguments.
#> 62                    Canonical package-wide notation for the number of pools containing individual i.
#> 63                  Use lowercase z_j only in plain-text descriptions when math mode is not available.

notations[notations$family == "pp/pbest", c("symbol", "meaning", "status")]
#>       symbol                                                meaning status
#> 33         d       family-local dimension or target order parameter  local
#> 34 k_{\\max}                           PP guaranteed-positive bound  local
#> 35       n_l                                    PP number of layers  local
#> 36         q family-local field order or pools-per-sample parameter  local

notations[
  notations$category == "display_alias",
  c("symbol", "alias_of", "notes")
]
#>              symbol         alias_of
#> 11       Lambda_max \\lambda_{\\max}
#> 12                M                J
#> 13      compression            N / J
#> 14  lines_per_point              w_i
#> 15          n_pools                J
#> 16        n_samples                N
#> 17  points_per_line              n_j
#> 18        pool_size              n_j
#> 19 pools_per_sample              w_i
#> 20  sample_coverage              w_i
#>                                                                                               notes
#> 11          Alias used in comparison outputs; the underlying quantity is the maximum column weight.
#> 12 Legacy output label only; do not use plain M as the package-wide symbol for the number of pools.
#> 13                       Alias used in design summaries; the canonical mathematical ratio is N / J.
#> 14                                      Alias used in projective-plane summaries for column weight.
#> 15                           Alias used in return metadata; the canonical mathematical symbol is J.
#> 16                           Alias used in return metadata; the canonical mathematical symbol is N.
#> 17                                          Alias used in projective-plane summaries for pool size.
#> 18                                                   Alias used widely in design summaries for n_j.
#> 19                                                   Alias used widely in design summaries for w_i.
#> 20                                                  Alias used in PP verification metadata for w_i.
```

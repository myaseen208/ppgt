# Compare Designs At Multiple Reference Sample Sizes

Generate three side-by-side comparison tables for the package
constructors of P-BEST, Tapestry, HYPER, and the package-specific
experimental
[`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md)
branch. The function reports each design at its chosen reference
configuration and at two common sample sizes.

## Usage

``` r
compare_all_designs_honest()
```

## Value

A named list with four components:

- `standard`:

  A data frame comparing the package's chosen standard constructor
  settings for each design.

- `n256`:

  A data frame comparing all four designs at \\N = 256\\.

- `n384`:

  A data frame comparing all four designs at \\N = 384\\.

- `summary`:

  A named list of short textual summaries describing the branch
  identified by the current implementation as most efficient, lowest
  column weight, or most flexible.

Each comparison data frame has four rows, one per design, and columns
describing the achieved number of individuals \\N\\, the number of pools
\\J\\, the ratio \\J/N\\, the maximum column weight \\\lambda\_{\max} =
\max_i \|\mathcal{J}\_i\|\\, and configuration labels. The
`Configuration` column and the entries of `summary` are package-specific
labels for this comparison helper; they are not paper claims or formal
optimality guarantees.

## Details

Let \\\mathbf{M} \in \\0,1\\^{J \times N}\\ denote a pooling matrix,
where \\N\\ is the number of individuals and \\J\\ is the number of
pools. This function reports the comparison metric \$\$\frac{J}{N}\$\$
together with the maximum column weight \$\$\lambda\_{\max} = \max\_{1
\le i \le N} \|\mathcal{J}\_i\| = \max\_{1 \le i \le N} \sum\_{j=1}^J
M\_{ji}.\$\$ These quantities are extracted directly from the matrices
returned by the package constructors, or from the corresponding
constructor objects when a design returns additional metadata.

The function produces three tables. The `standard` table compares the
specific constructor settings hard-coded as package reference points:
P-BEST at \\N = 384\\, Tapestry at \\N = 256\\, HYPER at \\N = 256\\,
and experimental HYPER-EC at \\N = 256\\. The `n256` table recomputes
all branches at the common target size \\N = 256\\, and the `n384` table
recomputes all branches at the common target size \\N = 384\\. The
purpose is not to prove a universally optimal design, but to make the
dependence of \\J/N\\ and \\\lambda\_{\max}\\ on the chosen operating
size explicit.

For a successful branch, the table entry `m_over_N` is the realized
ratio \\J/N\\. When a branch fails to generate, the corresponding table
entries are `NA`. The experimental HYPER-EC branch should be interpreted
as a package-specific extension rather than a faithful implementation of
the Hong et al. paper design. Likewise, the printed "summary" text is
package-specific reporting for the current helper and should not be
interpreted as a literature result.

## References

Shental, N., Levy, S., Wuvshet, V., et al. (2020). Efficient
high-throughput SARS-CoV-2 testing to detect asymptomatic carriers.
*Science Advances*, 6(37), eabc5961.

Ghosh, S., Agarwal, R., Rehan, M. A., et al. (2021). Tapestry: A
single-round smart pooling technique for COVID-19 testing. *Nature
Communications*, 12, Article 2995.

Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing
for SARS-CoV-2: A flexible, easy-to-implement, and highly efficient
method. *Nature Communications*, 13, Article 3626.

Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
information theory perspective. *Foundations and Trends in
Communications and Information Theory*, 15(3-4), 196-392.

## See also

[`compare_all_designs`](https://myaseen208.github.io/ppgt/reference/compare_all_designs.md),
[`pp_matrix`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md),
[`tapestry_matrix`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md),
[`hyper_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
[`hyper_ec_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md)

## Examples

``` r
comp <- compare_all_designs_honest()
#> ========================================
#> COMPREHENSIVE DESIGN COMPARISON
#> ========================================
#> 1. Generating designs at their STANDARD configurations...
#> 2. Generating all designs at N=256 (common comparison)...
#> 3. Generating all designs at N=384 (P-BEST optimal)...
#> 
#> ========================================
#> COMPARISON COMPLETE
#> ========================================
#> 
#> TABLE 1: STANDARD CONFIGURATIONS
#> (Each design at its package reference N)
#> ----------------------------------------
#>     Design N_Standard  M m_over_N Lambda_max    Configuration
#> 1   P-BEST        384 48    0.125          6 Standard (paper)
#> 2 Tapestry        256 45    0.176          3         Standard
#> 3    HYPER        256 44    0.172          2         Standard
#> 4 HYPER-EC        256 56    0.219          3         Standard
#> 
#> 
#> TABLE 2: ALL DESIGNS AT N=256
#> (Fair to HYPER, Tapestry, HYPER-EC)
#> ----------------------------------------
#>     Design   N  M m_over_N Lambda_max Configuration Rank_Efficiency
#> 1   P-BEST 256 48    0.188          6   Sub-optimal               3
#> 2 Tapestry 256 45    0.176          3       Optimal               1
#> 3    HYPER 256 44    0.172          2       Optimal               1
#> 4 HYPER-EC 256 56    0.219          3       Optimal               4
#> 
#> PACKAGE-SPECIFIC SUMMARY (efficiency): HYPER & Tapestry tied (m/N = 0.172)
#> 
#> 
#> TABLE 3: ALL DESIGNS AT N=384
#> (Fair to P-BEST)
#> ----------------------------------------
#>     Design   N  M m_over_N Lambda_max Configuration Rank_Efficiency
#> 1   P-BEST 384 48    0.125          6       Optimal               1
#> 2 Tapestry 384 96    0.250          3        Scaled               2
#> 3    HYPER 384 72    0.188          2        Scaled               2
#> 4 HYPER-EC 384 87    0.227          3        Scaled               4
#> 
#> PACKAGE-SPECIFIC SUMMARY (efficiency): P-BEST (m/N = 0.125)
#> 
#> 
#> PACKAGE-SPECIFIC FINDINGS:
#> ============================================
#> 1. Package-specific summary: P-BEST is most efficient at N=384 (0.125)
#> 2. Package-specific summary: HYPER & Tapestry tie at N=256 (0.172)
#> 3. Package-specific summary: HYPER-EC trades efficiency for error correction
#> 4. Package-specific summary: different designs excel at different N
#> 
comp
#> $standard
#>     Design N_Standard  M m_over_N Lambda_max    Configuration
#> 1   P-BEST        384 48    0.125          6 Standard (paper)
#> 2 Tapestry        256 45    0.176          3         Standard
#> 3    HYPER        256 44    0.172          2         Standard
#> 4 HYPER-EC        256 56    0.219          3         Standard
#> 
#> $n256
#>     Design   N  M m_over_N Lambda_max Configuration Rank_Efficiency
#> 1   P-BEST 256 48    0.188          6   Sub-optimal               3
#> 2 Tapestry 256 45    0.176          3       Optimal               1
#> 3    HYPER 256 44    0.172          2       Optimal               1
#> 4 HYPER-EC 256 56    0.219          3       Optimal               4
#> 
#> $n384
#>     Design   N  M m_over_N Lambda_max Configuration Rank_Efficiency
#> 1   P-BEST 384 48    0.125          6       Optimal               1
#> 2 Tapestry 384 96    0.250          3        Scaled               2
#> 3    HYPER 384 72    0.188          2        Scaled               2
#> 4 HYPER-EC 384 87    0.227          3        Scaled               4
#> 
#> $summary
#> $summary$best_at_n256
#> [1] "HYPER & Tapestry (tied)"
#> 
#> $summary$best_at_n384
#> [1] "P-BEST"
#> 
#> $summary$best_error_correction
#> [1] "P-BEST (Strong RS)"
#> 
#> $summary$lowest_lambda
#> [1] "HYPER (lambda=2)"
#> 
#> $summary$most_flexible
#> [1] "HYPER, Tapestry, HYPER-EC"
#> 
#> 
names(comp)
#> [1] "standard" "n256"     "n384"     "summary" 
comp$standard
#>     Design N_Standard  M m_over_N Lambda_max    Configuration
#> 1   P-BEST        384 48    0.125          6 Standard (paper)
#> 2 Tapestry        256 45    0.176          3         Standard
#> 3    HYPER        256 44    0.172          2         Standard
#> 4 HYPER-EC        256 56    0.219          3         Standard
comp$n256
#>     Design   N  M m_over_N Lambda_max Configuration Rank_Efficiency
#> 1   P-BEST 256 48    0.188          6   Sub-optimal               3
#> 2 Tapestry 256 45    0.176          3       Optimal               1
#> 3    HYPER 256 44    0.172          2       Optimal               1
#> 4 HYPER-EC 256 56    0.219          3       Optimal               4
comp$n384
#>     Design   N  M m_over_N Lambda_max Configuration Rank_Efficiency
#> 1   P-BEST 384 48    0.125          6       Optimal               1
#> 2 Tapestry 384 96    0.250          3        Scaled               2
#> 3    HYPER 384 72    0.188          2        Scaled               2
#> 4 HYPER-EC 384 87    0.227          3        Scaled               4
comp$summary
#> $best_at_n256
#> [1] "HYPER & Tapestry (tied)"
#> 
#> $best_at_n384
#> [1] "P-BEST"
#> 
#> $best_error_correction
#> [1] "P-BEST (Strong RS)"
#> 
#> $lowest_lambda
#> [1] "HYPER (lambda=2)"
#> 
#> $most_flexible
#> [1] "HYPER, Tapestry, HYPER-EC"
#> 
```

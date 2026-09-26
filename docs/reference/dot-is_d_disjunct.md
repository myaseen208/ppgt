# Check Whether a Pooling Matrix Is d-Disjunct

A binary matrix \\\mathbf{M}\\ (pools \\\times\\ samples) is
\\d\\-**disjunct** if for every set \\S\\ of \\d\\ columns and every
column \\c \notin S\\, the row support of \\c\\ is **not** a subset of
the union of row supports of \\S\\.

Equivalently: every sample can be uniquely identified even if any \\d\\
other samples are simultaneously positive. A \\d\\-disjunct matrix
guarantees zero false negatives for the COMP decoder when at most \\d\\
samples are truly positive.

## Usage

``` r
.is_d_disjunct(M, d)
```

## Arguments

- M:

  Pooling design matrix \\\mathbf{M} = (M\_{ji}) \in \\0,1\\^{J \times
  N}\\, where row \\j \in \\1,\ldots,J\\\\ indexes pool
  \\\mathcal{P}\_j\\ and column \\i \in \\1,\ldots,N\\\\ indexes
  individual \\i\\. Accepted classes: `matrix`, `dgCMatrix`, or any
  `Matrix` subclass.

- d:

  Positive integer giving the disjunctness level: every active set
  \\\mathcal{S}\\ with \\\|\mathcal{S}\| = d\\ must fail to cover any
  outside individual \\i^\star \notin \mathcal{S}\\.

## Value

A single `logical` scalar with the attributes listed above.

## Algorithm

**Exact** (when \\\binom{N}{d} \leq 10^5\\): for every size-\\d\\ active
set \\\mathcal{S} \subseteq \\1,\ldots,N\\\\ with \\\|\mathcal{S}\| =
d\\, compute the pooled union \$\$ \mathbf{u}\_{\mathcal{S}} =
\bigvee\_{i \in \mathcal{S}} \mathbf{M}\_{\cdot i}, \$\$ then for every
individual \\i^\star \notin \mathcal{S}\\ test whether
\\\mathbf{M}\_{\cdot i^\star} \leq \mathbf{u}\_{\mathcal{S}}\\
entry-wise. Equivalently, this checks whether \\\mathcal{J}\_{i^\star}
\subseteq \bigcup\_{i \in \mathcal{S}} \mathcal{J}\_{i}\\. Returns
`FALSE` immediately upon finding a violating pair \\(\mathcal{S},
i^\star)\\.

**Approximate** (when \\\binom{N}{d} \> 10^5\\): draw \\10^5\\ random
pairs \\(\mathcal{S}, i^\star)\\ with \\\|\mathcal{S}\| = d\\ and
\\i^\star \notin \mathcal{S}\\, and test the same containment condition.
Returns `FALSE` as soon as a violation is found; otherwise returns
`TRUE` (conservative, `attr(., "verified") = "approximate"`).

## Attributes on the Return Value

|             |                                                           |
|-------------|-----------------------------------------------------------|
| Attribute   | Content                                                   |
| `d`         | The \\d\\ tested                                          |
| `n_pools`   | \\J\\ = number of pools (rows of \\\mathbf{M}\\)          |
| `n_samples` | \\N\\ = number of individuals (columns of \\\mathbf{M}\\) |
| `verified`  | `"exact"` or `"approximate"`                              |

## References

Kautz, W. H., & Singleton, R. C. (1964). Nonrandom binary superimposed
codes. *IEEE Transactions on Information Theory*, 10(4), 363-377.

## Examples

``` r
# Identity matrix is 1-disjunct: each sample occupies a unique pool,
# so no other column's support can cover it.
I5 <- diag(5)
res <- ppgt:::.is_d_disjunct(I5, d = 1)
stopifnot(isTRUE(res), attr(res, "verified") == "exact")

# A matrix of all ones is NOT 1-disjunct:
# for any column c and set S = {s}, support(c) = all rows is a subset of support(s).
M_ones <- matrix(1L, nrow = 3, ncol = 4)
stopifnot(isFALSE(ppgt:::.is_d_disjunct(M_ones, d = 1)))

# 2-disjunct check on a small PP matrix
# \dontrun{
  M <- pp_matrix(q = 4, d = 3, nl = 5)   # 20 pools, 64 samples
  ppgt:::.is_d_disjunct(M, d = 2)        # exact (choose(64,2) = 2016)
#> [1] FALSE
#> attr(,"d")
#> [1] 2
#> attr(,"n_pools")
#> [1] 20
#> attr(,"n_samples")
#> [1] 64
#> attr(,"verified")
#> [1] "exact"
# }
```

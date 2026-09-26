# Check Whether a Pooling Matrix Is d-Separable

A binary pooling design matrix \\\mathbf{M} \in \\0,1\\^{J \times N}\\
(\\J\\ pools, \\N\\ individuals) is \\d\\-**separable** if no two
distinct active sets \\\mathcal{S}\_1, \mathcal{S}\_2 \subseteq
\\1,\ldots,N\\\\ with \\\|\mathcal{S}\_1\|, \|\mathcal{S}\_2\| \leq d\\
produce the same union activation vector \\\mathbf{u}\_\mathcal{S} =
\bigvee\_{i \in \mathcal{S}} \mathbf{m}\_i\\. Equivalently, the observed
pool-outcome vector \\\mathbf{z}\\ uniquely identifies the active set
whenever at most \\d\\ individuals are positive.

## Usage

``` r
.is_d_separable(M, d)
```

## Arguments

- M:

  Pooling design matrix \\\mathbf{M} = (M\_{ji}) \in \\0,1\\^{J \times
  N}\\, where row \\j \in \\1,\ldots,J\\\\ indexes pool
  \\\mathcal{P}\_j\\ and column \\i \in \\1,\ldots,N\\\\ indexes
  individual \\i\\. Accepted classes: `matrix`, `dgCMatrix`, or any
  `Matrix` subclass.

- d:

  Positive integer giving the largest active-set cardinality
  \\\|\mathcal{S}\| \leq d\\ to verify.

## Value

A single `logical` scalar with the attributes listed above.

## Algorithm

**Exact** (when \\\binom{N}{d} \leq 10^5\\): enumerate every non-empty
active set \\\mathcal{S} \subseteq \\1,\ldots,N\\\\ with
\\\|\mathcal{S}\| \in \\1,\ldots,d\\\\, compute the induced pool-status
vector \$\$ \mathbf{u}\_{\mathcal{S}} = \bigvee\_{i \in \mathcal{S}}
\mathbf{M}\_{\cdot i} \in \\0,1\\^{J}, \$\$ and store a fingerprint of
\\\mathbf{u}\_{\mathcal{S}}\\ in a hash table. A collision means that
two distinct active sets induce the same pool outcome pattern and are
therefore not distinguishable from \\\mathbf{z}\\.

**Approximate** (when \\\binom{N}{d} \> 10^5\\): draw \\10^5\\ random
active-set pairs \\(\mathcal{S}\_1, \mathcal{S}\_2)\\ with
\\\mathcal{S}\_1, \mathcal{S}\_2 \subseteq \\1,\ldots,N\\\\,
\\\|\mathcal{S}\_1\|, \|\mathcal{S}\_2\| \leq d\\, and \\\mathcal{S}\_1
\neq \mathcal{S}\_2\\, and test whether \\\mathbf{u}\_{\mathcal{S}\_1} =
\mathbf{u}\_{\mathcal{S}\_2}\\. Returns `FALSE` as soon as a collision
is found; returns `TRUE` if none are (conservative:
`attr(., "verified") = "approximate"`).

## Attributes on the Return Value

|             |                                                           |
|-------------|-----------------------------------------------------------|
| Attribute   | Content                                                   |
| `d`         | The \\d\\ tested                                          |
| `n_pools`   | \\J\\ = number of pools (rows of \\\mathbf{M}\\)          |
| `n_samples` | \\N\\ = number of individuals (columns of \\\mathbf{M}\\) |
| `verified`  | `"exact"` or `"approximate"`                              |

## References

D'yachkov, A. G., & Rykov, V. V. (1982). Bounds on the length of
disjunctive codes. *Problemy Peredachi Informatsii*, 18(3), 7-13.

## Examples

``` r
# Identity matrix: each individual i belongs to exactly one unique pool j.
# Any two distinct active sets have disjoint pool memberships → always separable.
I5 <- diag(5)
res <- ppgt:::.is_d_separable(I5, d = 2)
stopifnot(isTRUE(res), attr(res, "verified") == "exact")

# Two identical columns (individuals with identical pool membership) are
# indistinguishable → NOT 1-separable.
M_dup <- cbind(c(1, 0, 1), c(1, 0, 1), c(0, 1, 0))
stopifnot(isFALSE(ppgt:::.is_d_separable(M_dup, d = 1)))

# P-BEST: N = 384 individuals, J = 48 pools (approximate for d = 3)
if (FALSE) { # \dontrun{
  M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  ppgt:::.is_d_separable(M, d = 2)   # exact
  ppgt:::.is_d_separable(M, d = 3)   # approximate (choose(384,3) >> 1e5)
} # }
```

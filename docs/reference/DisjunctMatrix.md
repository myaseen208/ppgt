# Construct a \\d\\-Disjunct Pooling Design Matrix

Returns a binary pooling design matrix \\\mathbf{M} \in \\0,1\\^{J
\times N}\\ that is (or is designed to be) \\d\\-disjunct: for every
active set \\\mathcal{S} \subseteq \\1,\ldots,N\\\\ with
\\\|\mathcal{S}\| = d\\ and every outside individual \\i^\star \notin
\mathcal{S}\\, the pool set \\\mathcal{J}\_{i^\star}\\ is not contained
in \\\bigcup\_{i \in \mathcal{S}} \mathcal{J}\_i\\. The constructor
supports automatic method selection, a random Bernoulli design, and a
Reed-Solomon / Kautz-Singleton construction over \\GF(q)\\.

## Usage

``` r
DisjunctMatrix(
  N,
  d,
  M = NULL,
  method = c("auto", "random", "reed-solomon"),
  seed = 42L,
  verify = TRUE
)
```

## Arguments

- N:

  Positive integer; total number of individuals \\N\\.

- d:

  Positive integer giving the disjunctness order \\d\\.

- M:

  Positive integer or `NULL`. This argument requests the number of
  pools, corresponding to \\J\\ in the notation above. If `NULL`, the
  random method uses \\\lceil (d+1)^2 \log_2(N) \rceil\\; the
  Reed-Solomon method ignores this argument and uses the Kautz-Singleton
  value \\J = q\[d(k-1)+1\]\\.

- method:

  Character string specifying the construction rule: `"auto"`,
  `"random"`, or `"reed-solomon"`.

- seed:

  Integer seed used for the random Bernoulli construction and its retry
  sequence. The Reed-Solomon method records the supplied seed without
  using it in construction.

- verify:

  Logical; if `TRUE`, run `.is_d_disjunct` on the returned matrix and
  attach the result as `attr(result, "disjunct_verified")`.

## Value

A binary `dgCMatrix` with \\J = nrow(\mathbf{M})\\ rows and \\N =
ncol(\mathbf{M})\\ columns. The result carries the attributes `d`,
`method_used`, `seed`, and, when `verify = TRUE`, `disjunct_verified`.
The attribute `disjunct_verified` is a single logical value returned by
[`.is_d_disjunct()`](https://myaseen208.github.io/ppgt/reference/dot-is_d_disjunct.md)
and may itself carry the metadata attributes `d`, `n_pools`,
`n_samples`, and `verified` indicating whether the check was exact or
approximate.

## Details

Let \\N\\ denote the number of individuals, let \\J\\ denote the number
of pools, and let \\\mathbf{M} = (M\_{ji}) \in \\0,1\\^{J \times N}\\ be
the pooling design matrix. For each individual \\i\\, let
\\\mathcal{J}\_i = \\j : M\_{ji} = 1\\\\ denote the set of pools
containing individual \\i\\. Then \\\mathbf{M}\\ is \\d\\-disjunct if
\$\$ \mathcal{J}\_{i^\star} \nsubseteq \bigcup\_{i \in \mathcal{S}}
\mathcal{J}\_i \quad \text{for all } \mathcal{S} \subseteq
\\1,\ldots,N\\,\\ \|\mathcal{S}\| = d,\\ i^\star \notin \mathcal{S}.
\$\$ Equivalently, for every \\i^\star \notin \mathcal{S}\\ there exists
a pool \\j \in \mathcal{J}\_{i^\star}\\ with \\j \notin \bigcup\_{i \in
\mathcal{S}} \mathcal{J}\_i\\.

The Reed-Solomon branch uses the same Kautz-Singleton expansion as
[`SeparableMatrix()`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md).
Let \\q\\ be the smallest prime power for which there exists an integer
\\k \geq 1\\ satisfying \\q^k \geq N\\ and \\L = d(k-1)+1 \leq q+1\\.
Associate individual \\i\\ with the polynomial \\f_i(x)\\ over \\GF(q)\\
of degree at most \\k-1\\, then evaluate at \\L\\ projective points and
expand each q-ary symbol into a one-hot block of length \\q\\, giving
\\J = L q\\ pools.

When \\q \> d\\, the Kautz-Singleton construction yields a
\\d\\-disjunct matrix (Kautz and Singleton, 1964). The argument `M` is
ignored for this method because \\J\\ is fixed by \\q\\, \\k\\, and
\\d\\.

If the requested number of pools is not supplied, use the lower bound
\$\$ J \geq (d+1)^2 \log_2(N). \$\$ The random branch uses independent
Bernoulli entries \$\$ M\_{ji} \sim \mathrm{Bernoulli}(1/d). \$\$ Up to
10 seeds are tried until `.is_d_disjunct` returns `TRUE`, or the retry
budget is exhausted.

Uses the Reed-Solomon construction when the corresponding prime-power
parameter satisfies \\q \> d\\; otherwise falls back to the random
Bernoulli construction.

## Note

Every \\d\\-disjunct matrix is also \\d\\-separable, but the converse
need not hold. Longer derivations and comparisons belong in the
separable/disjunct vignette rather than in this constructor help page.

## Construction Methods

- `method = "reed-solomon"`:

  Uses the Reed-Solomon / Kautz-Singleton path and ignores the
  user-supplied `M` argument because \\J\\ is determined by the
  resulting field parameters.

- `method = "random"`:

  Uses an independent Bernoulli design with default size \\\lceil (d +
  1)^2 \log_2(N) \rceil\\ when `M` is not supplied.

- `method = "auto"`:

  Uses the Reed-Solomon path when its prime-power parameter satisfies
  \\q \> d\\; otherwise falls back to the random Bernoulli path.

## References

Kautz, W. H., & Singleton, R. C. (1964). Nonrandom binary superimposed
codes. *IEEE Transactions on Information Theory*, 10(4), 363-377.

D'yachkov, A. G., & Rykov, V. V. (1982). Bounds on the length of
disjunctive codes. *Problemy Peredachi Informatsii*, 18(3), 7-13.

## See also

[`SeparableMatrix`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md)
for the weaker separability guarantee,
[`.is_d_disjunct`](https://myaseen208.github.io/ppgt/reference/dot-is_d_disjunct.md)
for standalone verification,
[`.is_d_separable`](https://myaseen208.github.io/ppgt/reference/dot-is_d_separable.md)
for the related separability property, and
[`pp_matrix`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)
for the polynomial-pool construction.

## Examples

``` r
M <- DisjunctMatrix(N = 50L, d = 2L)
M
#> 20 x 50 sparse Matrix of class "dgCMatrix"
#>                                                                                
#>  [1,] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 . . . . . . . . . . . . . . . . . . . . .
#>  [2,] . . . . . . . . . . . . . . . . 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 . . . . .
#>  [3,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 1 1 1 1 1
#>  [4,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [5,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#>  [6,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#>  [7,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#>  [8,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>  [9,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
#> [10,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
#> [11,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
#> [12,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
#> [13,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
#> [14,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
#> [15,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
#> [16,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
#> [17,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
#> [18,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
#> [19,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
#> [20,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
#>                                
#>  [1,] . . . . . . . . . . . . .
#>  [2,] . . . . . . . . . . . . .
#>  [3,] 1 1 1 1 1 1 1 1 1 1 1 . .
#>  [4,] . . . . . . . . . . . 1 1
#>  [5,] . . . 1 . . . 1 . . . 1 .
#>  [6,] 1 . . . 1 . . . 1 . . . 1
#>  [7,] . 1 . . . 1 . . . 1 . . .
#>  [8,] . . 1 . . . 1 . . . 1 . .
#>  [9,] . . 1 1 . . . . 1 . . . .
#> [10,] . 1 . . 1 . . 1 . . . . .
#> [11,] 1 . . . . 1 . . . . 1 . 1
#> [12,] . . . . . . 1 . . 1 . 1 .
#> [13,] . . 1 . . 1 . 1 . . . . .
#> [14,] . 1 . . . . 1 . 1 . . . .
#> [15,] 1 . . 1 . . . . . 1 . 1 .
#> [16,] . . . . 1 . . . . . 1 . 1
#> [17,] . . . . . 1 . . 1 . . . 1
#> [18,] 1 . . . . . 1 1 . . . 1 .
#> [19,] . 1 . 1 . . . . . . 1 . .
#> [20,] . . 1 . 1 . . . . 1 . . .
attr(M, "d")
#> [1] 2
attr(M, "method_used")
#> [1] "reed-solomon"

M_rs <- DisjunctMatrix(N = 16L, d = 2L, method = "reed-solomon")
M_rs
#> 12 x 16 sparse Matrix of class "dgCMatrix"
#>                                      
#>  [1,] 1 1 1 1 . . . . . . . . . . . .
#>  [2,] . . . . 1 1 1 1 . . . . . . . .
#>  [3,] . . . . . . . . 1 1 1 1 . . . .
#>  [4,] . . . . . . . . . . . . 1 1 1 1
#>  [5,] 1 . . . 1 . . . 1 . . . 1 . . .
#>  [6,] . 1 . . . 1 . . . 1 . . . 1 . .
#>  [7,] . . 1 . . . 1 . . . 1 . . . 1 .
#>  [8,] . . . 1 . . . 1 . . . 1 . . . 1
#>  [9,] 1 . . . . 1 . . . . 1 . . . . 1
#> [10,] . 1 . . 1 . . . . . . 1 . . 1 .
#> [11,] . . 1 . . . . 1 1 . . . . 1 . .
#> [12,] . . . 1 . . 1 . . 1 . . 1 . . .
attr(M_rs, "disjunct_verified")
#> [1] TRUE
#> attr(,"d")
#> [1] 2
#> attr(,"n_pools")
#> [1] 12
#> attr(,"n_samples")
#> [1] 16
#> attr(,"verified")
#> [1] "exact"
```

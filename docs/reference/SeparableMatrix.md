# Construct a \\d\\-Separable Pooling Design Matrix

Returns a binary pooling design matrix \\\mathbf{M} \in \\0,1\\^{J
\times N}\\ that is (or is designed to be) \\d\\-**separable**: for any
two distinct active sets \\\mathcal{S}\_1, \mathcal{S}\_2 \subseteq
\\1,\ldots,N\\\\ with \\\|\mathcal{S}\_1\|, \|\mathcal{S}\_2\| \leq d\\,
the union pool-activation vectors satisfy \\\mathbf{u}\_{\mathcal{S}\_1}
\neq \mathbf{u}\_{\mathcal{S}\_2}\\. The constructor supports automatic
method selection, a random Bernoulli design, and a Reed-Solomon-style
construction over \\GF(q)\\.

## Usage

``` r
SeparableMatrix(
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

  Positive integer giving the separability order \\d\\.

- M:

  Positive integer or `NULL`. This argument requests the number of
  pools, corresponding to \\J\\ in the notation above. If `NULL`, the
  random method uses \\\lceil 2 d \log_2(N) \rceil\\; the Reed-Solomon
  method ignores this argument and uses the Kautz-Singleton value \\J =
  q\[d(k-1)+1\]\\.

- method:

  Character string specifying the construction rule: `"auto"`,
  `"random"`, or `"reed-solomon"`.

- seed:

  Integer seed used for the random Bernoulli construction and its retry
  sequence. The Reed-Solomon method records the supplied seed without
  using it in construction.

- verify:

  Logical; if `TRUE`, run `.is_d_separable` on the returned matrix and
  attach the result as `attr(result, "separable_verified")`.

## Value

A binary `dgCMatrix` with \\J = nrow(\mathbf{M})\\ rows and \\N =
ncol(\mathbf{M})\\ columns. The result carries the attributes `d`,
`method_used`, `seed`, and, when `verify = TRUE`, `separable_verified`.
The attribute `separable_verified` is a single logical value returned by
[`.is_d_separable()`](https://myaseen208.github.io/ppgt/reference/dot-is_d_separable.md)
and may itself carry the metadata attributes `d`, `n_pools`,
`n_samples`, and `verified` indicating whether the check was exact or
approximate.

## Details

Let \\N\\ denote the number of individuals, let \\J\\ denote the number
of pools, and let \\\mathbf{M} = (M\_{ji}) \in \\0,1\\^{J \times N}\\ be
the pooling design matrix, where \\M\_{ji} = 1\\ if and only if
individual \\i \in \\1,\ldots,N\\\\ belongs to pool \\\mathcal{P}\_j\\.
For an active set \\\mathcal{S} \subseteq \\1,\ldots,N\\\\ with
\\\|\mathcal{S}\| \leq d\\, define the induced pool-status vector \$\$
\mathbf{u}\_{\mathcal{S}} = \bigvee\_{i \in \mathcal{S}}
\mathbf{M}\_{\cdot i} \in \\0,1\\^{J}. \$\$ The matrix \\\mathbf{M}\\ is
\\d\\-separable if \$\$ \mathbf{u}\_{\mathcal{S}\_1} =
\mathbf{u}\_{\mathcal{S}\_2} \\\Longrightarrow\\ \mathcal{S}\_1 =
\mathcal{S}\_2 \quad \text{for all } \mathcal{S}\_1, \mathcal{S}\_2
\subseteq \\1,\ldots,N\\ \text{ with } \|\mathcal{S}\_1\|,
\|\mathcal{S}\_2\| \leq d. \$\$

## Reed-Solomon Construction (9LPR5XxA8tN6Y6t3TXvKwG01TAgq7elb-13-)

Let \\q\\ be the smallest prime power for which there exists an integer
\\k \geq 1\\ such that \\q^k \geq N\\ and \\L = d(k-1)+1 \leq q+1\\.
Associate individual \\i\\ with the polynomial \\f_i(x)\\ of degree at
most \\k-1\\ whose coefficient vector is the base-\\q\\ expansion of
\\i-1\\: \$\$ f_i(x) = c_0(i) + c_1(i)x + \cdots + c\_{k-1}(i)x^{k-1},
\quad i - 1 = \sum\_{\ell=0}^{k-1} c\_\ell(i) q^\ell. \$\$ Choose \\L =
d(k-1)+1\\ distinct evaluation points \\x_1, \ldots, x_L\\ in the
projective line over \\GF(q)\\ and index the pools by \\(\ell, v)\\ with
\\\ell \in \\1,\ldots,L\\\\ and \\v \in GF(q)\\. Then \\J = L q\\ and
\$\$ M\_{ji} = 1 \\\Longleftrightarrow\\ f_i(x\_\ell) = v, \quad j =
(\ell-1)q + v + 1. \$\$ Because any two Reed-Solomon codewords agree in
at most \\k-1\\ evaluation positions, the choice \\L \> d(k-1)\\
satisfies the Kautz-Singleton disjunctness criterion, hence the
resulting binary design is \\d\\-disjunct and therefore \\d\\-separable.
The argument `M` is ignored for this method because \\J\\ is fixed by
\\q\\, \\k\\, and \\d\\.

## Random Bernoulli Construction (9LPR5XxA8tN6Y6t3TXvKwG01TAgq7elb-40-)

If the requested number of pools is not supplied, use the theoretical
lower bound \$\$ J \geq 2 \binom{d}{1} \log_2(N) = 2 d \log_2(N). \$\$
The random construction uses independent Bernoulli entries \$\$ M\_{ji}
\sim \mathrm{Bernoulli}(p), \quad p = \frac{(d \log N)^{1/d}}{N}. \$\$
Up to 10 seeds are tried until `.is_d_separable` returns `TRUE`, or the
retry budget is exhausted.

## Automatic Selection (9LPR5XxA8tN6Y6t3TXvKwG01TAgq7elb-45-)

Uses the Reed-Solomon construction when \\d \leq 3\\; otherwise it falls
back to the random Bernoulli construction.

## References

D'yachkov, A. G., & Rykov, V. V. (1982). Bounds on the length of
disjunctive codes. *Problemy Peredachi Informatsii*, 18(3), 7-13.

Kautz, W. H., & Singleton, R. C. (1964). Nonrandom binary superimposed
codes. *IEEE Transactions on Information Theory*, 10(4), 363-377.

## See also

[`.is_d_separable`](https://myaseen208.github.io/ppgt/reference/dot-is_d_separable.md)
for standalone verification,
[`.is_d_disjunct`](https://myaseen208.github.io/ppgt/reference/dot-is_d_disjunct.md)
for the stronger disjunctness property,
[`DisjunctMatrix`](https://myaseen208.github.io/ppgt/reference/DisjunctMatrix.md)
for the stronger exported constructor, and
[`pp_matrix`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)
for the P-BEST polynomial-pool design.

## Examples

``` r
M <- SeparableMatrix(N = 50L, d = 2L)
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

M_rs <- SeparableMatrix(N = 16L, d = 2L, method = "reed-solomon")
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
attr(M_rs, "separable_verified")
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

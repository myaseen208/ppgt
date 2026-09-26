# Separable and Disjunct Matrix Families

## 1 Abstract

This vignette covers the exported constructors
[`SeparableMatrix()`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md)
and
[`DisjunctMatrix()`](https://myaseen208.github.io/ppgt/reference/DisjunctMatrix.md).

A binary matrix \\\mathbf{M} \in \\0,1\\^{J \times N}\\ is
\\d\\-separable if every two distinct active sets \\\mathcal{S}\_1,
\mathcal{S}\_2 \subseteq \\1, \ldots, N\\\\ with \\\|\mathcal{S}\_1\|
\leq d\\ and \\\|\mathcal{S}\_2\| \leq d\\ induce distinct pooled
outcomes. Equivalently, the pooled outcome vector
\\\widetilde{\mathbf{z}}\\ determines the active set uniquely whenever
at most \\d\\ individuals are positive.

A matrix is \\d\\-disjunct if for every active set \\\mathcal{S}
\subseteq \\1, \ldots, N\\\\ with \\\|\mathcal{S}\| \leq d\\ and every
outside individual \\i^\star \notin \mathcal{S}\\, there is at least one
pool containing \\i^\star\\ that contains none of the individuals in
\\\mathcal{S}\\. This is stronger than \\d\\-separability.

Using the pool index sets \\\mathcal{J}\_i = \\j : M\_{ji} = 1\\\\, the
\\d\\-disjunct condition can be written as \\ \mathcal{J}\_{i^\star}
\nsubseteq \bigcup\_{i \in \mathcal{S}} \mathcal{J}\_i \qquad \text{for
every } \|\mathcal{S}\| \le d \text{ and } i^\star \notin \mathcal{S}.
\\ That is the structural reason a \\d\\-disjunct matrix permits exact
elimination of nondefective individuals from pooled positives.

## 2 Union Maps And Distinguishability

For any active set \\\mathcal{S} \subseteq \\1, \ldots, N\\\\, define
the induced pool-status pattern \\ \mathbf{u}\_{\mathcal{S}} =
\bigvee\_{i \in \mathcal{S}} \mathbf{M}\_{\cdot i} \in \\0,1\\^J. \\
This notation is local to the separable/disjunct theory and is recorded
that way in
[`ppgt_notations()`](https://myaseen208.github.io/ppgt/reference/ppgt_notations.md).

The two core notions can then be written compactly:

- \\d\\-separable means that \\\mathbf{u}\_{\mathcal{S}\_1} \ne
  \mathbf{u}\_{\mathcal{S}\_2}\\ whenever \\\mathcal{S}\_1 \ne
  \mathcal{S}\_2\\ and \\\|\mathcal{S}\_1\|, \|\mathcal{S}\_2\| \le d\\.
- \\d\\-disjunct means that for every \\\mathcal{S} \subseteq
  \\1,\ldots,N\\\\ with \\\|\mathcal{S}\| \le d\\ and every outside
  individual \\i^\star \notin \mathcal{S}\\, \\\mathbf{M}\_{\cdot
  i^\star} \nleq \mathbf{u}\_{\mathcal{S}}\\ entrywise.

So \\d\\-separability asks whether bounded active sets have distinct
pooled signatures, while \\d\\-disjunctness asks for a stronger
one-vs-set isolation property.

## 3 Separable Construction Paths

[`SeparableMatrix()`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md)
supports:

- an automatic path
- a random construction
- a Reed-Solomon-style construction

## 4 A Separable Design

``` r
M_sep <- SeparableMatrix(N = 16L, d = 2L, method = "reed-solomon")
M_sep
  12 x 16 sparse Matrix of class "dgCMatrix"
                                       
   [1,] 1 1 1 1 . . . . . . . . . . . .
   [2,] . . . . 1 1 1 1 . . . . . . . .
   [3,] . . . . . . . . 1 1 1 1 . . . .
   [4,] . . . . . . . . . . . . 1 1 1 1
   [5,] 1 . . . 1 . . . 1 . . . 1 . . .
   [6,] . 1 . . . 1 . . . 1 . . . 1 . .
   [7,] . . 1 . . . 1 . . . 1 . . . 1 .
   [8,] . . . 1 . . . 1 . . . 1 . . . 1
   [9,] 1 . . . . 1 . . . . 1 . . . . 1
  [10,] . 1 . . 1 . . . . . . 1 . . 1 .
  [11,] . . 1 . . . . 1 1 . . . . 1 . .
  [12,] . . . 1 . . 1 . . 1 . . 1 . . .

dim(M_sep)
  [1] 12 16
M_sep[1:8, 1:16]
  8 x 16 sparse Matrix of class "dgCMatrix"
                                      
  [1,] 1 1 1 1 . . . . . . . . . . . .
  [2,] . . . . 1 1 1 1 . . . . . . . .
  [3,] . . . . . . . . 1 1 1 1 . . . .
  [4,] . . . . . . . . . . . . 1 1 1 1
  [5,] 1 . . . 1 . . . 1 . . . 1 . . .
  [6,] . 1 . . . 1 . . . 1 . . . 1 . .
  [7,] . . 1 . . . 1 . . . 1 . . . 1 .
  [8,] . . . 1 . . . 1 . . . 1 . . . 1
attr(M_sep, "method_used")
  [1] "reed-solomon"
attr(M_sep, "separable_verified")
  [1] TRUE
  attr(,"d")
  [1] 2
  attr(,"n_pools")
  [1] 12
  attr(,"n_samples")
  [1] 16
  attr(,"verified")
  [1] "exact"
```

## 5 Disjunct Construction Paths

[`DisjunctMatrix()`](https://myaseen208.github.io/ppgt/reference/DisjunctMatrix.md)
supports the same broad path structure, but targets the stronger
\\d\\-disjunct property.

## 6 A Disjunct Design

``` r
M_disj <- DisjunctMatrix(N = 16L, d = 2L, method = "reed-solomon")
M_disj
  12 x 16 sparse Matrix of class "dgCMatrix"
                                       
   [1,] 1 1 1 1 . . . . . . . . . . . .
   [2,] . . . . 1 1 1 1 . . . . . . . .
   [3,] . . . . . . . . 1 1 1 1 . . . .
   [4,] . . . . . . . . . . . . 1 1 1 1
   [5,] 1 . . . 1 . . . 1 . . . 1 . . .
   [6,] . 1 . . . 1 . . . 1 . . . 1 . .
   [7,] . . 1 . . . 1 . . . 1 . . . 1 .
   [8,] . . . 1 . . . 1 . . . 1 . . . 1
   [9,] 1 . . . . 1 . . . . 1 . . . . 1
  [10,] . 1 . . 1 . . . . . . 1 . . 1 .
  [11,] . . 1 . . . . 1 1 . . . . 1 . .
  [12,] . . . 1 . . 1 . . 1 . . 1 . . .

dim(M_disj)
  [1] 12 16
M_disj[1:8, 1:16]
  8 x 16 sparse Matrix of class "dgCMatrix"
                                      
  [1,] 1 1 1 1 . . . . . . . . . . . .
  [2,] . . . . 1 1 1 1 . . . . . . . .
  [3,] . . . . . . . . 1 1 1 1 . . . .
  [4,] . . . . . . . . . . . . 1 1 1 1
  [5,] 1 . . . 1 . . . 1 . . . 1 . . .
  [6,] . 1 . . . 1 . . . 1 . . . 1 . .
  [7,] . . 1 . . . 1 . . . 1 . . . 1 .
  [8,] . . . 1 . . . 1 . . . 1 . . . 1
attr(M_disj, "method_used")
  [1] "reed-solomon"
attr(M_disj, "disjunct_verified")
  [1] TRUE
  attr(,"d")
  [1] 2
  attr(,"n_pools")
  [1] 12
  attr(,"n_samples")
  [1] 16
  attr(,"verified")
  [1] "exact"
```

## 7 Relationship Between The Two

Every `d`-disjunct design is also `d`-separable, but the converse need
not hold.

The inclusion is strict in general. A \\d\\-separable design only
requires distinct active sets of size at most \\d\\ to induce distinct
pooled outcomes, whereas a \\d\\-disjunct design requires every outside
individual to retain a pool not covered by any active set of size at
most \\d\\.

Equivalently, if \\ \mathbf{u}\_{\mathcal{S}} = \bigvee\_{i \in
\mathcal{S}} \mathbf{M}\_{\cdot i}, \\ then \\d\\-disjunctness requires
an outside individual \\i^\star\\ to have at least one row with \\ M\_{j
i^\star} = 1 \qquad\text{and}\qquad u\_{\mathcal{S}, j} = 0. \\ That
extra uncovered row is what COMP exploits when a disjunct design is used
as a decoder-friendly construction.

``` r
Matrix::colSums(M_sep)
   [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
Matrix::colSums(M_disj)
   [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
```

## 8 Overlap Diagnostics For Separable And Disjunct Designs

The sample co-occurrence matrix \\ \mathbf{M}^\top \mathbf{M} \\ helps
explain why these constructions work. Its diagonal is the column-weight
vector \\w_i\\, and its off-diagonal entries measure how many pools are
shared by pairs of individuals. Small off-diagonal overlap leaves more
rows available to separate one active pattern from another.

``` r
sep_overlap <- Matrix::crossprod(M_sep)
sep_overlap[1:12, 1:12]
  12 x 12 sparse Matrix of class "dsCMatrix"
                               
   [1,] 3 1 1 1 1 1 . . 1 . 1 .
   [2,] 1 3 1 1 1 1 . . . 1 . 1
   [3,] 1 1 3 1 . . 1 1 1 . 1 .
   [4,] 1 1 1 3 . . 1 1 . 1 . 1
   [5,] 1 1 . . 3 1 1 1 1 . . 1
   [6,] 1 1 . . 1 3 1 1 . 1 1 .
   [7,] . . 1 1 1 1 3 1 . 1 1 .
   [8,] . . 1 1 1 1 1 3 1 . . 1
   [9,] 1 . 1 . 1 . . 1 3 1 1 1
  [10,] . 1 . 1 . 1 1 . 1 3 1 1
  [11,] 1 . 1 . . 1 1 . 1 1 3 1
  [12,] . 1 . 1 1 . . 1 1 1 1 3

disj_overlap <- Matrix::crossprod(M_disj)
disj_overlap[1:12, 1:12]
  12 x 12 sparse Matrix of class "dsCMatrix"
                               
   [1,] 3 1 1 1 1 1 . . 1 . 1 .
   [2,] 1 3 1 1 1 1 . . . 1 . 1
   [3,] 1 1 3 1 . . 1 1 1 . 1 .
   [4,] 1 1 1 3 . . 1 1 . 1 . 1
   [5,] 1 1 . . 3 1 1 1 1 . . 1
   [6,] 1 1 . . 1 3 1 1 . 1 1 .
   [7,] . . 1 1 1 1 3 1 . 1 1 .
   [8,] . . 1 1 1 1 1 3 1 . . 1
   [9,] 1 . 1 . 1 . . 1 3 1 1 1
  [10,] . 1 . 1 . 1 1 . 1 3 1 1
  [11,] 1 . 1 . . 1 1 . 1 1 3 1
  [12,] . 1 . 1 1 . . 1 1 1 1 3
```

## 9 Verification Behavior

Both constructors attach verification metadata when verification is
enabled.

``` r
attributes(attr(M_sep, "separable_verified"))
  $d
  [1] 2
  
  $n_pools
  [1] 12
  
  $n_samples
  [1] 16
  
  $verified
  [1] "exact"
attributes(attr(M_disj, "disjunct_verified"))
  $d
  [1] 2
  
  $n_pools
  [1] 12
  
  $n_samples
  [1] 16
  
  $verified
  [1] "exact"
```

## 10 Package Scope

The package also contains internal verification helpers for separability
and disjunctness. This vignette focuses on the exported constructors and
their returned sparse matrices.

## 11 References

- Kautz, W. H., & Singleton, R. C. (1964). Nonrandom binary superimposed
  codes. *IEEE Transactions on Information Theory*, 10(4), 363-377.
- D’yachkov, A. G., & Rykov, V. V. (1982). Bounds on the length of
  disjunctive codes. *Problemy Peredachi Informatsii*, 18(3), 7-13.

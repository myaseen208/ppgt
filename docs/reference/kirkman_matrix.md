# Generate A Kirkman Triple System Pooling Matrix

Construct the incidence matrix of a Kirkman triple system (KTS): a
resolvable Steiner triple system for pooled testing, in which the
triples (blocks) can be partitioned into parallel classes ("days") that
each cover every individual exactly once.

## Usage

``` r
kirkman_matrix(v)
```

## Arguments

- v:

  Integer scalar giving the number of individuals. A Kirkman triple
  system exists if and only if \\v \equiv 3 \pmod 6\\ and \\v \geq 3\\
  (Ray-Chaudhuri & Wilson, 1971); any other value is an error.

## Value

A sparse binary matrix of class `dgCMatrix` with dimensions \\J \times
N\\, where \\N = v\\ and \\J = v(v-1)/6\\. Each row has size three and
each column sum is \\(v-1)/2\\. The attached `design` attribute is a
named list with fields `type`, `v`, `b`, `r`, `k`, and `parallel_class`
(an integer vector of length \\J\\ giving each row's parallel-class
label, \\1,\ldots,(v-1)/2\\). Use
[`is_kts`](https://myaseen208.github.io/ppgt/reference/is_kts.md) to
independently verify that a returned matrix satisfies all four defining
KTS properties.

## Details

A Kirkman triple system is a Steiner triple system whose blocks can be
partitioned into parallel classes. In the returned incidence matrix
\\\mathbf{M} \in \\0,1\\^{J \times N}\\, the columns index the \\N=v\\
individuals and the rows index the triples. The defining combinatorial
relations are \$\$n_j = \|\mathcal{P}\_j\| = 3 \quad \text{for all }
j,\$\$ \$\$w_i = \|\mathcal{J}\_i\| = \frac{v-1}{2} \quad \text{for all
} i,\$\$ every unordered pair of distinct individuals appears together
in exactly one pool (so \\J = v(v-1)/6\\), and the rows can be
partitioned into \\(v-1)/2\\ parallel classes of \\v/3\\ rows each,
every class covering all \\v\\ individuals exactly once. The necessary
and sufficient existence condition is \$\$v \equiv 3 \pmod 6, \quad v
\geq 3.\$\$

## Construction

Three dedicated, deterministic (non-search) constructions are used
whenever \\v\\ matches their pattern:

- \\v = 3^n\\ (\\n = 1, 2, 3, 4, \ldots\\, i.e. \\v \in \\3, 9, 27, 81,
  243, \ldots\\\\): the affine geometry \\\mathrm{AG}(n,3)\\
  construction. Points are the \\3^n\\ vectors of \\\mathrm{GF}(3)^n\\.
  For each of the \\(3^n-1)/2\\ canonical line directions \\d\\ (nonzero
  vectors of \\\mathrm{GF}(3)^n\\ taken up to the equivalence \\d \sim
  -d\\), the lines \\\\p, p+d, p+2d\\\\ for \\p\\ ranging over
  \\\mathrm{GF}(3)^n\\ partition all \\3^n\\ points into \\3^{n-1}\\
  parallel lines, giving one parallel class per direction (Colbourn &
  Rosa, 1999, Section 3.3). This is the standard direct construction of
  a resolvable Steiner triple system from an affine space over
  \\\mathrm{GF}(3)\\; it covers \\v = 3\\ (the trivial system) and \\v =
  9\\ as special cases, in addition to \\v = 27\\ and \\v = 81\\.

- \\v = 15\\: the construction packs the 35 lines of
  \\\mathrm{PG}(3,2)\\ into 7 pairwise line-disjoint spreads (Colbourn &
  Dinitz, 2007, "Kirkman systems"), the classical solution to Kirkman's
  1850 schoolgirl problem.

For every other admissible \\v\\ (e.g. \\v = 21, 33, 39, \ldots\\, and
also \\v = 63\\, since no deterministic construction is currently
implemented for it despite \\63 = 2^6 - 1\\ suggesting a
\\\mathrm{PG}(5,2)\\-spread approach analogous to \\v = 15\\ – packing
all 651 lines of \\\mathrm{PG}(5,2)\\ into 31 disjoint spreads is a
substantially harder combinatorial search than the \\v = 15\\ case and
was not solved within the current implementation), the constructor falls
back to a general randomized search: it generates a pool of candidate
parallel classes (uniformly random partitions of the \\v\\ points into
triples) and performs an exact-cover backtracking search for \\(v-1)/2\\
of them that are pairwise pair-disjoint and jointly cover every pair
exactly once. Every constructed system, from every code path, is
independently re-verified with
[`is_kts`](https://myaseen208.github.io/ppgt/reference/is_kts.md) before
being returned; the function errors instead of ever returning an object
that fails that check.

This search-based general fallback is empirically unreliable once \\v\\
moves beyond the dedicated cases above: it is not guaranteed to find a
system within its computational budget even for \\v = 21\\, and is not
tuned for larger \\v\\ such as 33 or 63. This is a limitation of the
current search heuristic, not evidence that a KTS fails to exist – by
Ray-Chaudhuri & Wilson (1971), KTS(v) exists for every \\v \equiv 3
\pmod 6\\. If the search exhausts its attempt budget, `kirkman_matrix()`
stops with an informative error naming the values of \\v\\ that are
guaranteed fast: any power of 3 (\\3, 9, 27, 81, \ldots\\), and 15.

## References

Ray-Chaudhuri, D. K., & Wilson, R. M. (1971). Solution of Kirkman's
schoolgirl problem. *Proceedings of Symposia in Pure Mathematics*, 19,
187-203.

Colbourn, C. J., & Rosa, A. (1999). *Triple Systems*. Oxford University
Press.

Colbourn, C. J., & Dinitz, J. H. (Eds.). (2007). *Handbook of
Combinatorial Designs* (2nd ed.). Chapman & Hall/CRC.

## See also

[`is_kts`](https://myaseen208.github.io/ppgt/reference/is_kts.md) to
verify the defining properties of a returned (or any other) putative KTS
incidence matrix;
[`bibd_matrix`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md),
[`pg_matrix`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md),
[`compare_designs`](https://myaseen208.github.io/ppgt/reference/compare_designs.md)

Other pooling_designs:
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md),
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`additional_designs`](https://myaseen208.github.io/ppgt/reference/additional_designs.md),
[`array_matrix()`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
[`array_matrix_3d()`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md),
[`bibd_matrix()`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md),
[`compare_designs()`](https://myaseen208.github.io/ppgt/reference/compare_designs.md),
[`dorfman_matrix()`](https://myaseen208.github.io/ppgt/reference/dorfman_matrix.md),
[`dorfman_optimal_g()`](https://myaseen208.github.io/ppgt/reference/dorfman_optimal_g.md),
[`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md),
[`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
[`hypercube_matrix()`](https://myaseen208.github.io/ppgt/reference/hypercube_matrix.md),
[`is_kts()`](https://myaseen208.github.io/ppgt/reference/is_kts.md),
[`list_designs()`](https://myaseen208.github.io/ppgt/reference/list_designs.md),
[`pbest_clinical_matrix()`](https://myaseen208.github.io/ppgt/reference/pbest_clinical_matrix.md),
[`pg_matrix()`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md),
[`regular_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/regular_pooling_matrix.md),
[`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)

## Examples

``` r
M9 <- kirkman_matrix(9)
M9
#> 12 x 9 sparse Matrix of class "dgCMatrix"
#>                        
#>  [1,] 1 1 1 . . . . . .
#>  [2,] . . . 1 1 1 . . .
#>  [3,] . . . . . . 1 1 1
#>  [4,] 1 . . 1 . . 1 . .
#>  [5,] . 1 . . 1 . . 1 .
#>  [6,] . . 1 . . 1 . . 1
#>  [7,] 1 . . . 1 . . . 1
#>  [8,] . 1 . . . 1 1 . .
#>  [9,] . . 1 1 . . . 1 .
#> [10,] 1 . . . . 1 . 1 .
#> [11,] . 1 . 1 . . . . 1
#> [12,] . . 1 . 1 . 1 . .
Matrix::rowSums(M9)
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3
Matrix::colSums(M9)
#> [1] 4 4 4 4 4 4 4 4 4
is_kts(M9)
#> [1] TRUE
#> attr(,"reasons")
#> character(0)

# v = 15: the classic Kirkman schoolgirl problem
M15 <- kirkman_matrix(15)
is_kts(M15)
#> [1] TRUE
#> attr(,"reasons")
#> character(0)
attr(M15, "design")$parallel_class
#>  [1] 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 7 7 7 7 7

# v = 27 and v = 81: AG(n,3), also fast and deterministic
is_kts(kirkman_matrix(27))
#> [1] TRUE
#> attr(,"reasons")
#> character(0)
is_kts(kirkman_matrix(81))
#> [1] TRUE
#> attr(,"reasons")
#> character(0)

# Invalid v (not == 3 mod 6) is an error, not a silent NULL/warning
tryCatch(kirkman_matrix(10), error = function(e) conditionMessage(e))
#> [1] "A Kirkman triple system does not exist for v = 10: v must satisfy v == 3 (mod 6) (e.g. 3, 9, 15, 21, 27, ...)."
```

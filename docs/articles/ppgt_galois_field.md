# GF(q) Arithmetic and Finite-Field Objects

## 1 Abstract

The Polynomial Pools design family uses arithmetic in `GF(q)`, where `q`
must be a prime power. This vignette covers the package finite-field
utilities:
[`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
[`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
and [`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md).

Let \\q = p^n\\, where \\p\\ is prime and \\n \geq 1\\. The package
finite-field constructor builds a `galois_field` object encoding
addition and multiplication in \\GF(q)\\.

For prime fields, arithmetic is ordinary modular arithmetic. For
extension fields, the package represents field elements by integer
labels and uses lookup tables induced by an irreducible polynomial
model. Those lookup tables are what the PP constructor ultimately uses
when it evaluates polynomial relations over \\\mathrm{GF}(q)\\.

## 2 Prime And Prime-Power Checks

``` r
is_prime(7)
  [1] TRUE
is_prime(9)
  [1] FALSE

pp8 <- is_prime_power(8)
pp8
  $is_prime_power
  [1] TRUE
  
  $p
  [1] 2
  
  $n
  [1] 3

pp10 <- is_prime_power(10)
pp10
  $is_prime_power
  [1] FALSE
  
  $p
  [1] NA
  
  $n
  [1] NA
```

The PP construction requires \\q = p^n\\ for a prime \\p\\ and an
integer \\n \geq 1\\.

## 3 Creating `GF(q)`

``` r
gf7 <- gf(7)
gf7
  Galois Field Prime field GF(7) = Z/7Z 
    Order q = 7 elements
    Characteristic p = 7 
    Extension degree n = 1

gf8 <- gf(8)
gf8
  Galois Field GF(2^3) with irreducible x^3 + x + 1 
    Order q = 8 elements
    Characteristic p = 2 
    Extension degree n = 3
```

For `GF(7)`, the package uses modular arithmetic. For `GF(8)`, the
package uses an extension-field construction compatible with the
package’s finite-field implementation choices for PP.

## 4 Supported Fields

The exported [`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md)
constructor supports:

- all prime fields \\GF(p)\\
- binary extensions through `GF(256)` for the implemented irreducible
  polynomials
- selected nonbinary extensions: `GF(9)`, `GF(25)`, `GF(27)`, and
  `GF(49)`

## 5 Arithmetic Examples

``` r
gf_add(3, 5, gf7)
  [1] 1
gf_mult(3, 5, gf7)
  [1] 1

gf_add(5, 3, gf8)
  [1] 6
gf_mult(2, 4, gf8)
  [1] 3
gf_mult(5, 3, gf8)
  [1] 4
gf_pow(3, 7, gf8)
  [1] 1
```

These operations are used inside
[`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)
when the construction evaluates polynomials over `GF(q)`.

## 6 Lookup Tables

The field objects also expose lookup tables. Printing a small slice is
useful for understanding the representation.

``` r
gf8$add[1:4, 1:4]
       [,1] [,2] [,3] [,4]
  [1,]    0    1    2    3
  [2,]    1    0    3    2
  [3,]    2    3    0    1
  [4,]    3    2    1    0
gf8$mult[1:4, 1:4]
       [,1] [,2] [,3] [,4]
  [1,]    0    0    0    0
  [2,]    0    1    2    3
  [3,]    0    2    4    6
  [4,]    0    3    6    5
```

## 7 Prime Field And Extension Field Contrast

``` r
gf7$add[1:4, 1:4]
       [,1] [,2] [,3] [,4]
  [1,]    0    1    2    3
  [2,]    1    2    3    4
  [3,]    2    3    4    5
  [4,]    3    4    5    6
gf7$mult[1:4, 1:4]
       [,1] [,2] [,3] [,4]
  [1,]    0    0    0    0
  [2,]    0    1    2    3
  [3,]    0    2    4    6
  [4,]    0    3    6    2

gf9 <- gf(9)
gf9
  Galois Field GF(3^2) with irreducible x^2 + 2*x + 2 
    Order q = 9 elements
    Characteristic p = 3 
    Extension degree n = 2
gf9$add[1:4, 1:4]
       [,1] [,2] [,3] [,4]
  [1,]    0    1    2    3
  [2,]    1    2    0    4
  [3,]    2    0    1    5
  [4,]    3    4    5    6
gf9$mult[1:4, 1:4]
       [,1] [,2] [,3] [,4]
  [1,]    0    0    0    0
  [2,]    0    1    2    3
  [3,]    0    2    1    6
  [4,]    0    3    6    4
```

## 8 Relationship To Polynomial Pools

In the PP construction, the matrix \\\mathbf{M}\\ is built by evaluating
field-valued polynomial expressions and then converting those
evaluations into binary pool assignments. The field arithmetic in this
vignette is therefore the algebraic foundation for
[`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md).

## 9 References

- Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
  Cambridge University Press.
- Tan, Y. H. I. (2020). Pooling matrix designs for group testing. *SIAM
  Undergraduate Research Online*, 13, 1-21.

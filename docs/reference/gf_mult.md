# Field Multiplication in GF(q)

Computes the product \\a \otimes b\\ in a Galois Field using precomputed
lookup tables for O(1) performance.

## Usage

``` r
gf_mult(a, b, gf)
```

## Arguments

- a:

  Integer scalar representing the first field element. It must lie in
  \\\\0,1,\ldots,q-1\\\\, where `q = gf$q`.

- b:

  Integer scalar representing the second field element. It must lie in
  \\\\0,1,\ldots,q-1\\\\, where `q = gf$q`.

- gf:

  A `"galois_field"` object created by
  [`gf`](https://myaseen208.github.io/ppgt/reference/gf.md). It contains
  the precomputed \\q \times q\\ multiplication table used to evaluate
  \\a \otimes b\\.

## Value

An integer in \\\\0, 1, \ldots, q-1\\\\ representing \\a \otimes b\\.

## Multiplication Rules by Field Type

### Prime Fields GF(p)

\$\$a \otimes b = (a \times b) \mod p\$\$

### Extension Fields GF(p^n)

Polynomial multiplication modulo the irreducible polynomial.

## Example in GF(8)

Compute \\5 \otimes 3\\ in GF(8):

\\5 = x^2 + 1\\ and \\3 = x + 1\\

\$\$(x^2 + 1)(x + 1) = x^3 + x^2 + x + 1\$\$

Reduce mod \\x^3 + x + 1\\: \$\$x^3 + x^2 + x + 1 \equiv x^2 \pmod{x^3 +
x + 1}\$\$ (since \\x^3 \equiv x + 1\\)

So \\5 \otimes 3 = 4\\.

## Key Property

Every non-zero element has a multiplicative inverse: \$\$\forall a \neq
0, \exists a^{-1}: a \otimes a^{-1} = 1\$\$

## References

Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
Cambridge University Press.

Roman, S. (2006). *Field Theory* (2nd ed.). Springer.

## See also

[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_pow`](https://myaseen208.github.io/ppgt/reference/gf_pow.md)

Other galois_field:
[`galois_field_arithmetic`](https://myaseen208.github.io/ppgt/reference/galois_field_arithmetic.md),
[`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md),
[`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
[`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`print.galois_field()`](https://myaseen208.github.io/ppgt/reference/print.galois_field.md)

## Examples

``` r
# GF(8) multiplication
gf8 <- gf(8)
gf_mult(5, 3, gf8)  # (x^2+1)(x+1) = x^2 = 4 (mod x^3+x+1)
#> [1] 4
gf_mult(2, 4, gf8)  # x * x^2 = x^3 = x+1 = 3 (mod x^3+x+1)
#> [1] 3
gf_mult(0, 5, gf8)  # 0 x anything = 0
#> [1] 0
gf_mult(1, 5, gf8)  # 1 is identity: 1 x 5 = 5
#> [1] 5

# GF(7) multiplication (modular)
gf7 <- gf(7)
gf_mult(3, 5, gf7)  # (3 x 5) mod 7 = 15 mod 7 = 1
#> [1] 1
gf_mult(2, 4, gf7)  # (2 x 4) mod 7 = 8 mod 7 = 1
#> [1] 1

# Note: 3 and 5 are inverses in GF(7)!

# Commutative property: a x b = b x a
gf_mult(3, 5, gf8) == gf_mult(5, 3, gf8)  # TRUE
#> [1] TRUE
```

# Field Addition in GF(q)

Computes the sum \\a \oplus b\\ in a Galois Field using precomputed
lookup tables for O(1) performance.

## Usage

``` r
gf_add(a, b, gf)
```

## Arguments

- a:

  Integer scalar representing the first field element. It must lie in
  \\\\0,1,\ldots,q-1\\\\, where \\q = gf\$q\\.

- b:

  Integer scalar representing the second field element. It must lie in
  \\\\0,1,\ldots,q-1\\\\, where \\q = gf\$q\\.

- gf:

  A `"galois_field"` object created by
  [`gf`](https://myaseen208.github.io/ppgt/reference/gf.md). It contains
  the precomputed \\q \times q\\ addition table used to evaluate \\a
  \oplus b\\.

## Value

An integer in \\\\0, 1, \ldots, q-1\\\\ representing \\a \oplus b\\.

## Addition Rules by Field Type

### Prime Fields GF(p)

\$\$a \oplus b = (a + b) \mod p\$\$

### Binary Extension Fields GF(2^n)

\$\$a \oplus b = a \text{ XOR } b\$\$

This is coefficient-wise addition mod 2, which is XOR.

### General Extension Fields GF(p^n), p \> 2

Coefficient-wise addition modulo p.

## Key Property

In any \\GF(2^n)\\: \\a \oplus a = 0\\ (every element is its own
additive inverse)

## References

Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
Cambridge University Press.

Roman, S. (2006). *Field Theory* (2nd ed.). Springer.

## See also

[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_mult`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`gf_pow`](https://myaseen208.github.io/ppgt/reference/gf_pow.md)

Other galois_field:
[`galois_field_arithmetic`](https://myaseen208.github.io/ppgt/reference/galois_field_arithmetic.md),
[`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md),
[`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
[`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`print.galois_field()`](https://myaseen208.github.io/ppgt/reference/print.galois_field.md)

## Examples

``` r
# GF(8) - Binary extension field (XOR addition)
gf8 <- gf(8)
gf_add(5, 3, gf8)  # 101 XOR 011 = 110 = 6
#> [1] 6
gf_add(7, 7, gf8)  # 111 XOR 111 = 000 = 0 (a + a = 0)
#> [1] 0
gf_add(0, 5, gf8)  # 0 is identity: 0 + 5 = 5
#> [1] 5

# GF(7) - Prime field (modular addition)
gf7 <- gf(7)
gf_add(3, 5, gf7)  # (3 + 5) mod 7 = 1
#> [1] 1
gf_add(6, 6, gf7)  # (6 + 6) mod 7 = 5
#> [1] 5

# GF(9) - Extension field GF(3^2)
gf9 <- gf(9)
gf_add(5, 7, gf9)  # Addition in GF(3^2)
#> [1] 0

# Commutative property: a + b = b + a
gf_add(3, 5, gf8) == gf_add(5, 3, gf8)  # TRUE
#> [1] TRUE
```

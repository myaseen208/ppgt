# Field Exponentiation in GF(q)

Computes \\a^n\\ in a Galois Field using the efficient
square-and-multiply algorithm.

## Usage

``` r
gf_pow(a, n, gf)
```

## Arguments

- a:

  Integer scalar representing the base field element. It must lie in
  \\\\0,1,\ldots,q-1\\\\, where `q = gf$q`.

- n:

  Non-negative integer scalar giving the exponent. The implementation
  applies repeated squaring to compute \\a^n\\.

- gf:

  A `"galois_field"` object created by
  [`gf`](https://myaseen208.github.io/ppgt/reference/gf.md). Its
  multiplication table determines the field product used in the power
  recursion.

## Value

An integer in \\\\0, 1, \ldots, q-1\\\\ representing \\a^n\\.

## Algorithm

Uses **square-and-multiply** (binary exponentiation):

1.  Write \\n\\ in binary: \\n = \sum b_i 2^i\\

2.  Compute \\a^{2^i}\\ by repeated squaring

3.  Multiply together the powers where \\b_i = 1\\

Time complexity: \\O(\log n)\\ multiplications.

## Fermat's Little Theorem

For any \\a \neq 0\\ in \\GF(q)\\: \$\$a^{q-1} = 1\$\$

This is fundamental to finite field theory and implies: \$\$a^{-1} =
a^{q-2}\$\$

## Example in GF(8)

Since \\\|GF(8)^\*\| = 7\\: \$\$a^7 = 1 \text{ for all } a \neq 0\$\$

## References

Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
Cambridge University Press.

Menezes, A. J., van Oorschot, P. C., & Vanstone, S. A. (1996). *Handbook
of Applied Cryptography*. CRC Press.

## See also

[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult`](https://myaseen208.github.io/ppgt/reference/gf_mult.md)

Other galois_field:
[`galois_field_arithmetic`](https://myaseen208.github.io/ppgt/reference/galois_field_arithmetic.md),
[`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
[`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`print.galois_field()`](https://myaseen208.github.io/ppgt/reference/print.galois_field.md)

## Examples

``` r
# GF(8) exponentiation
gf8 <- gf(8)
gf_pow(2, 3, gf8)  # x^3 = x + 1 = 3 (mod x^3+x+1)
#> [1] 3
gf_pow(2, 7, gf8)  # 2^7 = 1 (Fermat's Little Theorem)
#> [1] 1

# Verify Fermat's theorem for all non-zero elements
all(sapply(1:7, function(a) gf_pow(a, 7, gf8)) == 1)  # TRUE
#> [1] TRUE

# GF(7) exponentiation
gf7 <- gf(7)
gf_pow(3, 6, gf7)  # 3^6 = 729 mod 7 = 1 (Fermat)
#> [1] 1
gf_pow(3, 3, gf7)  # 3^3 = 27 mod 7 = 6
#> [1] 6

# Special cases
gf_pow(5, 0, gf8)  # a^0 = 1
#> [1] 1
gf_pow(5, 1, gf8)  # a^1 = a = 5
#> [1] 5
gf_pow(0, 5, gf8)  # 0^n = 0
#> [1] 0

# Compute multiplicative inverse using Fermat: a^(-1) = a^(q-2)
a <- 5
a_inv <- gf_pow(a, 6, gf8)  # 5^(-1) in GF(8)
gf_mult(a, a_inv, gf8)      # Should be 1
#> [1] 1
```

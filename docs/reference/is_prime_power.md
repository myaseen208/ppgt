# Check if Integer is a Prime Power

Tests whether \\q = p^n\\ for some prime \\p\\ and positive integer \\n
\geq 1\\. This is the necessary and sufficient condition for the
existence of a Galois Field of order \\q\\.

## Usage

``` r
is_prime_power(q)
```

## Arguments

- q:

  Integer scalar to test as a candidate field order. It plays the role
  of the finite-field size \\q\\ and must satisfy \\q \ge 1\\.

## Value

A named list with three components:

- is_prime_power:

  Logical; `TRUE` if `q` is a prime power

- p:

  Integer; the prime base (characteristic) if prime power, else `NA`

- n:

  Integer; the exponent (extension degree) if prime power, else `NA`

## Mathematical Definition

A **prime power** is any positive integer expressible as: \$\$q =
p^n\$\$ where \\p\\ is prime and \\n \geq 1\\.

## Valid Prime Powers (Galois Field Orders)

The following are valid orders for Galois Fields: \$\$2, 3, 4, 5, 7, 8,
9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 53,
59, 61, 64, \ldots\$\$

## Invalid Values (NOT Prime Powers)

These **cannot** be Galois Field orders: \$\$6, 10, 12, 14, 15, 18, 20,
21, 22, 24, 26, 28, 30, \ldots\$\$

## Why This Matters

The Polynomial Pools algorithm requires \\q\\ to be a prime power
because:

1.  **Field existence**: Galois Fields only exist for prime power orders

2.  **Division**: We need every non-zero element to have an inverse

3.  **Polynomial theory**: Root-counting theorems require field
    structure

For \\q = 6 = 2 \times 3\\, we have \\2 \times 3 = 0 \pmod{6}\\ (zero
divisors), so neither 2 nor 3 has a multiplicative inverse. This
violates the field axioms.

## References

Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
Cambridge University Press.

Wan, Z.-X. (2003). *Lectures on Finite Fields and Galois Rings*. World
Scientific.

## See also

[`is_prime`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md)

Other galois_field:
[`galois_field_arithmetic`](https://myaseen208.github.io/ppgt/reference/galois_field_arithmetic.md),
[`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md),
[`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
[`print.galois_field()`](https://myaseen208.github.io/ppgt/reference/print.galois_field.md)

## Examples

``` r
# Prime powers with decomposition
is_prime_power(8)
#> $is_prime_power
#> [1] TRUE
#> 
#> $p
#> [1] 2
#> 
#> $n
#> [1] 3
#> 
is_prime_power(9)
#> $is_prime_power
#> [1] TRUE
#> 
#> $p
#> [1] 3
#> 
#> $n
#> [1] 2
#> 
is_prime_power(49)
#> $is_prime_power
#> [1] TRUE
#> 
#> $p
#> [1] 7
#> 
#> $n
#> [1] 2
#> 
is_prime_power(64)
#> $is_prime_power
#> [1] TRUE
#> 
#> $p
#> [1] 2
#> 
#> $n
#> [1] 6
#> 
is_prime_power(81)
#> $is_prime_power
#> [1] TRUE
#> 
#> $p
#> [1] 3
#> 
#> $n
#> [1] 4
#> 

# Primes are prime powers with n = 1
is_prime_power(7)
#> $is_prime_power
#> [1] TRUE
#> 
#> $p
#> [1] 7
#> 
#> $n
#> [1] 1
#> 
is_prime_power(13)
#> $is_prime_power
#> [1] TRUE
#> 
#> $p
#> [1] 13
#> 
#> $n
#> [1] 1
#> 

# NOT prime powers (cannot be Galois Field orders)
is_prime_power(6)
#> $is_prime_power
#> [1] FALSE
#> 
#> $p
#> [1] NA
#> 
#> $n
#> [1] NA
#> 
is_prime_power(10)
#> $is_prime_power
#> [1] FALSE
#> 
#> $p
#> [1] NA
#> 
#> $n
#> [1] NA
#> 
is_prime_power(12)
#> $is_prime_power
#> [1] FALSE
#> 
#> $p
#> [1] NA
#> 
#> $n
#> [1] NA
#> 
is_prime_power(15)
#> $is_prime_power
#> [1] FALSE
#> 
#> $p
#> [1] NA
#> 
#> $n
#> [1] NA
#> 

# Extract components
result <- is_prime_power(27)
result
#> $is_prime_power
#> [1] TRUE
#> 
#> $p
#> [1] 3
#> 
#> $n
#> [1] 3
#> 
result$is_prime_power
#> [1] TRUE
result$p
#> [1] 3
result$n
#> [1] 3
```

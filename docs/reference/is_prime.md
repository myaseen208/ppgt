# Check if Integer is Prime

Tests whether a positive integer is a prime number using optimized trial
division.

## Usage

``` r
is_prime(x)
```

## Arguments

- x:

  A positive integer to test.

## Value

A logical value:

- TRUE:

  if `x` is a prime number

- FALSE:

  if `x` is not prime (composite, 0, 1, or negative)

## Definition

A **prime number** \\p \> 1\\ has exactly two positive divisors: 1 and
itself.

## Algorithm

Uses trial division checking divisibility only up to \\\sqrt{x}\\,
skipping even numbers after checking 2.

## Role in Galois Fields

Prime numbers are the building blocks of Galois Fields. Every finite
field has order \\q = p^n\\ where \\p\\ is prime.

## References

Hardy, G. H., & Wright, E. M. (1979). *An Introduction to the Theory of
Numbers* (5th ed.). Oxford University Press.

Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
Cambridge University Press.

## See also

[`is_prime_power`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md)

Other galois_field:
[`galois_field_arithmetic`](https://myaseen208.github.io/ppgt/reference/galois_field_arithmetic.md),
[`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md),
[`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`print.galois_field()`](https://myaseen208.github.io/ppgt/reference/print.galois_field.md)

## Examples

``` r
# Prime numbers
is_prime(2)
#> [1] TRUE
is_prime(7)
#> [1] TRUE
is_prime(97)
#> [1] TRUE
is_prime(104729)
#> [1] TRUE

# Non-prime numbers
is_prime(1)
#> [1] FALSE
is_prime(4)
#> [1] FALSE
is_prime(9)
#> [1] FALSE
is_prime(100)
#> [1] FALSE

# First 10 primes
primes <- sapply(2:30, is_prime)
primes
#>  [1]  TRUE  TRUE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE  TRUE
#> [13] FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE
#> [25] FALSE FALSE FALSE  TRUE FALSE
which(primes) + 1
#>  [1]  2  3  5  7 11 13 17 19 23 29
```

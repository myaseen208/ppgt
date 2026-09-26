# Print a Galois Field Object

Prints a compact summary of a `"galois_field"` object created by
[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md).

## Usage

``` r
# S3 method for class 'galois_field'
print(x, ...)
```

## Arguments

- x:

  An object of class `"galois_field"` returned by
  [`gf`](https://myaseen208.github.io/ppgt/reference/gf.md).

- ...:

  Unused arguments passed through the generic
  [`print()`](https://rdrr.io/r/base/print.html) interface.

## Value

The input object `x`, invisibly.

## Details

Let \\x = \mathrm{GF}(q)\\ be the field object returned by
[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md), where \$\$ q
= p^n \$\$ with prime characteristic \\p\\ and extension degree \\n \ge
1\\. The printed summary reports:

- the field description,

- the order \\q\\,

- the characteristic \\p\\, and

- the extension degree \\n\\.

This is a display method only; it does not modify the object or its
lookup tables.

## References

Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
Cambridge University Press.

## See also

[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`gf_pow`](https://myaseen208.github.io/ppgt/reference/gf_pow.md)

Other galois_field:
[`galois_field_arithmetic`](https://myaseen208.github.io/ppgt/reference/galois_field_arithmetic.md),
[`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md),
[`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
[`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md)

## Examples

``` r
gf8 <- gf(8)
print(gf8)
#> Galois Field GF(2^3) with irreducible x^3 + x + 1 
#>   Order q = 8 elements
#>   Characteristic p = 2 
#>   Extension degree n = 3 

gf7 <- gf(7)
print(gf7)
#> Galois Field Prime field GF(7) = Z/7Z 
#>   Order q = 7 elements
#>   Characteristic p = 7 
#>   Extension degree n = 1 
```

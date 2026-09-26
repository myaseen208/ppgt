# Create a Galois Field GF(q)

Constructs a Galois Field (finite field) of order \\q\\ with precomputed
addition and multiplication lookup tables for O(1) field operations.

## Usage

``` r
gf(q)
```

## Arguments

- q:

  Prime power (the order of the field). Must be \\p^n\\ for some prime
  \\p\\ and positive integer \\n\\.

## Value

An object of class `"galois_field"` containing:

- q:

  Integer; order of the field (number of elements)

- p:

  Integer; prime characteristic

- n:

  Integer; extension degree (\\q = p^n\\)

- add:

  Integer matrix (\\q \times q\\); addition lookup table.
  `add[i+1, j+1]` gives \\i \oplus j\\.

- mult:

  Integer matrix (\\q \times q\\); multiplication lookup table.
  `mult[i+1, j+1]` gives \\i \otimes j\\.

- description:

  Character; human-readable description of the field

## Mathematical Construction

### Prime Fields (n = 1)

When \\q = p\\ is prime, the field is: \$\$GF(p) \cong \mathbb{Z}\_p =
\\0, 1, 2, \ldots, p-1\\\$\$

with modular arithmetic: \$\$a \oplus b = (a + b) \mod p\$\$ \$\$a
\otimes b = (a \times b) \mod p\$\$

**Example: GF(7)** \$\$3 \oplus 5 = 8 \mod 7 = 1\$\$ \$\$3 \otimes 5 =
15 \mod 7 = 1\$\$

### Extension Fields (n \> 1)

When \\q = p^n\\ with \\n \> 1\\, we construct: \$\$GF(p^n) \cong
GF(p)\[x\] / (f(x))\$\$

where \\f(x)\\ is an **irreducible polynomial** of degree \\n\\ over
\\GF(p)\\.

Elements are polynomials of degree \\\< n\\ with coefficients in
\\GF(p)\\: \$\$a\_{n-1}x^{n-1} + a\_{n-2}x^{n-2} + \cdots + a_1 x +
a_0\$\$

These are encoded as integers: \$\$\sum\_{i=0}^{n-1} a_i \cdot p^i \in
\\0, 1, \ldots, q-1\\\$\$

### GF(8) = GF(2^3) Construction (Used in P-BEST)

Irreducible polynomial: \\f(x) = x^3 + x + 1\\

|         |        |                 |
|---------|--------|-----------------|
| Integer | Binary | Polynomial      |
| 0       | 000    | 0               |
| 1       | 001    | 1               |
| 2       | 010    | \\x\\           |
| 3       | 011    | \\x + 1\\       |
| 4       | 100    | \\x^2\\         |
| 5       | 101    | \\x^2 + 1\\     |
| 6       | 110    | \\x^2 + x\\     |
| 7       | 111    | \\x^2 + x + 1\\ |

**Addition in GF(2^n)**: Coefficient-wise XOR \$\$5 \oplus 3 =
(x^2 + 1) + (x + 1) = x^2 + x = 6\$\$ \$\$101_2 \oplus 011_2 = 110_2\$\$

**Multiplication in GF(2^n)**: Polynomial multiplication mod \\f(x)\\
\$\$2 \otimes 4 = x \cdot x^2 = x^3 \equiv x + 1 = 3 \pmod{x^3 + x +
1}\$\$

### GF(9) = GF(3^2) Construction

Irreducible polynomial: \\f(x) = x^2 + 2x + 2\\ over \\GF(3)\\ (Conway
polynomial)

Elements are \\a_1 x + a_0\\ with \\a_i \in \\0, 1, 2\\\\, encoded as
\\3a_1 + a_0\\.

## Supported Fields and Conway Polynomials

All extension fields use the **Conway polynomial** from Frank Luebeck's
database, consistent with the choice made by GaloisFields.jl.

|         |         |                               |
|---------|---------|-------------------------------|
| Field   | Order   | Conway polynomial             |
| GF(4)   | \\2^2\\ | \\x^2 + x + 1\\               |
| GF(8)   | \\2^3\\ | \\x^3 + x + 1\\               |
| GF(16)  | \\2^4\\ | \\x^4 + x + 1\\               |
| GF(32)  | \\2^5\\ | \\x^5 + x^2 + 1\\             |
| GF(64)  | \\2^6\\ | \\x^6 + x^4 + x^3 + x + 1\\   |
| GF(128) | \\2^7\\ | \\x^7 + x + 1\\               |
| GF(256) | \\2^8\\ | \\x^8 + x^4 + x^3 + x^2 + 1\\ |
| GF(9)   | \\3^2\\ | \\x^2 + 2x + 2\\              |
| GF(27)  | \\3^3\\ | \\x^3 + 2x + 1\\              |
| GF(25)  | \\5^2\\ | \\x^2 + 4x + 2\\              |
| GF(49)  | \\7^2\\ | \\x^2 + 6x + 3\\              |

- **Prime fields**: GF(2), GF(3), GF(5), GF(7), GF(11), GF(13), ...

- **Binary extension fields**: GF(4), GF(8), GF(16), GF(32), GF(64),
  GF(128), GF(256)

- **Other extension fields**: GF(9), GF(25), GF(27), GF(49)

## References

Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
Cambridge University Press.

## See also

[`is_prime_power`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`gf_add`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`gf_pow`](https://myaseen208.github.io/ppgt/reference/gf_pow.md),
[`print.galois_field`](https://myaseen208.github.io/ppgt/reference/print.galois_field.md)

Other galois_field:
[`galois_field_arithmetic`](https://myaseen208.github.io/ppgt/reference/galois_field_arithmetic.md),
[`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md),
[`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
[`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`print.galois_field()`](https://myaseen208.github.io/ppgt/reference/print.galois_field.md)

## Examples

``` r
# GF(8) - Used in P-BEST design
gf8 <- gf(8)
print(gf8)
#> Galois Field GF(2^3) with irreducible x^3 + x + 1 
#>   Order q = 8 elements
#>   Characteristic p = 2 
#>   Extension degree n = 3 

# Verify addition table (XOR)
gf8$add[6, 4]  # Element 5 + 3 = 6 (index is value + 1)
#> [1] 6

# Verify multiplication
gf8$mult[3, 5]  # Element 2 x 4 = 3
#> [1] 3

# GF(9) = GF(3^2) for larger designs
gf9 <- gf(9)
print(gf9)
#> Galois Field GF(3^2) with irreducible x^2 + 2*x + 2 
#>   Order q = 9 elements
#>   Characteristic p = 3 
#>   Extension degree n = 2 

# Prime field GF(7)
gf7 <- gf(7)
gf7$add[4, 6]   # 3 + 5 = 8 mod 7 = 1
#> [1] 1
gf7$mult[4, 6]  # 3 x 5 = 15 mod 7 = 1
#> [1] 1

# Error: 6 is not a prime power
if (FALSE) { # \dontrun{
gf6 <- gf(6)  # Error!
} # }
```

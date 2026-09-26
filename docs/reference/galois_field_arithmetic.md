# Galois Field (Finite Field) Arithmetic

Complete implementation of Galois Field (finite field) arithmetic
required for the Polynomial Pools group testing algorithm.

## Mathematical Background

A **Galois Field** (or **finite field**) \\GF(q)\\ is an algebraic
structure containing exactly \\q\\ elements with well-defined addition,
subtraction, multiplication, and division operations satisfying the
field axioms.

### Fundamental Theorem

Finite fields exist if and only if the order \\q\\ is a **prime power**:
\$\$q = p^n\$\$ where \\p\\ is a prime number (the **characteristic**)
and \\n \geq 1\\ is a positive integer (the **extension degree**).

### Why Galois Fields for Group Testing?

The Polynomial Pools algorithm constructs pooling matrices by evaluating
polynomials over \\GF(q)\\. The field structure ensures:

1.  **Closure**: All operations remain within the field

2.  **Invertibility**: Every non-zero element has a multiplicative
    inverse

3.  **Polynomial roots**: Polynomials of degree \\d\\ have at most \\d\\
    roots

These properties guarantee the **intersection bound**: any two samples
share at most \\d-1\\ pools, which is essential for the detection
guarantees.

### Relationship to Polynomial Pools and P-BEST

- **Galois Field**: Provides the arithmetic foundation

- **Polynomial Pools**: Uses GF(q) to construct pooling matrices

- **P-BEST**: A specific PP design using GF(8) = GF(2^3)

## References

Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
Cambridge University Press.

## See also

[`pp_matrix`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)
for using Galois Fields in pooling matrix construction,
[`pp_design`](https://myaseen208.github.io/ppgt/reference/pp_design.md)
for exploring valid configurations

Other galois_field:
[`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
[`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
[`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md),
[`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
[`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`print.galois_field()`](https://myaseen208.github.io/ppgt/reference/print.galois_field.md)

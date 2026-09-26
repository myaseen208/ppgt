#' @title Galois Field (Finite Field) Arithmetic
#'
#' @description
#' Complete implementation of Galois Field (finite field) arithmetic required
#' for the Polynomial Pools group testing algorithm.
#'
#' @details
#' # Mathematical Background
#'
#' A **Galois Field** (or **finite field**) \eqn{GF(q)} is an algebraic structure
#' containing exactly \eqn{q} elements with well-defined addition, subtraction,
#' multiplication, and division operations satisfying the field axioms.
#'
#' ## Fundamental Theorem
#'
#' Finite fields exist if and only if the order \eqn{q} is a **prime power**:
#' \deqn{q = p^n}
#' where \eqn{p} is a prime number (the **characteristic**) and \eqn{n \geq 1}
#' is a positive integer (the **extension degree**).
#'
#' ## Why Galois Fields for Group Testing?
#'
#' The Polynomial Pools algorithm constructs pooling matrices by evaluating
#' polynomials over \eqn{GF(q)}. The field structure ensures:
#'
#' 1. **Closure**: All operations remain within the field
#' 2. **Invertibility**: Every non-zero element has a multiplicative inverse
#' 3. **Polynomial roots**: Polynomials of degree \eqn{d} have at most \eqn{d} roots
#'
#' These properties guarantee the **intersection bound**: any two samples share

#' at most \eqn{d-1} pools, which is essential for the detection guarantees.
#'
#' ## Relationship to Polynomial Pools and P-BEST
#'
#' - **Galois Field**: Provides the arithmetic foundation
#' - **Polynomial Pools**: Uses GF(q) to construct pooling matrices
#' - **P-BEST**: A specific PP design using GF(8) = GF(2^3)
#'
#' @name galois_field_arithmetic
#' @family galois_field
#'
#' @references
#' Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
#' Cambridge University Press.
#'
#' @seealso
#' \code{\link{pp_matrix}} for using Galois Fields in pooling matrix construction,
#' \code{\link{pp_design}} for exploring valid configurations
NULL


#' @title Check if Integer is Prime
#'
#' @description
#' Tests whether a positive integer is a prime number using optimized
#' trial division.
#'
#' @details
#' # Definition
#'
#' A **prime number** \eqn{p > 1} has exactly two positive divisors: 1 and itself.
#'
#' # Algorithm
#'
#' Uses trial division checking divisibility only up to \eqn{\sqrt{x}},
#' skipping even numbers after checking 2.
#'
#' # Role in Galois Fields
#'
#' Prime numbers are the building blocks of Galois Fields. Every finite field
#' has order \eqn{q = p^n} where \eqn{p} is prime.
#'
#' @param x A positive integer to test.
#'
#' @return A logical value:
#' \describe{
#'   \item{TRUE}{if \code{x} is a prime number}
#'   \item{FALSE}{if \code{x} is not prime (composite, 0, 1, or negative)}
#' }
#'
#' @examples
#' # Prime numbers
#' is_prime(2)
#' is_prime(7)
#' is_prime(97)
#' is_prime(104729)
#'
#' # Non-prime numbers
#' is_prime(1)
#' is_prime(4)
#' is_prime(9)
#' is_prime(100)
#'
#' # First 10 primes
#' primes <- sapply(2:30, is_prime)
#' primes
#' which(primes) + 1
#'
#' @references
#' Hardy, G. H., & Wright, E. M. (1979). \emph{An Introduction to the Theory of
#' Numbers} (5th ed.). Oxford University Press.
#'
#' Lidl, R., & Niederreiter, H. (1997). \emph{Finite Fields} (2nd ed.).
#' Cambridge University Press.
#'
#' @seealso
#' \code{\link{is_prime_power}}, \code{\link{gf}}
#'
#' @family galois_field
#' @export
is_prime <- function(x) {
  x <- as.integer(x)
  if (x < 2L) return(FALSE)
  if (x == 2L) return(TRUE)
  if (x %% 2L == 0L) return(FALSE)
  i <- 3L
  while (i * i <= x) {
    if (x %% i == 0L) return(FALSE)
    i <- i + 2L
  }
  TRUE
}


#' @title Check if Integer is a Prime Power
#'
#' @description
#' Tests whether \eqn{q = p^n} for some prime \eqn{p} and positive integer
#' \eqn{n \geq 1}. This is the necessary and sufficient condition for the
#' existence of a Galois Field of order \eqn{q}.
#'
#' @details
#' # Mathematical Definition
#'
#' A **prime power** is any positive integer expressible as:
#' \deqn{q = p^n}
#' where \eqn{p} is prime and \eqn{n \geq 1}.
#'
#' # Valid Prime Powers (Galois Field Orders)
#'
#' The following are valid orders for Galois Fields:
#' \deqn{2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37, 41, 43, 47, 49, 53, 59, 61, 64, \ldots}
#'
#' # Invalid Values (NOT Prime Powers)
#'
#' These **cannot** be Galois Field orders:
#' \deqn{6, 10, 12, 14, 15, 18, 20, 21, 22, 24, 26, 28, 30, \ldots}
#'
#' # Why This Matters
#'
#' The Polynomial Pools algorithm requires \eqn{q} to be a prime power because:
#'
#' 1. **Field existence**: Galois Fields only exist for prime power orders
#' 2. **Division**: We need every non-zero element to have an inverse
#' 3. **Polynomial theory**: Root-counting theorems require field structure
#'
#' For \eqn{q = 6 = 2 \times 3}, we have \eqn{2 \times 3 = 0 \pmod{6}} (zero
#' divisors), so neither 2 nor 3 has a multiplicative inverse. This violates
#' the field axioms.
#'
#' @param q Integer scalar to test as a candidate field order. It plays the role
#'   of the finite-field size \eqn{q} and must satisfy \eqn{q \ge 1}.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{is_prime_power}{Logical; \code{TRUE} if \code{q} is a prime power}
#'   \item{p}{Integer; the prime base (characteristic) if prime power, else \code{NA}}
#'   \item{n}{Integer; the exponent (extension degree) if prime power, else \code{NA}}
#' }
#'
#' @examples
#' # Prime powers with decomposition
#' is_prime_power(8)
#' is_prime_power(9)
#' is_prime_power(49)
#' is_prime_power(64)
#' is_prime_power(81)
#'
#' # Primes are prime powers with n = 1
#' is_prime_power(7)
#' is_prime_power(13)
#'
#' # NOT prime powers (cannot be Galois Field orders)
#' is_prime_power(6)
#' is_prime_power(10)
#' is_prime_power(12)
#' is_prime_power(15)
#'
#' # Extract components
#' result <- is_prime_power(27)
#' result
#' result$is_prime_power
#' result$p
#' result$n
#'
#' @references
#' Lidl, R., & Niederreiter, H. (1997). \emph{Finite Fields} (2nd ed.).
#' Cambridge University Press.
#'
#' Wan, Z.-X. (2003). \emph{Lectures on Finite Fields and Galois Rings}. World
#' Scientific.
#'
#' @seealso
#' \code{\link{is_prime}}, \code{\link{gf}}
#'
#' @family galois_field
#' @export
is_prime_power <- function(q) {
  q <- as.integer(q)
  if (q < 2L) {
    return(list(is_prime_power = FALSE, p = NA_integer_, n = NA_integer_))
  }
  for (p in 2L:q) {
    if (!is_prime(p)) next
    n <- 0L
    temp <- q
    while (temp %% p == 0L) {
      temp <- temp %/% p
      n <- n + 1L
    }
    if (temp == 1L && n >= 1L) {
      return(list(is_prime_power = TRUE, p = p, n = n))
    }
    if (p * p > q) break
  }
  if (is_prime(q)) {
    return(list(is_prime_power = TRUE, p = q, n = 1L))
  }
  list(is_prime_power = FALSE, p = NA_integer_, n = NA_integer_)
}


#' @title Create a Galois Field GF(q)
#'
#' @description
#' Constructs a Galois Field (finite field) of order \eqn{q} with precomputed
#' addition and multiplication lookup tables for O(1) field operations.
#'
#' @details
#' # Mathematical Construction
#'
#' ## Prime Fields (n = 1)
#'
#' When \eqn{q = p} is prime, the field is:
#' \deqn{GF(p) \cong \mathbb{Z}_p = \{0, 1, 2, \ldots, p-1\}}
#'
#' with modular arithmetic:
#' \deqn{a \oplus b = (a + b) \mod p}
#' \deqn{a \otimes b = (a \times b) \mod p}
#'
#' **Example: GF(7)**
#' \deqn{3 \oplus 5 = 8 \mod 7 = 1}
#' \deqn{3 \otimes 5 = 15 \mod 7 = 1}
#'
#' ## Extension Fields (n > 1)
#'
#' When \eqn{q = p^n} with \eqn{n > 1}, we construct:
#' \deqn{GF(p^n) \cong GF(p)[x] / (f(x))}
#'
#' where \eqn{f(x)} is an **irreducible polynomial** of degree \eqn{n} over \eqn{GF(p)}.
#'
#' Elements are polynomials of degree \eqn{< n} with coefficients in \eqn{GF(p)}:
#' \deqn{a_{n-1}x^{n-1} + a_{n-2}x^{n-2} + \cdots + a_1 x + a_0}
#'
#' These are encoded as integers:
#' \deqn{\sum_{i=0}^{n-1} a_i \cdot p^i \in \{0, 1, \ldots, q-1\}}
#'
#' ## GF(8) = GF(2^3) Construction (Used in P-BEST)
#'
#' Irreducible polynomial: \eqn{f(x) = x^3 + x + 1}
#'
#' | Integer | Binary | Polynomial |
#' |---------|--------|------------|
#' | 0 | 000 | 0 |
#' | 1 | 001 | 1 |
#' | 2 | 010 | \eqn{x} |
#' | 3 | 011 | \eqn{x + 1} |
#' | 4 | 100 | \eqn{x^2} |
#' | 5 | 101 | \eqn{x^2 + 1} |
#' | 6 | 110 | \eqn{x^2 + x} |
#' | 7 | 111 | \eqn{x^2 + x + 1} |
#'
#' **Addition in GF(2^n)**: Coefficient-wise XOR
#' \deqn{5 \oplus 3 = (x^2 + 1) + (x + 1) = x^2 + x = 6}
#' \deqn{101_2 \oplus 011_2 = 110_2}
#'
#' **Multiplication in GF(2^n)**: Polynomial multiplication mod \eqn{f(x)}
#' \deqn{2 \otimes 4 = x \cdot x^2 = x^3 \equiv x + 1 = 3 \pmod{x^3 + x + 1}}
#'
#' ## GF(9) = GF(3^2) Construction
#'
#' Irreducible polynomial: \eqn{f(x) = x^2 + 2x + 2} over \eqn{GF(3)} (Conway polynomial)
#'
#' Elements are \eqn{a_1 x + a_0} with \eqn{a_i \in \{0, 1, 2\}}, encoded as \eqn{3a_1 + a_0}.
#'
#' # Supported Fields and Conway Polynomials
#'
#' All extension fields use the **Conway polynomial** from Frank Luebeck's database,
#' consistent with the choice made by GaloisFields.jl.
#'
#' | Field | Order | Conway polynomial |
#' |-------|-------|-------------------|
#' | GF(4)   | \eqn{2^2} | \eqn{x^2 + x + 1} |
#' | GF(8)   | \eqn{2^3} | \eqn{x^3 + x + 1} |
#' | GF(16)  | \eqn{2^4} | \eqn{x^4 + x + 1} |
#' | GF(32)  | \eqn{2^5} | \eqn{x^5 + x^2 + 1} |
#' | GF(64)  | \eqn{2^6} | \eqn{x^6 + x^4 + x^3 + x + 1} |
#' | GF(128) | \eqn{2^7} | \eqn{x^7 + x + 1} |
#' | GF(256) | \eqn{2^8} | \eqn{x^8 + x^4 + x^3 + x^2 + 1} |
#' | GF(9)   | \eqn{3^2} | \eqn{x^2 + 2x + 2} |
#' | GF(27)  | \eqn{3^3} | \eqn{x^3 + 2x + 1} |
#' | GF(25)  | \eqn{5^2} | \eqn{x^2 + 4x + 2} |
#' | GF(49)  | \eqn{7^2} | \eqn{x^2 + 6x + 3} |
#'
#' - **Prime fields**: GF(2), GF(3), GF(5), GF(7), GF(11), GF(13), ...
#' - **Binary extension fields**: GF(4), GF(8), GF(16), GF(32), GF(64), GF(128), GF(256)
#' - **Other extension fields**: GF(9), GF(25), GF(27), GF(49)
#'
#' @param q Prime power (the order of the field). Must be \eqn{p^n} for some
#'   prime \eqn{p} and positive integer \eqn{n}.
#'
#' @return An object of class \code{"galois_field"} containing:
#' \describe{
#'   \item{q}{Integer; order of the field (number of elements)}
#'   \item{p}{Integer; prime characteristic}
#'   \item{n}{Integer; extension degree (\eqn{q = p^n})}
#'   \item{add}{Integer matrix (\eqn{q \times q}); addition lookup table.
#'     \code{add[i+1, j+1]} gives \eqn{i \oplus j}.}
#'   \item{mult}{Integer matrix (\eqn{q \times q}); multiplication lookup table.
#'     \code{mult[i+1, j+1]} gives \eqn{i \otimes j}.}
#'   \item{description}{Character; human-readable description of the field}
#' }
#'
#' @examples
#' # GF(8) - Used in P-BEST design
#' gf8 <- gf(8)
#' print(gf8)
#'
#' # Verify addition table (XOR)
#' gf8$add[6, 4]  # Element 5 + 3 = 6 (index is value + 1)
#'
#' # Verify multiplication
#' gf8$mult[3, 5]  # Element 2 x 4 = 3
#'
#' # GF(9) = GF(3^2) for larger designs
#' gf9 <- gf(9)
#' print(gf9)
#'
#' # Prime field GF(7)
#' gf7 <- gf(7)
#' gf7$add[4, 6]   # 3 + 5 = 8 mod 7 = 1
#' gf7$mult[4, 6]  # 3 x 5 = 15 mod 7 = 1
#'
#' # Error: 6 is not a prime power
#' \dontrun{
#' gf6 <- gf(6)  # Error!
#' }
#'
#' @references
#' Lidl, R., & Niederreiter, H. (1997). *Finite Fields* (2nd ed.).
#' Cambridge University Press.
#'
#' @seealso
#' \code{\link{is_prime_power}}, \code{\link{gf_add}},
#' \code{\link{gf_mult}}, \code{\link{gf_pow}}, \code{\link{print.galois_field}}
#'
#' @family galois_field
#' @export
gf <- function(q) {

  pp <- is_prime_power(q)
  if (!pp$is_prime_power) {
    stop("q must be a prime power (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, ...).\n",
         "Got q = ", q, " which is NOT a prime power.\n",
         "Hint: ", q, " = ", paste(factorize(q), collapse = " x "))
  }

  p <- pp$p
  n <- pp$n
  q <- as.integer(q)

  add_table <- matrix(0L, nrow = q, ncol = q)
  mult_table <- matrix(0L, nrow = q, ncol = q)

  if (n == 1L) {
    # Prime field GF(p): modular arithmetic
    for (i in 0L:(q - 1L)) {
      for (j in 0L:(q - 1L)) {
        add_table[i + 1L, j + 1L] <- (i + j) %% p
        mult_table[i + 1L, j + 1L] <- (i * j) %% p
      }
    }
    irr_desc <- paste0("Prime field GF(", p, ") = Z/", p, "Z")

  } else if (p == 2L) {
    # Binary extension field GF(2^n)
    irr_polys <- list(
      "2" = list(poly = 7L, name = "x^2 + x + 1"),
      "3" = list(poly = 11L, name = "x^3 + x + 1"),
      "4" = list(poly = 19L, name = "x^4 + x + 1"),
      "5" = list(poly = 37L, name = "x^5 + x^2 + 1"),
      "6" = list(poly = 91L, name = "x^6 + x^4 + x^3 + x + 1"),
      "7" = list(poly = 131L, name = "x^7 + x + 1"),
      "8" = list(poly = 285L, name = "x^8 + x^4 + x^3 + x^2 + 1")
    )

    if (!(as.character(n) %in% names(irr_polys))) {
      stop("GF(2^", n, ") not implemented. Supported: n = 2, 3, 4, 5, 6, 7, 8")
    }

    irr_info <- irr_polys[[as.character(n)]]
    irr <- irr_info$poly
    irr_desc <- paste0("GF(2^", n, ") with irreducible ", irr_info$name)

    for (i in 0L:(q - 1L)) {
      for (j in 0L:(q - 1L)) {
        add_table[i + 1L, j + 1L] <- bitwXor(i, j)
        mult_table[i + 1L, j + 1L] <- gf2n_mult_internal(i, j, n, irr)
      }
    }

  } else {
    # General extension field GF(p^n) with p > 2
    irr_polys <- list(
      "3_2" = list(coeffs = c(2L, 2L, 1L), name = "x^2 + 2*x + 2"),
      "3_3" = list(coeffs = c(1L, 2L, 0L, 1L), name = "x^3 + 2*x + 1"),
      "5_2" = list(coeffs = c(2L, 4L, 1L), name = "x^2 + 4*x + 2"),
      "7_2" = list(coeffs = c(3L, 6L, 1L), name = "x^2 + 6*x + 3")
    )

    key <- paste0(p, "_", n)
    if (!(key %in% names(irr_polys))) {
      stop("GF(", p, "^", n, ") not implemented.\n",
           "Supported extension fields: GF(3^2)=GF(9), GF(3^3)=GF(27), ",
           "GF(5^2)=GF(25), GF(7^2)=GF(49)")
    }

    irr_info <- irr_polys[[key]]
    irr <- irr_info$coeffs
    irr_desc <- paste0("GF(", p, "^", n, ") with irreducible ", irr_info$name)

    for (i in 0L:(q - 1L)) {
      for (j in 0L:(q - 1L)) {
        add_table[i + 1L, j + 1L] <- gfpn_add_internal(i, j, p, n)
        mult_table[i + 1L, j + 1L] <- gfpn_mult_internal(i, j, p, n, irr)
      }
    }
  }

  structure(
    list(
      q = q, p = p, n = n,
      add = add_table,
      mult = mult_table,
      description = irr_desc
    ),
    class = "galois_field"
  )
}


#' @title Field Addition in GF(q)
#'
#' @description
#' Computes the sum \eqn{a \oplus b} in a Galois Field using precomputed
#' lookup tables for O(1) performance.
#'
#' @details
#' # Addition Rules by Field Type
#'
#' ## Prime Fields GF(p)
#' \deqn{a \oplus b = (a + b) \mod p}
#'
#' ## Binary Extension Fields GF(2^n)
#' \deqn{a \oplus b = a \text{ XOR } b}
#'
#' This is coefficient-wise addition mod 2, which is XOR.
#'
#' ## General Extension Fields GF(p^n), p > 2
#' Coefficient-wise addition modulo p.
#'
#' # Key Property
#'
#' In any \eqn{GF(2^n)}: \eqn{a \oplus a = 0} (every element is its own additive inverse)
#'
#' @param a Integer scalar representing the first field element. It must lie in
#'   \eqn{\{0,1,\ldots,q-1\}}, where \code{q = gf$q}.
#' @param b Integer scalar representing the second field element. It must lie in
#'   \eqn{\{0,1,\ldots,q-1\}}, where \code{q = gf$q}.
#' @param gf A \code{"galois_field"} object created by \code{\link{gf}}. It
#'   contains the precomputed \eqn{q \times q} addition table used to evaluate
#'   \eqn{a \oplus b}.
#'
#' @return An integer in \eqn{\{0, 1, \ldots, q-1\}} representing \eqn{a \oplus b}.
#'
#' @examples
#' # GF(8) - Binary extension field (XOR addition)
#' gf8 <- gf(8)
#' gf_add(5, 3, gf8)  # 101 XOR 011 = 110 = 6
#' gf_add(7, 7, gf8)  # 111 XOR 111 = 000 = 0 (a + a = 0)
#' gf_add(0, 5, gf8)  # 0 is identity: 0 + 5 = 5
#'
#' # GF(7) - Prime field (modular addition)
#' gf7 <- gf(7)
#' gf_add(3, 5, gf7)  # (3 + 5) mod 7 = 1
#' gf_add(6, 6, gf7)  # (6 + 6) mod 7 = 5
#'
#' # GF(9) - Extension field GF(3^2)
#' gf9 <- gf(9)
#' gf_add(5, 7, gf9)  # Addition in GF(3^2)
#'
#' # Commutative property: a + b = b + a
#' gf_add(3, 5, gf8) == gf_add(5, 3, gf8)  # TRUE
#'
#' @references
#' Lidl, R., & Niederreiter, H. (1997). \emph{Finite Fields} (2nd ed.).
#' Cambridge University Press.
#'
#' Roman, S. (2006). \emph{Field Theory} (2nd ed.). Springer.
#'
#' @seealso
#' \code{\link{gf}}, \code{\link{gf_mult}}, \code{\link{gf_pow}}
#'
#' @family galois_field
#' @export
gf_add <- function(a, b, gf) {
  gf$add[a + 1L, b + 1L]
}


#' @title Field Multiplication in GF(q)
#'
#' @description
#' Computes the product \eqn{a \otimes b} in a Galois Field using precomputed
#' lookup tables for O(1) performance.
#'
#' @details
#' # Multiplication Rules by Field Type
#'
#' ## Prime Fields GF(p)
#' \deqn{a \otimes b = (a \times b) \mod p}
#'
#' ## Extension Fields GF(p^n)
#' Polynomial multiplication modulo the irreducible polynomial.
#'
#' # Example in GF(8)
#'
#' Compute \eqn{5 \otimes 3} in GF(8):
#'
#' \eqn{5 = x^2 + 1} and \eqn{3 = x + 1}
#'
#' \deqn{(x^2 + 1)(x + 1) = x^3 + x^2 + x + 1}
#'
#' Reduce mod \eqn{x^3 + x + 1}:
#' \deqn{x^3 + x^2 + x + 1 \equiv x^2 \pmod{x^3 + x + 1}}
#' (since \eqn{x^3 \equiv x + 1})
#'
#' So \eqn{5 \otimes 3 = 4}.
#'
#' # Key Property
#'
#' Every non-zero element has a multiplicative inverse:
#' \deqn{\forall a \neq 0, \exists a^{-1}: a \otimes a^{-1} = 1}
#'
#' @param a Integer scalar representing the first field element. It must lie in
#'   \eqn{\{0,1,\ldots,q-1\}}, where \code{q = gf$q}.
#' @param b Integer scalar representing the second field element. It must lie in
#'   \eqn{\{0,1,\ldots,q-1\}}, where \code{q = gf$q}.
#' @param gf A \code{"galois_field"} object created by \code{\link{gf}}. It
#'   contains the precomputed \eqn{q \times q} multiplication table used to
#'   evaluate \eqn{a \otimes b}.
#'
#' @return An integer in \eqn{\{0, 1, \ldots, q-1\}} representing \eqn{a \otimes b}.
#'
#' @examples
#' # GF(8) multiplication
#' gf8 <- gf(8)
#' gf_mult(5, 3, gf8)  # (x^2+1)(x+1) = x^2 = 4 (mod x^3+x+1)
#' gf_mult(2, 4, gf8)  # x * x^2 = x^3 = x+1 = 3 (mod x^3+x+1)
#' gf_mult(0, 5, gf8)  # 0 x anything = 0
#' gf_mult(1, 5, gf8)  # 1 is identity: 1 x 5 = 5
#'
#' # GF(7) multiplication (modular)
#' gf7 <- gf(7)
#' gf_mult(3, 5, gf7)  # (3 x 5) mod 7 = 15 mod 7 = 1
#' gf_mult(2, 4, gf7)  # (2 x 4) mod 7 = 8 mod 7 = 1
#'
#' # Note: 3 and 5 are inverses in GF(7)!
#'
#' # Commutative property: a x b = b x a
#' gf_mult(3, 5, gf8) == gf_mult(5, 3, gf8)  # TRUE
#'
#' @references
#' Lidl, R., & Niederreiter, H. (1997). \emph{Finite Fields} (2nd ed.).
#' Cambridge University Press.
#'
#' Roman, S. (2006). \emph{Field Theory} (2nd ed.). Springer.
#'
#' @seealso
#' \code{\link{gf}}, \code{\link{gf_add}}, \code{\link{gf_pow}}
#'
#' @family galois_field
#' @export
gf_mult <- function(a, b, gf) {
  gf$mult[a + 1L, b + 1L]
}


#' @title Field Exponentiation in GF(q)
#'
#' @description
#' Computes \eqn{a^n} in a Galois Field using the efficient square-and-multiply
#' algorithm.
#'
#' @details
#' # Algorithm
#'
#' Uses **square-and-multiply** (binary exponentiation):
#'
#' 1. Write \eqn{n} in binary: \eqn{n = \sum b_i 2^i}
#' 2. Compute \eqn{a^{2^i}} by repeated squaring
#' 3. Multiply together the powers where \eqn{b_i = 1}
#'
#' Time complexity: \eqn{O(\log n)} multiplications.
#'
#' # Fermat's Little Theorem
#'
#' For any \eqn{a \neq 0} in \eqn{GF(q)}:
#' \deqn{a^{q-1} = 1}
#'
#' This is fundamental to finite field theory and implies:
#' \deqn{a^{-1} = a^{q-2}}
#'
#' # Example in GF(8)
#'
#' Since \eqn{|GF(8)^*| = 7}:
#' \deqn{a^7 = 1 \text{ for all } a \neq 0}
#'
#' @param a Integer scalar representing the base field element. It must lie in
#'   \eqn{\{0,1,\ldots,q-1\}}, where \code{q = gf$q}.
#' @param n Non-negative integer scalar giving the exponent. The implementation
#'   applies repeated squaring to compute \eqn{a^n}.
#' @param gf A \code{"galois_field"} object created by \code{\link{gf}}. Its
#'   multiplication table determines the field product used in the power
#'   recursion.
#'
#' @return An integer in \eqn{\{0, 1, \ldots, q-1\}} representing \eqn{a^n}.
#'
#' @examples
#' # GF(8) exponentiation
#' gf8 <- gf(8)
#' gf_pow(2, 3, gf8)  # x^3 = x + 1 = 3 (mod x^3+x+1)
#' gf_pow(2, 7, gf8)  # 2^7 = 1 (Fermat's Little Theorem)
#'
#' # Verify Fermat's theorem for all non-zero elements
#' all(sapply(1:7, function(a) gf_pow(a, 7, gf8)) == 1)  # TRUE
#'
#' # GF(7) exponentiation
#' gf7 <- gf(7)
#' gf_pow(3, 6, gf7)  # 3^6 = 729 mod 7 = 1 (Fermat)
#' gf_pow(3, 3, gf7)  # 3^3 = 27 mod 7 = 6
#'
#' # Special cases
#' gf_pow(5, 0, gf8)  # a^0 = 1
#' gf_pow(5, 1, gf8)  # a^1 = a = 5
#' gf_pow(0, 5, gf8)  # 0^n = 0
#'
#' # Compute multiplicative inverse using Fermat: a^(-1) = a^(q-2)
#' a <- 5
#' a_inv <- gf_pow(a, 6, gf8)  # 5^(-1) in GF(8)
#' gf_mult(a, a_inv, gf8)      # Should be 1
#'
#' @references
#' Lidl, R., & Niederreiter, H. (1997). \emph{Finite Fields} (2nd ed.).
#' Cambridge University Press.
#'
#' Menezes, A. J., van Oorschot, P. C., & Vanstone, S. A. (1996). \emph{Handbook
#' of Applied Cryptography}. CRC Press.
#'
#' @seealso
#' \code{\link{gf}}, \code{\link{gf_add}}, \code{\link{gf_mult}}
#'
#' @family galois_field
#' @export
gf_pow <- function(a, n, gf) {
  if (n == 0L) return(1L)
  if (n == 1L) return(a)
  if (a == 0L) return(0L)

  result <- 1L
  base <- a
  exp <- as.integer(n)

  while (exp > 0L) {
    if (exp %% 2L == 1L) {
      result <- gf_mult(result, base, gf)
    }
    base <- gf_mult(base, base, gf)
    exp <- exp %/% 2L
  }
  result
}


# =============================================================================
# Internal helper functions
# =============================================================================

# Internal: GF(2^n) multiplication
gf2n_mult_internal <- function(a, b, n, irr) {
  if (a == 0L || b == 0L) return(0L)
  result <- 0L
  q <- bitwShiftL(1L, n)
  for (i in 0:(n - 1L)) {
    if (bitwAnd(b, bitwShiftL(1L, i)) != 0L) {
      result <- bitwXor(result, a)
    }
    a <- bitwShiftL(a, 1L)
    if (bitwAnd(a, q) != 0L) {
      a <- bitwXor(a, irr)
    }
  }
  as.integer(bitwAnd(result, q - 1L))
}

# Internal: GF(p^n) addition for p > 2
gfpn_add_internal <- function(a, b, p, n) {
  result <- 0L
  for (i in 0:(n - 1L)) {
    ai <- (a %/% (p^i)) %% p
    bi <- (b %/% (p^i)) %% p
    ri <- (ai + bi) %% p
    result <- result + ri * (p^i)
  }
  as.integer(result)
}

# Internal: GF(p^n) multiplication for p > 2
gfpn_mult_internal <- function(a, b, p, n, irr) {
  if (a == 0L || b == 0L) return(0L)

  a_coeffs <- integer(n)
  b_coeffs <- integer(n)
  for (i in 0:(n - 1L)) {
    a_coeffs[i + 1L] <- (a %/% (p^i)) %% p
    b_coeffs[i + 1L] <- (b %/% (p^i)) %% p
  }

  prod_coeffs <- integer(2 * n - 1)
  for (i in 1:n) {
    for (j in 1:n) {
      prod_coeffs[i + j - 1] <- (prod_coeffs[i + j - 1] + a_coeffs[i] * b_coeffs[j]) %% p
    }
  }

  for (i in (2 * n - 1):(n + 1)) {
    if (prod_coeffs[i] != 0L) {
      coef <- prod_coeffs[i]
      for (j in 1:length(irr)) {
        idx <- i - n + j - 1
        if (idx >= 1 && idx <= length(prod_coeffs)) {
          prod_coeffs[idx] <- (prod_coeffs[idx] - coef * irr[j]) %% p
          if (prod_coeffs[idx] < 0) prod_coeffs[idx] <- prod_coeffs[idx] + p
        }
      }
    }
  }

  result <- 0L
  for (i in 1:n) {
    result <- result + prod_coeffs[i] * (p^(i - 1))
  }
  as.integer(result)
}

# Internal: Integer factorization for error messages
factorize <- function(n) {
  factors <- integer(0)
  d <- 2L
  while (n > 1) {
    while (n %% d == 0) {
      factors <- c(factors, d)
      n <- n %/% d
    }
    d <- d + 1L
    if (d * d > n && n > 1) {
      factors <- c(factors, n)
      break
    }
  }
  factors
}


#' @title Print a Galois Field Object
#'
#' @description
#' Prints a compact summary of a \code{"galois_field"} object created by
#' \code{\link{gf}}.
#'
#' @details
#' Let \eqn{x = \mathrm{GF}(q)} be the field object returned by
#' \code{\link{gf}}, where
#' \deqn{
#'   q = p^n
#' }
#' with prime characteristic \eqn{p} and extension degree \eqn{n \ge 1}. The
#' printed summary reports:
#'
#' - the field description,
#' - the order \eqn{q},
#' - the characteristic \eqn{p}, and
#' - the extension degree \eqn{n}.
#'
#' This is a display method only; it does not modify the object or its lookup
#' tables.
#'
#' @param x An object of class \code{"galois_field"} returned by
#'   \code{\link{gf}}.
#' @param ... Unused arguments passed through the generic \code{print()}
#'   interface.
#'
#' @return The input object \code{x}, invisibly.
#'
#' @examples
#' gf8 <- gf(8)
#' print(gf8)
#'
#' gf7 <- gf(7)
#' print(gf7)
#'
#' @references
#' Lidl, R., & Niederreiter, H. (1997). \emph{Finite Fields} (2nd ed.).
#' Cambridge University Press.
#'
#' @seealso
#' \code{\link{gf}}, \code{\link{gf_add}}, \code{\link{gf_mult}},
#' \code{\link{gf_pow}}
#'
#' @family galois_field
#' @export
print.galois_field <- function(x, ...) {
  cat("Galois Field", x$description, "\n")
  cat("  Order q =", x$q, "elements\n")
  cat("  Characteristic p =", x$p, "\n")
  cat("  Extension degree n =", x$n, "\n")
  invisible(x)
}

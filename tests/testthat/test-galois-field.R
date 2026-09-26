test_that("gf(8) has exact lookup-table values used by PP examples", {
  gf8 <- gf(8)

  expect_equal(gf8$q, 8L)
  expect_equal(gf8$p, 2L)
  expect_equal(gf8$n, 3L)
  expect_equal(gf8$description, "GF(2^3) with irreducible x^3 + x + 1")

  # Exact rows from the current lookup tables.
  expect_equal(gf8$add[3, ], c(2L, 3L, 0L, 1L, 6L, 7L, 4L, 5L))
  expect_equal(gf8$mult[3, ], c(0L, 2L, 4L, 6L, 3L, 1L, 7L, 5L))

  # Known worked examples in the docs.
  expect_equal(gf_add(5L, 3L, gf8), 6L)
  expect_equal(gf_mult(2L, 4L, gf8), 3L)
  expect_equal(gf_mult(5L, 3L, gf8), 4L)
  expect_equal(gf_pow(2L, 7L, gf8), 1L)
})

test_that("prime and prime-power helpers cover exact outputs and edge cases", {
  expect_true(is_prime(2L))
  expect_true(is_prime(97L))
  expect_false(is_prime(1L))
  expect_false(is_prime(100L))

  expect_equal(
    is_prime_power(8L),
    list(is_prime_power = TRUE, p = 2L, n = 3L)
  )
  expect_equal(
    is_prime_power(49L),
    list(is_prime_power = TRUE, p = 7L, n = 2L)
  )
  expect_equal(
    is_prime_power(6L),
    list(is_prime_power = FALSE, p = NA_integer_, n = NA_integer_)
  )
  expect_equal(
    is_prime_power(1L),
    list(is_prime_power = FALSE, p = NA_integer_, n = NA_integer_)
  )
})

test_that("prime and extension fields satisfy exact algebraic identities", {
  gf7 <- gf(7)
  gf9 <- gf(9)

  expect_equal(gf_add(3L, 5L, gf7), 1L)
  expect_equal(gf_mult(3L, 5L, gf7), 1L)
  expect_equal(gf_pow(3L, 6L, gf7), 1L)

  expect_equal(gf_add(5L, 7L, gf9), gf9$add[6, 8])
  expect_equal(gf_mult(4L, 8L, gf9), gf9$mult[5, 9])
  expect_equal(gf_pow(5L, 0L, gf9), 1L)
  expect_equal(gf_pow(0L, 5L, gf9), 0L)
})

test_that("gf() and print.galois_field() handle invalid and display paths", {
  expect_error(gf(6L), "prime power")
  expect_equal(factorize(60L), c(2L, 2L, 3L, 5L))

  gf8 <- gf(8)
  expect_output(print(gf8), "Galois Field GF\\(2\\^3\\) with irreducible x\\^3 \\+ x \\+ 1")
  expect_invisible(print(gf8))
})

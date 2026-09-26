# Comprehensive test suite for ppgt package

library(testthat)

# =============================================================================
# GALOIS FIELD TESTS
# =============================================================================

test_that("is_prime correctly identifies primes", {
  expect_true(is_prime(2))
  expect_true(is_prime(3))
  expect_true(is_prime(7))
  expect_true(is_prime(97))
  expect_true(is_prime(104729))  # 10000th prime
  
  expect_false(is_prime(0))
  expect_false(is_prime(1))
  expect_false(is_prime(4))
  expect_false(is_prime(9))
  expect_false(is_prime(100))
})

test_that("is_prime_power identifies prime powers correctly", {
  # Valid prime powers
  expect_true(is_prime_power(2)$is_prime_power)
  expect_true(is_prime_power(4)$is_prime_power)
  expect_true(is_prime_power(8)$is_prime_power)
  expect_true(is_prime_power(9)$is_prime_power)
  expect_true(is_prime_power(27)$is_prime_power)
  expect_true(is_prime_power(49)$is_prime_power)
  expect_true(is_prime_power(64)$is_prime_power)
  
  # Check decomposition
  expect_equal(is_prime_power(8)$p, 2)
  expect_equal(is_prime_power(8)$n, 3)
  expect_equal(is_prime_power(9)$p, 3)
  expect_equal(is_prime_power(9)$n, 2)
  expect_equal(is_prime_power(49)$p, 7)
  expect_equal(is_prime_power(49)$n, 2)
  
  # Invalid (not prime powers)
  expect_false(is_prime_power(6)$is_prime_power)
  expect_false(is_prime_power(10)$is_prime_power)
  expect_false(is_prime_power(12)$is_prime_power)
  expect_false(is_prime_power(15)$is_prime_power)
  expect_false(is_prime_power(18)$is_prime_power)
})

test_that("GF(8) arithmetic is correct", {
  gf8 <- gf(8)
  
  expect_equal(gf8$q, 8)
  expect_equal(gf8$p, 2)
  expect_equal(gf8$n, 3)
  
  # Addition (XOR)
  expect_equal(gf_add(0, 0, gf8), 0)
  expect_equal(gf_add(5, 0, gf8), 5)
  expect_equal(gf_add(5, 5, gf8), 0)  # a + a = 0 in GF(2^n)
  expect_equal(gf_add(5, 3, gf8), 6)  # 101 XOR 011 = 110
  expect_equal(gf_add(7, 7, gf8), 0)
  
  # Multiplication
  expect_equal(gf_mult(0, 5, gf8), 0)
  expect_equal(gf_mult(1, 5, gf8), 5)
  expect_equal(gf_mult(2, 2, gf8), 4)  # x * x = x²
  expect_equal(gf_mult(2, 4, gf8), 3)  # x * x² = x³ = x + 1
  expect_equal(gf_mult(5, 3, gf8), 4)
  
  # Commutative
  expect_equal(gf_add(3, 5, gf8), gf_add(5, 3, gf8))
  expect_equal(gf_mult(3, 5, gf8), gf_mult(5, 3, gf8))
  
  # Fermat's Little Theorem: a^7 = 1 for all a ≠ 0
  for (a in 1:7) {
    expect_equal(gf_pow(a, 7, gf8), 1)
  }
})

test_that("GF(9) arithmetic is correct", {
  gf9 <- gf(9)
  
  expect_equal(gf9$q, 9)
  expect_equal(gf9$p, 3)
  expect_equal(gf9$n, 2)
  
  # Basic operations
  expect_equal(gf_mult(0, 5, gf9), 0)
  expect_equal(gf_mult(1, 5, gf9), 5)
  
  # Fermat: a^8 = 1 for all a ≠ 0 in GF(9)
  for (a in 1:8) {
    expect_equal(gf_pow(a, 8, gf9), 1)
  }
})

test_that("Prime field GF(7) works", {
  gf7 <- gf(7)
  
  expect_equal(gf7$q, 7)
  expect_equal(gf7$p, 7)
  expect_equal(gf7$n, 1)
  
  # Modular arithmetic
  expect_equal(gf_add(3, 5, gf7), 1)   # (3+5) mod 7
  expect_equal(gf_mult(3, 5, gf7), 1)  # (15) mod 7
  expect_equal(gf_pow(3, 6, gf7), 1)   # Fermat
})

test_that("Invalid field order throws error", {
  expect_error(gf(6), "prime power")
  expect_error(gf(10), "prime power")
  expect_error(gf(12), "prime power")
})

# =============================================================================
# PP MATRIX TESTS
# =============================================================================

test_that("PP matrix has correct dimensions", {
  # q=4, d=3, nl=5: 20 pools × 64 samples
  M1 <- pp_matrix(q = 4, d = 3, nl = 5)
  expect_equal(nrow(M1), 20)
  expect_equal(ncol(M1), 64)
  
  # P-BEST: 48 pools × 384 samples
  M2 <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  expect_equal(nrow(M2), 48)
  expect_equal(ncol(M2), 384)
  
  # GF(9): 63 pools × 729 samples
  M3 <- pp_matrix(q = 9, d = 3, nl = 7)
  expect_equal(nrow(M3), 63)
  expect_equal(ncol(M3), 729)
})

test_that("PP matrix has uniform structure", {
  M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  
  # Uniform pool sizes
  pool_sizes <- as.vector(Matrix::rowSums(M))
  expect_equal(length(unique(pool_sizes)), 1)
  
  # Uniform sample coverage
  sample_cov <- as.vector(Matrix::colSums(M))
  expect_equal(length(unique(sample_cov)), 1)
  expect_equal(unique(sample_cov), 6)  # nl = 6
})

test_that("PP matrix stores correct parameters", {
  M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  params <- attr(M, "pp")
  
  expect_equal(params$q, 8)
  expect_equal(params$d, 3)
  expect_equal(params$nl, 6)
  expect_equal(params$N, 384)
  expect_equal(params$M, 48)
  expect_equal(params$pool_size, 48L)
  expect_equal(params$k_max, 2)
})

test_that("pp_design returns valid configurations", {
  configs <- pp_design(N_min = 300, N_max = 500, k_min = 2)
  
  expect_true(nrow(configs) > 0)
  expect_true(all(configs$k_max >= 2))
  expect_true(all(configs$compression > 1))
  
  # Check that returned configs can generate valid matrices
  row1 <- configs[1, ]
  M <- pp_matrix(q = row1$q, d = row1$d, nl = row1$nl)
  expect_equal(nrow(M), row1$M)
})

test_that("Different q values produce valid matrices", {
  for (q in c(3, 4, 5, 7, 8, 9)) {
    M <- pp_matrix(q = q, d = 3, nl = 3)
    expect_equal(nrow(M), 3 * q)
    expect_equal(ncol(M), q^3)
    
    # Check uniformity
    pool_sizes <- as.vector(Matrix::rowSums(M))
    expect_equal(length(unique(pool_sizes)), 1)
  }
})

# =============================================================================
# DECODING TESTS
# =============================================================================

test_that("COMP identifies positives correctly (Julia example)", {
  # Julia example_1.jl: q=4, d=3, positives at 2 and 13
  M <- pp_matrix(q = 4, d = 3, nl = 5)
  x <- rep(0L, 64)
  x[c(2, 13)] <- 1L
  y <- as.integer(as.vector(M %*% x) > 0)
  
  result <- pp_comp(M, y)
  expect_equal(sort(result$candidates), c(2, 13))
})

test_that("COMP has no false negatives", {
  M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  
  set.seed(42)
  for (k in 1:3) {
    true_pos <- sample(384, k)
    x <- rep(0L, 384)
    x[true_pos] <- 1L
    y <- as.integer(as.vector(M %*% x) > 0)
    
    result <- pp_comp(M, y)
    
    # All true positives must be in candidates
    expect_true(all(true_pos %in% result$candidates))
  }
})

test_that("pp_decode identifies positives exactly when k <= k_max", {
  M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  
  # k = 2 (within k_max = 2)
  true_pos <- c(72, 142)
  x <- rep(0L, 384)
  x[true_pos] <- 1L
  y <- as.integer(as.vector(M %*% x) > 0)
  
  result <- pp_decode(M, y, verbose = FALSE)
  expect_equal(sort(result$positives), sort(true_pos))
  expect_equal(result$error, 0)
})

test_that("pp_validate computes correct metrics", {
  # Perfect match
  m1 <- pp_validate(c(1, 2, 3), c(1, 2, 3))
  expect_equal(m1$TP, 3)
  expect_equal(m1$FP, 0)
  expect_equal(m1$FN, 0)
  expect_equal(m1$sensitivity, 1)
  expect_equal(m1$precision, 1)
  expect_true(m1$exact_match)
  
  # False positives
  m2 <- pp_validate(c(1, 2, 3, 4), c(1, 2, 3))
  expect_equal(m2$TP, 3)
  expect_equal(m2$FP, 1)
  expect_equal(m2$FN, 0)
  expect_false(m2$exact_match)
  
  # False negatives
  m3 <- pp_validate(c(1, 2), c(1, 2, 3))
  expect_equal(m3$TP, 2)
  expect_equal(m3$FP, 0)
  expect_equal(m3$FN, 1)
})

test_that("pp_gpsr produces valid output", {
  A <- matrix(c(1, 1, 0, 0, 1, 1, 1, 0, 1), nrow = 3, byrow = TRUE)
  y <- c(1, 0.8, 0.5)
  
  result <- pp_gpsr(A, y)
  
  expect_equal(length(result$x), 3)
  expect_true(all(result$x >= 0))
  expect_true(is.logical(result$converged))
})

# =============================================================================
# INTEGRATION TESTS
# =============================================================================

test_that("Complete workflow works", {
  # 1. Design
  configs <- pp_design(N_min = 50, N_max = 100)
  expect_true(nrow(configs) > 0)
  
  # 2. Generate matrix
  M <- pp_matrix(q = 4, d = 3, nl = 5)
  
  # 3. Verify
  info <- pp_verify(M)
  expect_true(info$uniform_pools)
  expect_true(info$uniform_coverage)
  
  # 4. Simulate and decode
  x <- rep(0L, 64)
  x[c(10, 20)] <- 1L
  y <- as.integer(as.vector(M %*% x) > 0)
  
  result <- pp_decode(M, y, verbose = FALSE)
  
  # 5. Validate
  metrics <- pp_validate(result$positives, c(10, 20))
  expect_true(metrics$sensitivity == 1)  # No false negatives
})

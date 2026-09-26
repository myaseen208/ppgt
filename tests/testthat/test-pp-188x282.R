test_that("pp_matrix(q = 47, d = 3, nl = 4, N = 282) matches the vignette claims", {
  # Protects the numeric claims made in the ppgt_pp_188x282 vignette
  # ("Building a 188 x 282 pooling matrix"): 282 specimens in 188 pools,
  # 6 specimens per pool, 4 pools per specimen, and a pairwise overlap
  # bounded by d - 1 = 2.
  M4 <- pp_matrix(q = 47, d = 3, nl = 4, N = 282)

  expect_equal(dim(M4), c(188L, 282L))
  expect_true(all(Matrix::rowSums(M4) == 6L))
  expect_true(all(Matrix::colSums(M4) == 4L))

  total_incidences <- sum(M4)
  expect_equal(total_incidences, 188L * 6L)
  expect_equal(total_incidences, 282L * 4L)

  overlap <- Matrix::crossprod(M4)
  Matrix::diag(overlap) <- 0L
  expect_true(max(overlap) <= 2L)  # d - 1 bound

  expect_equal(attr(M4, "pp")$k_max, 1L)
})

test_that("pp_decode recovers a single simulated positive in the 188 x 282 design", {
  M4 <- pp_matrix(q = 47, d = 3, nl = 4, N = 282)

  Y_tilde <- integer(282L)
  Y_tilde[137L] <- 1L
  z_tilde <- as.integer(as.vector(M4 %*% Y_tilde) > 0L)

  decoded <- pp_decode(M4, z_tilde, verbose = FALSE)
  expect_identical(decoded$positives, 137L)
  expect_equal(decoded$method, "COMP (exact)")
})

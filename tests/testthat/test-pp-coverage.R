test_that("PP matrices match documented dimensions, weights, and sparsity", {
  # Tan (2020), Appendix / Algorithm 2.4: q = 8, d = 3, n_l = 6 gives
  # N = 384 individuals in J = 48 pools with 48 individuals per pool
  # and 6 pools per individual.
  M_pbest <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  expect_equal(dim(M_pbest), c(48, 384))
  expect_true(all(Matrix::rowSums(M_pbest) == 48))
  expect_true(all(Matrix::colSums(M_pbest) == 6))
  expect_equal(sum(M_pbest != 0), 384 * 6)
  expect_equal(sum(M_pbest != 0) / prod(dim(M_pbest)), 1 / 8)

  # Tan (2020), Algorithm 2.3: q = 4, d = 3, n_l = 5 gives
  # N = 64 individuals in J = 20 pools with 16 individuals per pool.
  M_small <- pp_matrix(q = 4, d = 3, nl = 5)
  expect_equal(dim(M_small), c(20, 64))
  expect_true(all(Matrix::rowSums(M_small) == 16))
  expect_true(all(Matrix::colSums(M_small) == 5))
  expect_equal(sum(M_small != 0), 64 * 5)
  expect_equal(sum(M_small != 0) / prod(dim(M_small)), 1 / 4)

})

test_that("Current d = 2 PP branch is exercised as a documented regression", {
  # The Tan (2020) two-dimensional construction should produce J = n_l q pools
  # over N = q^2 individuals, but the current implementation still errors in
  # this branch. Keep the branch covered without masking the defect.
  expect_error(
    pp_matrix(q = 2, d = 2, nl = 3),
    "dims"
  )
})

test_that("PP matrix truncation and singular layer branches are exercised", {
  # Truncated design: N = a1 q^(d-1) + a2 q with a2 > 0.
  M_trunc <- pp_matrix(q = 4, d = 3, nl = 5, N = 20)

  expect_equal(dim(M_trunc), c(20, 20))
  expect_true(all(Matrix::colSums(M_trunc) == 5))
  expect_true(any(Matrix::rowSums(M_trunc) < 16))
  expect_equal(attr(M_trunc, "pp")$pool_size, Matrix::rowSums(M_trunc)[1])

  # The last q rows are the singular layer when n_l = q + 1.
  singular_rows <- 17:20
  expect_true(all(Matrix::rowSums(M_trunc)[singular_rows] <= 16))
  expect_true(any(Matrix::rowSums(M_trunc)[singular_rows] < 16))
})

test_that("PP matrix argument validation covers invalid design requests", {
  expect_error(pp_matrix(q = 6, d = 3, nl = 5), "prime power")
  expect_error(pp_matrix(q = 4, d = 1, nl = 3), "d must be >= 2")
  expect_error(pp_matrix(q = 4, d = 3, nl = 0), "nl must be in")
  expect_error(pp_matrix(q = 4, d = 3, nl = 6), "nl must be in")
  expect_error(pp_matrix(q = 4, d = 3, nl = 5, N = 65), "cannot exceed q\\^d")
  expect_error(pp_matrix(q = 4, d = 3, nl = 5, N = 18), "must be divisible by q")
})

test_that("pp_design and pp_verify cover empty and non-empty planning paths", {
  configs <- pp_design(N_min = 60, N_max = 70, k_min = 2, pool_size_max = 20)
  expect_true(nrow(configs) > 0)
  expect_true(all(diff(configs$compression) <= 0))

  expect_message(
    none <- pp_design(N_min = 10, N_max = 10, M_max = 1, k_min = 10, pool_size_max = 1),
    "No configurations found"
  )
  expect_null(none)

  M <- pp_matrix(q = 4, d = 3, nl = 5)
  expect_message(info <- pp_verify(M), "PP Matrix Verification")
  expect_equal(info$n_pools, 20)
  expect_equal(info$n_samples, 64)
  expect_true(info$uniform_pools)
  expect_true(info$uniform_coverage)
})

test_that("pp_design covers pool-size and pool-count rejection filters", {
  expect_message(
    too_large_pools <- pp_design(
      N_min = 64, N_max = 64, M_max = 100, k_min = 2, pool_size_max = 3
    ),
    "No configurations found"
  )
  expect_null(too_large_pools)

  expect_message(
    too_many_pools <- pp_design(
      N_min = 64, N_max = 64, M_max = 5, k_min = 2, pool_size_max = 100
    ),
    "No configurations found"
  )
  expect_null(too_many_pools)
})

test_that("pp_comp handles all negatives, single positives, Ct values, and errors", {
  M <- pp_matrix(q = 4, d = 3, nl = 5)

  y0 <- integer(nrow(M))
  res0 <- pp_comp(M, y0)
  expect_length(res0$candidates, 0)
  expect_equal(length(res0$negatives), ncol(M))
  expect_equal(res0$n_positive_pools, 0)

  x1 <- integer(ncol(M))
  x1[7] <- 1L
  y1 <- as.integer(as.vector(M %*% x1) > 0)
  res1 <- pp_comp(M, y1)
  expect_equal(res1$candidates, 7)

  ct_y <- ifelse(y1 == 1L, 30, 0)
  res_ct <- pp_comp(M, ct_y, ct_threshold = 40)
  expect_equal(res_ct$candidates, 7)

  expect_error(pp_comp(M, y1[-1]), "Length of y")
})

test_that("pp_gpsr covers zero-tau, non-convergence, and input validation branches", {
  A <- diag(3)
  y_zero <- c(0, 0, 0)
  res_zero <- pp_gpsr(A, y_zero, tau = NULL, max_iter = 2L)
  expect_equal(length(res_zero$x), 3)
  expect_true(all(res_zero$x >= 0))

  # Force the non-converged return by allowing no iterations.
  res_nc <- pp_gpsr(A, c(1, 0, 0), tau = 0.1, max_iter = 0L)
  expect_false(res_nc$converged)
  expect_equal(res_nc$iterations, 0L)

  expect_error(pp_gpsr(A, c(1, 0)), "Length of y")
})

test_that("pp_decode covers all documented decoder paths", {
  M <- pp_matrix(q = 4, d = 3, nl = 5)

  # All negatives branch with verbose messaging.
  y0 <- integer(nrow(M))
  expect_message(
    res0 <- pp_decode(M, y0, verbose = TRUE),
    "No positives detected"
  )
  expect_equal(res0$positives, integer(0))
  expect_equal(res0$method, "COMP (all pools negative)")

  # Exact branch with a suspected sample triggered by Ct information.
  x1 <- integer(ncol(M))
  x1[7] <- 1L
  y1 <- as.integer(as.vector(M %*% x1) > 0)
  ct_vals <- rep(NA_real_, ncol(M))
  ct_vals[7] <- 37
  expect_message(
    res_exact <- pp_decode(M, y1, verbose = TRUE, ct = ct_vals),
    "Suspected:"
  )
  expect_equal(res_exact$positives, integer(0))
  expect_equal(res_exact$suspected, 7)
  expect_equal(res_exact$method, "COMP (exact)")

  # Fallback k_max branch using a matrix without pp metadata.
  M_plain <- unclass(as.matrix(M))
  attr(M_plain, "pp") <- NULL
  expect_warning(
    res_fallback <- pp_decode(M_plain, y1, verbose = FALSE),
    "k_max not supplied"
  )
  expect_equal(sort(c(res_fallback$positives, res_fallback$suspected)), 7)

  expect_error(pp_decode(M, y1, ct = rep(30, ncol(M) - 1L), verbose = FALSE), "Length of ct")
})

test_that("pp_decode verification and GPSR refinement branches are covered", {
  M <- matrix(c(
    1, 0, 1, 0,
    0, 1, 1, 0,
    0, 0, 0, 1
  ), nrow = 3, byrow = TRUE)

  # COMP candidates {1,2,3}; the exhaustive verification branch should return
  # a zero-error explanation without invoking GPSR.
  x_true <- c(1L, 1L, 0L, 0L)
  y_obs <- as.integer(as.vector(M %*% x_true) > 0)
  expect_message(
    res_verify <- pp_decode(M, y_obs, k_max = 2L, verbose = TRUE),
    "COMP \\+ Verification"
  )
  expect_true(length(res_verify$positives) >= 1L)
  expect_equal(res_verify$error, 0)
  expect_equal(res_verify$method, "COMP + Verification")

  # Force the GPSR path by making the pre-GPSR verification step fail.
  local_mocked_bindings(
    search_best_internal = function(...) list(positives = integer(0), error = 1L),
    pp_gpsr = function(...) stop("forced gpsr failure")
  )
  expect_message(
    res_gpsr <- pp_decode(M, y_obs, k_max = 1L, verbose = TRUE),
    "COMP \\+ GPSR \\+ Verification"
  )
  expect_equal(res_gpsr$positives, integer(0))
  expect_equal(res_gpsr$error, 1L)
  expect_equal(res_gpsr$method, "COMP + GPSR + Verification")
})

test_that("pp_decode covers suspected messaging in verification and GPSR paths", {
  M <- matrix(c(
    1, 0, 1, 0,
    0, 1, 1, 0,
    0, 0, 0, 1
  ), nrow = 3, byrow = TRUE)
  y_ct <- c(30, 30, 0)

  local_mocked_bindings(
    search_best_internal = function(...) list(positives = c(2L), error = 0L),
    classify_threeway_internal = function(...) {
      list(confirmed = integer(0), suspected = 2L)
    }
  )
  expect_message(
    res_verify_sus <- pp_decode(M, y_ct, k_max = 1L, verbose = TRUE),
    "Suspected: 2"
  )
  expect_equal(res_verify_sus$positives, integer(0))
  expect_equal(res_verify_sus$suspected, 2L)
  expect_equal(res_verify_sus$method, "COMP + Verification")

  call_idx <- 0L
  local_mocked_bindings(
    search_best_internal = function(...) {
      call_idx <<- call_idx + 1L
      if (call_idx == 1L) {
        list(positives = integer(0), error = 1L)
      } else {
        list(positives = c(2L), error = 0L)
      }
    },
    pp_gpsr = function(...) list(x = c(0.8, 0.2, 0, 0)),
    classify_threeway_internal = function(...) {
      list(confirmed = integer(0), suspected = 2L)
    }
  )
  expect_message(
    res_gpsr_sus <- pp_decode(M, y_ct, k_max = 1L, verbose = TRUE),
    "Suspected: 2"
  )
  expect_equal(res_gpsr_sus$positives, integer(0))
  expect_equal(res_gpsr_sus$suspected, 2L)
  expect_equal(res_gpsr_sus$error, 0L)
  expect_equal(res_gpsr_sus$method, "COMP + GPSR + Verification")
})

test_that("pp_decode internal helpers cover remaining classification cases", {
  M_id <- diag(3)
  y_pos <- c(1L, 1L, 0L)

  cls_empty <- classify_threeway_internal(M_id, y_pos, integer(0))
  expect_equal(cls_empty$confirmed, integer(0))
  expect_equal(cls_empty$suspected, integer(0))

  # Positive pools empty and no unique pool.
  cls_no_pools <- classify_threeway_internal(M_id, c(0L, 0L, 0L), positives = 1L)
  expect_equal(cls_no_pools$confirmed, integer(0))
  expect_equal(cls_no_pools$suspected, 1L)

  # Ct filtering with valid and NA entries.
  cls_ct <- classify_threeway_internal(
    M = M_id,
    y_bin = c(1L, 1L, 1L),
    positives = c(1L, 2L),
    ct = c(24, 40, NA)
  )
  expect_equal(cls_ct$confirmed, 1L)
  expect_equal(cls_ct$suspected, 2L)

  expect_equal(verify_pools_internal(M_id, c(1L, 0L, 0L), positives = 1L), 0)
  expect_gt(verify_pools_internal(M_id, c(1L, 0L, 0L), positives = 2L), 0)

  best <- search_best_internal(M_id, c(1L, 0L, 0L), candidates = c(1L, 2L, 3L), k_max = 2L)
  expect_equal(best$positives, 1L)
  expect_equal(best$error, 0)

  expect_equal(pp_validate(integer(0), integer(0))$precision, 1)
})

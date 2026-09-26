test_that("SeparableMatrix reed-solomon branch is d-separable for d = 2", {
  mat <- SeparableMatrix(N = 16L, d = 2L, method = "reed-solomon", verify = TRUE)
  params <- .sm_rs_params(16L, 2L)

  expect_s4_class(mat, "dgCMatrix")
  expect_equal(dim(mat), c(params$J, 16L))
  expect_equal(attr(mat, "d"), 2L)
  expect_equal(attr(mat, "method_used"), "reed-solomon")
  expect_true(all(Matrix::colSums(mat) == params$L))
  expect_true(isTRUE(as.logical(attr(mat, "separable_verified"))))
  expect_equal(attr(attr(mat, "separable_verified"), "verified"), "exact")
})

test_that("SeparableMatrix reed-solomon branch is d-separable for d = 3", {
  mat <- SeparableMatrix(N = 27L, d = 3L, method = "reed-solomon", verify = TRUE)
  params <- .sm_rs_params(27L, 3L)

  expect_equal(dim(mat), c(params$J, 27L))
  expect_equal(attr(mat, "method_used"), "reed-solomon")
  expect_true(all(Matrix::colSums(mat) == params$L))
  expect_true(isTRUE(as.logical(attr(mat, "separable_verified"))))
  expect_equal(attr(attr(mat, "separable_verified"), "verified"), "exact")
})

test_that("SeparableMatrix reed-solomon ignores user-supplied M", {
  mat <- SeparableMatrix(N = 16L, d = 2L, M = 999L, method = "reed-solomon", verify = TRUE)
  params <- .sm_rs_params(16L, 2L)

  expect_equal(nrow(mat), params$J)
  expect_true(isTRUE(as.logical(attr(mat, "separable_verified"))))
})

test_that("DisjunctMatrix reed-solomon branch is d-disjunct for d = 2", {
  mat <- DisjunctMatrix(N = 16L, d = 2L, method = "reed-solomon", verify = TRUE)
  params <- .sm_rs_params(16L, 2L)

  expect_s4_class(mat, "dgCMatrix")
  expect_equal(dim(mat), c(params$J, 16L))
  expect_equal(attr(mat, "d"), 2L)
  expect_equal(attr(mat, "method_used"), "reed-solomon")
  expect_true(all(Matrix::colSums(mat) == params$L))
  expect_true(isTRUE(as.logical(attr(mat, "disjunct_verified"))))
  expect_equal(attr(attr(mat, "disjunct_verified"), "verified"), "exact")
})

test_that("DisjunctMatrix auto uses reed-solomon when q > d", {
  mat <- DisjunctMatrix(N = 16L, d = 2L, method = "auto", verify = TRUE)

  expect_equal(attr(mat, "method_used"), "reed-solomon")
  expect_true(isTRUE(as.logical(attr(mat, "disjunct_verified"))))
})

test_that("SeparableMatrix validates arguments and random auto path", {
  expect_error(SeparableMatrix(N = 0L, d = 1L), "N must be a positive integer")
  expect_error(SeparableMatrix(N = 5L, d = 0L), "d must be a positive integer")
  expect_error(SeparableMatrix(N = 5L, d = 5L), "d must be less than N")
  expect_error(SeparableMatrix(N = 5L, d = 1L, M = 0L), "M must be a positive integer")

  mat <- SeparableMatrix(N = 20L, d = 4L, method = "auto", seed = 7L, verify = FALSE)
  expect_equal(attr(mat, "method_used"), "random")
  expect_equal(attr(mat, "seed"), 7L)
  expect_null(attr(mat, "separable_verified"))
  expect_equal(dim(mat), c(as.integer(ceiling(2 * 4 * log2(20))), 20L))
})

test_that("SeparableMatrix random branch retries and warns on failed verification", {
  fake_mat <- Matrix::Matrix(0L, nrow = 3, ncol = 5, sparse = TRUE)
  calls <- integer(0)

  local_mocked_bindings(
    .build_random_sep_matrix = function(N, d, M, seed) {
      calls <<- c(calls, seed)
      fake_mat
    },
    .is_d_separable = function(mat, d) FALSE
  )

  expect_warning(
    mat <- SeparableMatrix(N = 5L, d = 2L, M = 3L, method = "random", seed = 11L, verify = TRUE),
    "random construction did not verify as 2-separable after 10 attempts"
  )

  expect_equal(calls, 11:20)
  expect_equal(as.matrix(mat), as.matrix(fake_mat))
  expect_false(isTRUE(as.logical(attr(mat, "separable_verified"))))
  expect_equal(attr(mat, "method_used"), "random")
  expect_equal(attr(mat, "seed"), 20L)
})

test_that("DisjunctMatrix validates arguments and random branches", {
  expect_error(DisjunctMatrix(N = 0L, d = 1L), "N must be a positive integer")
  expect_error(DisjunctMatrix(N = 5L, d = 0L), "d must be a positive integer")
  expect_error(DisjunctMatrix(N = 5L, d = 5L), "d must be less than N")
  expect_error(DisjunctMatrix(N = 5L, d = 1L, M = 0L), "M must be a positive integer")
  
  local_mocked_bindings(
    .sm_rs_params = function(N, d) list(q = 4L),
    .build_random_disjunct_matrix = function(N, d, M, seed) {
      Matrix::Matrix(0L, nrow = M, ncol = N, sparse = TRUE)
    }
  )
  expect_error(
    DisjunctMatrix(N = 64L, d = 4L, method = "reed-solomon", verify = FALSE),
    "reed-solomon construction requires q > d"
  )

  mat_auto <- DisjunctMatrix(N = 64L, d = 4L, method = "auto", seed = 9L, verify = FALSE)
  expect_equal(attr(mat_auto, "method_used"), "random")
  expect_equal(attr(mat_auto, "seed"), 9L)
  expect_null(attr(mat_auto, "disjunct_verified"))
  expect_equal(dim(mat_auto), c(as.integer(ceiling((4 + 1)^2 * log2(64))), 64L))

  mat_random <- DisjunctMatrix(N = 32L, d = 2L, M = 15L, method = "random", seed = 4L, verify = FALSE)
  expect_equal(dim(mat_random), c(15L, 32L))
  expect_equal(attr(mat_random, "method_used"), "random")
  expect_equal(attr(mat_random, "seed"), 4L)
})

test_that("DisjunctMatrix random branch retries and warns on failed verification", {
  fake_mat <- Matrix::Matrix(0L, nrow = 4, ncol = 6, sparse = TRUE)
  calls <- integer(0)

  local_mocked_bindings(
    .build_random_disjunct_matrix = function(N, d, M, seed) {
      calls <<- c(calls, seed)
      fake_mat
    },
    .is_d_disjunct = function(mat, d) FALSE
  )

  expect_warning(
    mat <- DisjunctMatrix(N = 6L, d = 2L, M = 4L, method = "random", seed = 21L, verify = TRUE),
    "random construction did not verify as 2-disjunct after 10 attempts"
  )

  expect_equal(calls, 21:30)
  expect_equal(as.matrix(mat), as.matrix(fake_mat))
  expect_false(isTRUE(as.logical(attr(mat, "disjunct_verified"))))
  expect_equal(attr(mat, "method_used"), "random")
  expect_equal(attr(mat, "seed"), 30L)
})

test_that("SeparableMatrix random branch stops on first verified success", {
  fake_mat <- Matrix::Matrix(1L, nrow = 3, ncol = 5, sparse = TRUE)
  calls <- integer(0)
  verify_calls <- 0L

  local_mocked_bindings(
    .build_random_sep_matrix = function(N, d, M, seed) {
      calls <<- c(calls, seed)
      fake_mat
    },
    .is_d_separable = function(mat, d) {
      verify_calls <<- verify_calls + 1L
      verify_calls >= 2L
    }
  )

  mat <- SeparableMatrix(N = 5L, d = 2L, M = 3L, method = "random", seed = 31L, verify = TRUE)

  expect_equal(calls, 31:32)
  expect_true(isTRUE(as.logical(attr(mat, "separable_verified"))))
  expect_equal(attr(mat, "seed"), 32L)
  expect_equal(attr(mat, "method_used"), "random")
})

test_that("DisjunctMatrix random branch stops on first verified success", {
  fake_mat <- Matrix::Matrix(1L, nrow = 4, ncol = 6, sparse = TRUE)
  calls <- integer(0)
  verify_calls <- 0L

  local_mocked_bindings(
    .build_random_disjunct_matrix = function(N, d, M, seed) {
      calls <<- c(calls, seed)
      fake_mat
    },
    .is_d_disjunct = function(mat, d) {
      verify_calls <<- verify_calls + 1L
      verify_calls >= 2L
    }
  )

  mat <- DisjunctMatrix(N = 6L, d = 2L, M = 4L, method = "random", seed = 41L, verify = TRUE)

  expect_equal(calls, 41:42)
  expect_true(isTRUE(as.logical(attr(mat, "disjunct_verified"))))
  expect_equal(attr(mat, "seed"), 42L)
  expect_equal(attr(mat, "method_used"), "random")
})

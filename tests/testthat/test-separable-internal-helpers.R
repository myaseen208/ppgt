test_that("prime and prime-power helpers classify boundary cases exactly", {
  expect_false(.sm_is_prime(1L))
  expect_true(.sm_is_prime(2L))
  expect_true(.sm_is_prime(3L))
  expect_false(.sm_is_prime(9L))

  expect_false(.sm_is_prime_power(1L))
  expect_true(.sm_is_prime_power(2L))
  expect_true(.sm_is_prime_power(27L))
  expect_false(.sm_is_prime_power(12L))
  expect_false(.sm_is_prime_power(18L))

  expect_equal(.sm_next_prime_power(0), 2L)
  expect_equal(.sm_next_prime_power(1.2), 2L)
  expect_equal(.sm_next_prime_power(6.1), 7L)
  expect_equal(.sm_next_prime_power(8), 8L)
})

test_that("random disjunct builder respects shape and Bernoulli edge probability", {
  mat <- .build_random_disjunct_matrix(N = 4L, d = 1L, J = 3L, seed = 99L)

  expect_s4_class(mat, "dgCMatrix")
  expect_equal(dim(mat), c(3L, 4L))
  expect_true(all(as.matrix(mat) == 1L))
})

test_that(".is_d_separable covers exact positive and exact collision cases", {
  M_sep <- Matrix::Matrix(c(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  res_true <- .is_d_separable(M_sep, 2L)
  expect_true(isTRUE(as.logical(res_true)))
  expect_equal(attr(res_true, "verified"), "exact")

  M_collision <- Matrix::Matrix(c(
    1, 1,
    0, 0
  ), nrow = 2, byrow = TRUE, sparse = TRUE)
  res_false <- .is_d_separable(M_collision, 1L)
  expect_false(isTRUE(as.logical(res_false)))
  expect_equal(attr(res_false, "verified"), "exact")
})

test_that(".is_d_disjunct covers exact positive and exact collision cases", {
  M_disj <- Matrix::Matrix(c(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  res_true <- .is_d_disjunct(M_disj, 1L)
  expect_true(isTRUE(as.logical(res_true)))
  expect_equal(attr(res_true, "verified"), "exact")

  M_cover <- Matrix::Matrix(c(
    1, 1,
    0, 0
  ), nrow = 2, byrow = TRUE, sparse = TRUE)
  res_false <- .is_d_disjunct(M_cover, 1L)
  expect_false(isTRUE(as.logical(res_false)))
  expect_equal(attr(res_false, "verified"), "exact")
})

test_that("separability checker covers approximate warning and early false return", {
  M <- Matrix::Matrix(c(
    1, 1, 0,
    0, 0, 1
  ), nrow = 2, byrow = TRUE, sparse = TRUE)
  draw_idx <- 0L

  local_mocked_bindings(.EXACT_THRESHOLD = 0L)
  local_mocked_bindings(
    seq_len = function(n) 1L,
    sample.int = function(n, size) {
      if (n == 1L) {
        return(1L)
      }
      draw_idx <<- draw_idx + 1L
      if (draw_idx == 1L) 1L else 2L
    },
    .package = "base"
  )

  expect_warning(
    res <- .is_d_separable(M, 1L),
    "Monte Carlo approximation"
  )
  expect_false(isTRUE(as.logical(res)))
  expect_equal(attr(res, "verified"), "approximate")
})

test_that("separability checker covers identical-draw skip and approximate true return", {
  M <- Matrix::Matrix(c(
    1, 0,
    0, 1
  ), nrow = 2, byrow = TRUE, sparse = TRUE)

  local_mocked_bindings(.EXACT_THRESHOLD = 0L)
  local_mocked_bindings(
    seq_len = function(n) 1L,
    sample.int = function(n, size) 1L,
    .package = "base"
  )

  expect_warning(
    res <- .is_d_separable(M, 1L),
    "Monte Carlo approximation"
  )
  expect_true(isTRUE(as.logical(res)))
  expect_equal(attr(res, "verified"), "approximate")
})

test_that("disjunctness checker covers approximate warning and early false return", {
  M <- Matrix::Matrix(c(
    1, 1, 0,
    0, 0, 1
  ), nrow = 2, byrow = TRUE, sparse = TRUE)

  local_mocked_bindings(.EXACT_THRESHOLD = 0L)
  local_mocked_bindings(
    seq_len = function(n) 1L,
    sample.int = function(n, size) 1L,
    sample = function(x, size) x[[1L]],
    .package = "base"
  )

  expect_warning(
    res <- .is_d_disjunct(M, 1L),
    "Monte Carlo approximation"
  )
  expect_false(isTRUE(as.logical(res)))
  expect_equal(attr(res, "verified"), "approximate")
})

test_that("disjunctness checker covers approximate true return", {
  M <- Matrix::Matrix(c(
    1, 0,
    0, 1
  ), nrow = 2, byrow = TRUE, sparse = TRUE)

  local_mocked_bindings(.EXACT_THRESHOLD = 0L)
  local_mocked_bindings(
    seq_len = function(n) 1L,
    sample.int = function(n, size) 1L,
    sample = function(x, size) x[[length(x)]],
    .package = "base"
  )

  expect_warning(
    res <- .is_d_disjunct(M, 1L),
    "Monte Carlo approximation"
  )
  expect_true(isTRUE(as.logical(res)))
  expect_equal(attr(res, "verified"), "approximate")
})

test_that("checker helpers validate invalid inputs", {
  empty_mat <- Matrix::Matrix(integer(0), nrow = 1, ncol = 0, sparse = TRUE)
  one_col <- Matrix::Matrix(c(1L, 0L), nrow = 2, ncol = 1, sparse = TRUE)

  expect_error(.is_d_separable(empty_mat, 0L), "d must be a positive integer")
  expect_error(.is_d_separable(empty_mat, 1L), "at least one column")

  expect_error(.is_d_disjunct(one_col, 0L), "d must be a positive integer")
  expect_error(.is_d_disjunct(one_col, 1L), "more columns than d")
})

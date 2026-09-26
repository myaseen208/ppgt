test_that("regular_pooling_matrix uses pp_matrix when exact PP construction exists", {
  A <- regular_pooling_matrix(
    n_pools = 188, n_samples = 282, pool_size = 6,
    pools_per_sample = 4, seed = 4L
  )
  B <- pp_matrix(q = 47, d = 3, nl = 4, N = 282)

  expect_equal(dim(A), c(188, 282))
  expect_true(all(Matrix::rowSums(A) == 6))
  expect_true(all(Matrix::colSums(A) == 4))
  expect_false(is.null(attr(A, "pp")))
  expect_equal(attr(A, "pp")$q, 47)
  expect_equal(attr(A, "pp")$d, 3)
  expect_equal(attr(A, "pp")$nl, 4)
  expect_equal(attr(A, "pp")$N, 282)
  expect_equal(as.matrix(A), as.matrix(B))
})

test_that("regular_pooling_matrix falls back for 94 x 188", {
  A <- regular_pooling_matrix(94, 188, 8, 4, seed = 4L)
  expect_equal(dim(A), c(94, 188))
  expect_true(all(Matrix::rowSums(A) == 8))
  expect_true(all(Matrix::colSums(A) == 4))
  expect_true(is.null(attr(A, "pp")))
})

test_that("regular_pooling_matrix falls back for 186 x 372", {
  A <- regular_pooling_matrix(186, 372, 8, 4, seed = 4L)
  expect_equal(dim(A), c(186, 372))
  expect_true(all(Matrix::rowSums(A) == 8))
  expect_true(all(Matrix::colSums(A) == 4))
  expect_true(is.null(attr(A, "pp")))
})

test_that("method pp errors when exact PP construction does not exist", {
  expect_error(
    regular_pooling_matrix(94, 188, 8, 4, method = "pp"),
    "No exact pp_matrix"
  )
  expect_error(
    regular_pooling_matrix(186, 372, 8, 4, method = "pp"),
    "No exact pp_matrix"
  )
})

test_that("method regular skips pp_matrix", {
  A <- regular_pooling_matrix(
    188, 282, 6, 4, seed = 4L, method = "regular"
  )
  expect_equal(dim(A), c(188, 282))
  expect_true(all(Matrix::rowSums(A) == 6))
  expect_true(all(Matrix::colSums(A) == 4))
  expect_true(is.null(attr(A, "pp")))
})

test_that("overlap-optimized fallback improves overlap quality", {
  A <- regular_pooling_matrix(
    94, 188, 8, 4, seed = 4L, method = "overlap_optimized"
  )
  score <- score_pooling_matrix(A)
  expect_equal(score$max_col_overlap, 1)
  expect_equal(score$n_col_pairs_overlap_gt_1, 0)
  expect_false(score$has_pp_attr)
  expect_equal(unname(score[["row_sum_table"]]), unname(table(rep(8, 94))))
  expect_equal(unname(score[["col_sum_table"]]), unname(table(rep(4, 188))))
})

test_that("regular_pooling_matrix is reproducible without changing RNG state", {
  set.seed(123L)
  state_before <- .Random.seed
  A1 <- regular_pooling_matrix(12L, 18L, 6L, 4L, seed = 9L)
  expect_identical(.Random.seed, state_before)
  A2 <- regular_pooling_matrix(12L, 18L, 6L, 4L, seed = 9L)
  expect_identical(A1, A2)
})

test_that("regular_pooling_matrix rejects incompatible inputs", {
  expect_error(
    regular_pooling_matrix(94L, 188L, 7L, 4L),
    "Incompatible marginals"
  )
  expect_error(
    regular_pooling_matrix(3L, 4L, 4L, 4L),
    "cannot exceed `n_pools`"
  )
  expect_error(
    regular_pooling_matrix(4L, 3L, 4L, 3L),
    "cannot exceed `n_samples`"
  )
  expect_error(
    regular_pooling_matrix(12, 18, 6, 4, max_tries = 0),
    "max_tries"
  )
})


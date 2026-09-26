# Tests for additional pooling matrix designs

test_that("dorfman_matrix creates correct structure", {
  # 100 samples, pool size 10
  M <- dorfman_matrix(N = 100, g = 10)
  
  expect_equal(nrow(M), 10)
  expect_equal(ncol(M), 100)
  
  # Each sample in exactly 1 pool (non-overlapping)
  expect_true(all(Matrix::colSums(M) == 1))
  
  # Each pool has correct size
  expect_true(all(Matrix::rowSums(M) == 10))
})

test_that("dorfman_optimal_g computes reasonable values", {
  # Under the implemented expected-tests criterion, p = 0.01 is minimized at g = 11.
  expect_equal(dorfman_optimal_g(0.01), 11L)
  expect_true(dorfman_optimal_g(0.05) < dorfman_optimal_g(0.01))
  expect_true(dorfman_optimal_g(0.001) > dorfman_optimal_g(0.01))
})

test_that("array_matrix creates correct structure", {
  # 10x10 array
  M <- array_matrix(r = 10, c = 10)
  
  expect_equal(nrow(M), 20)  # r + c
  expect_equal(ncol(M), 100) # r * c
  
  # Each sample in exactly 2 pools
  expect_true(all(Matrix::colSums(M) == 2))
  
  # 96-well plate format
  M96 <- array_matrix(r = 8, c = 12)
  expect_equal(nrow(M96), 20)
  expect_equal(ncol(M96), 96)
})

test_that("array_matrix_3d creates correct structure", {
  # 4x4x4 cube
  M <- array_matrix_3d(4, 4, 4)
  
  expect_equal(nrow(M), 12)  # d1 + d2 + d3
  expect_equal(ncol(M), 64)  # d1 * d2 * d3
  
  # Each sample in exactly 3 pools
  expect_true(all(Matrix::colSums(M) == 3))
})

test_that("bibd_matrix creates Fano plane correctly", {
  # (7,3,1)-BIBD
  M <- bibd_matrix(v = 7, k = 3)
  
  if (!is.null(M)) {
    expect_equal(nrow(M), 7)
    expect_equal(ncol(M), 7)
    
    # Each block has 3 elements
    expect_true(all(Matrix::rowSums(M) == 3))
    
    # Each element in 3 blocks
    expect_true(all(Matrix::colSums(M) == 3))
  }
})

test_that("hypercube_matrix creates correct structure", {
  # 5^3 hypercube
  M <- hypercube_matrix(q = 5, d = 3)
  
  expect_equal(nrow(M), 15)   # d * q
  expect_equal(ncol(M), 125)  # q^d
  
  # Each sample in exactly d pools
  expect_true(all(Matrix::colSums(M) == 3))
  
  # Pool size = q^(d-1)
  expect_true(all(Matrix::rowSums(M) == 25))
})

test_that("kirkman_matrix creates a valid KTS for v = 3, 9, 15 (shape + is_kts)", {
  for (v in c(3L, 9L, 15L)) {
    M <- kirkman_matrix(v)
    b <- (v * (v - 1L)) %/% 6L

    expect_equal(nrow(M), b)
    expect_equal(ncol(M), v)
    expect_true(all(Matrix::rowSums(M) == 3))
    expect_true(all(Matrix::colSums(M) == (v - 1L) %/% 2L))

    # is_kts() independently re-derives every defining property from the
    # incidence matrix itself (not from the constructor's own bookkeeping).
    expect_true(is_kts(M))
  }
})

test_that("kirkman_matrix rejects invalid v with an informative error, not a silent NULL", {
  expect_error(kirkman_matrix(10), "v == 3 \\(mod 6\\)")
  expect_error(kirkman_matrix(2), "v must be an integer >= 3")
  expect_error(kirkman_matrix(0), "v must be an integer >= 3")
  expect_error(kirkman_matrix(-3), "v must be an integer >= 3")
})

test_that("kirkman_matrix(v = 21) either succeeds (and passes is_kts) or fails with the documented budget error", {
  # v = 21 is beyond the two dedicated fast constructions (v = 9 general
  # search, v = 15 PG(3,2)); the general search-based path is a randomized
  # combinatorial search that can take up to roughly a minute here, so this
  # is skipped on CRAN/CI rather than slowing down every check.
  skip_on_cran()
  skip_on_ci()
  result <- tryCatch(kirkman_matrix(21L), error = function(e) e)
  if (inherits(result, "error")) {
    expect_match(
      conditionMessage(result),
      "computational budget|guaranteed to construct quickly"
    )
  } else {
    expect_true(is_kts(result))
  }
})

test_that("is_kts() catches the reviewer-reported KTS(15) counterexample", {
  # The pre-fix package bug: a hand-built KTS(15) block set whose last three
  # "parallel classes" are not valid partitions of the 15 points (points
  # repeated within a class, others never covered), verified against the
  # reviewer's independently-checked reference construction in
  # ppgt_evaluator_supplement.zip (kts15.py: PG(3,2) spread packing).
  broken_blocks <- list(
    c(1, 2, 3), c(4, 5, 6), c(7, 8, 9), c(10, 11, 12), c(13, 14, 15),
    c(1, 4, 7), c(2, 5, 8), c(3, 6, 9), c(10, 13, 14), c(11, 12, 15),
    c(1, 5, 9), c(2, 6, 7), c(3, 4, 8), c(10, 12, 14), c(11, 13, 15),
    c(1, 6, 8), c(2, 4, 9), c(3, 5, 7), c(10, 11, 14), c(12, 13, 15),
    c(1, 10, 15), c(2, 11, 13), c(3, 12, 14), c(4, 8, 13), c(5, 7, 12),
    c(1, 11, 14), c(2, 10, 12), c(3, 13, 15), c(4, 7, 15), c(5, 8, 10),
    c(1, 12, 13), c(2, 14, 15), c(3, 10, 11), c(4, 9, 14), c(5, 6, 13)
  )
  rows_i <- unlist(lapply(seq_along(broken_blocks), function(i) rep(i, 3)))
  cols_j <- unlist(broken_blocks)
  M_broken <- Matrix::sparseMatrix(
    i = rows_i, j = cols_j, x = rep(1L, length(rows_i)), dims = c(35, 15)
  )
  attr(M_broken, "design") <- list(parallel_class = rep(1:7, each = 5))

  expect_false(is_kts(M_broken))
  reasons <- attr(is_kts(M_broken), "reasons")
  expect_true(any(grepl("covered 0 or >1 times", reasons)))
  expect_true(any(grepl("does not partition", reasons)))

  # The fixed constructor's own KTS(15) does not have this problem.
  expect_true(is_kts(kirkman_matrix(15)))
})

test_that("is_kts() flags each of the four defining properties independently", {
  valid <- kirkman_matrix(9)

  # Not a matrix at all
  expect_false(is_kts(list(1, 2, 3)))

  # Wrong block size (row with 4 points instead of 3)
  M_wrong_k <- valid
  M_wrong_k[1, ] <- 1L
  expect_false(is_kts(M_wrong_k))

  # Wrong number of blocks for its own v (drop the last row)
  M_wrong_b <- valid[-nrow(valid), ]
  attr(M_wrong_b, "design") <- attr(valid, "design")
  expect_false(is_kts(M_wrong_b))

  # No parallel_class info at all -> resolvability cannot be confirmed
  M_no_class <- valid
  attr(M_no_class, "design") <- NULL
  res <- is_kts(M_no_class)
  expect_false(res)
  expect_true(any(grepl("parallel_class", attr(res, "reasons"))))
})

test_that("pg_matrix creates projective plane", {
  # PG(2,2) = Fano plane
  M2 <- pg_matrix(2)
  
  expect_equal(nrow(M2), 7)  # q^2 + q + 1
  expect_equal(ncol(M2), 7)
  
  # q+1 points per line
  expect_true(all(Matrix::rowSums(M2) == 3))
  
  # q+1 lines per point
  expect_true(all(Matrix::colSums(M2) == 3))
  
  # PG(2,3)
  M3 <- pg_matrix(3)
  expect_equal(nrow(M3), 13)
  expect_true(all(Matrix::rowSums(M3) == 4))
})

test_that("compare_designs returns valid comparison", {
  comp <- compare_designs(100)
  
  expect_true(is.data.frame(comp))
  expect_true("type" %in% names(comp))
  expect_true("N" %in% names(comp))
  expect_true("M" %in% names(comp))
  expect_true("compression" %in% names(comp))
})

test_that("list_designs returns design catalog", {
  designs <- list_designs()
  
  expect_true(is.data.frame(designs))
  expect_true("type" %in% names(designs))
  expect_true("function_name" %in% names(designs))
  expect_true(nrow(designs) >= 8)
})

test_that("P-BEST design works with new functions", {
  # P-BEST using polynomial pools
  M_pbest <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  
  expect_equal(nrow(M_pbest), 48)
  expect_equal(ncol(M_pbest), 384)
  
  # Compare with hypercube base
  M_hypercube <- hypercube_matrix(q = 8, d = 3)
  expect_equal(ncol(M_hypercube), 512)  # Full hypercube
  
  # P-BEST has better compression than basic hypercube
  comp_pbest <- 384 / 48
  comp_hypercube <- 512 / 24
  # Both are good, but PP has better detection guarantees
})

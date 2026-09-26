#' @title Construct a Pooling Design Matrix by Family
#'
#' @description
#' Dispatch to one of the existing pooling-design constructors in \pkg{ppgt}
#' using a user-facing family label. The function is a convenience wrapper over
#' the package's structured constructors and returns only the pooling matrix
#' \eqn{\mathbf{M}}. It does not claim that \eqn{d} alone uniquely determines a
#' pooling matrix across design families.
#'
#' @param family Character scalar selecting the constructor family. Supported
#'   values are \code{"pp"}, \code{"pbest"}, \code{"hyper"},
#'   \code{"separable"}, \code{"disjunct"}, \code{"dorfman"},
#'   \code{"array"}, \code{"bibd"}, \code{"hypercube"},
#'   \code{"kirkman"}, and \code{"pg"}.
#' @param d Optional integer scalar. This argument is passed only to families
#'   for which \eqn{d} is a native parameter, such as Polynomial Pools,
#'   separable designs, disjunct designs, and hypercube designs. For other
#'   families it is ignored unless explicitly documented for the selected
#'   constructor.
#' @param N Optional integer scalar giving the number of individuals. When the
#'   selected family uses \eqn{N} directly, this value is forwarded to the
#'   underlying constructor. For design families whose native interface uses a
#'   different symbol such as \eqn{v}, \code{PoolMatrix()} uses \code{N} only as
#'   a convenience alias when that mapping is exact.
#' @param q Optional integer scalar. This argument is forwarded only to families
#'   where \eqn{q} is a native parameter, such as Polynomial Pools, HYPER,
#'   hypercube designs, and projective-plane designs.
#' @param J Optional integer scalar giving the requested number of pools. This
#'   is forwarded only to families where the native constructor accepts a free
#'   row-count argument. For Polynomial Pools, \code{J} may be used only to
#'   infer \eqn{n_l = J/q} when that ratio is an integer.
#' @param method Optional character scalar giving a family-specific construction
#'   mode. Examples include \code{"pp"} or \code{"clinical"} for
#'   \code{family = "pbest"}, \code{"auto"}, \code{"random"}, or
#'   \code{"reed-solomon"} for \code{family = "separable"} and
#'   \code{family = "disjunct"}, and \code{"2d"} or \code{"3d"} for
#'   \code{family = "array"}.
#' @param verify Logical scalar. This is forwarded to families that expose an
#'   explicit verification option, namely \code{"separable"} and
#'   \code{"disjunct"}.
#' @param seed Optional integer scalar. This is forwarded only to families whose
#'   current constructor uses randomness, such as random separable/disjunct
#'   designs and the package-specific clinical P-BEST helper.
#' @param ... Additional family-specific arguments. Supported names are:
#'   \describe{
#'     \item{\code{nl}}{Number of PP layers for \code{family = "pp"} or the
#'       PP-based \code{"pbest"} path.}
#'     \item{\code{reorder}}{Logical reorder flag for \code{family = "hyper"}.}
#'     \item{\code{g}}{Pool size for \code{family = "dorfman"}.}
#'     \item{\code{r}, \code{c}}{Grid dimensions for the two-dimensional array
#'       path.}
#'     \item{\code{d1}, \code{d2}, \code{d3}}{Grid dimensions for the
#'       three-dimensional array path.}
#'     \item{\code{v}, \code{k}, \code{lambda}}{BIBD or Kirkman parameters.}
#'   }
#' @param verbose Logical scalar; if \code{TRUE}, emit a short message
#'   describing the chosen dispatch path.
#'
#' @return A binary pooling matrix returned by the selected constructor, usually
#'   of class \code{dgCMatrix}. The matrix is returned directly, and carries an
#'   additional \code{"poolmatrix"} attribute containing the selected family,
#'   dispatch method, and the key arguments used by the wrapper. The return
#'   value is the stage-1 design matrix \eqn{\mathbf{M}} only.
#'
#' @details
#' Let \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} denote a pooling design matrix,
#' where \eqn{J} is the number of pools and \eqn{N} is the number of
#' individuals. The wrapper
#' \code{PoolMatrix()} is a constructor-only dispatcher over existing matrix
#' constructors. It returns \eqn{\mathbf{M}} and stops there. In particular, it
#' does not handle adaptive versus non-adaptive workflow selection, retesting
#' scenarios, stage-2 protocol logic, or decoder behavior. Those concerns belong
#' in separate functions or in package documentation, not in this wrapper.
#'
#' The function also does not solve a universal inverse problem of the form
#' "given \eqn{d}, construct the unique matrix." The symbol \eqn{d} has
#' different meanings across the supported families, and some families do not
#' use \eqn{d} at all.
#'
#' The dispatch rules are:
#' \describe{
#'   \item{\code{family = "pp"}}{Dispatch to \code{\link{pp_matrix}} using
#'   \eqn{q}, \eqn{d}, \eqn{n_l}, and optional \eqn{N}.}
#'   \item{\code{family = "pbest"}}{Dispatch either to the standard PP-based
#'   48-pool family through \code{\link{pp_matrix}} with
#'   \eqn{(q,d,n_l,N) = (8,3,6,384)}, or to the package-specific clinical helper
#'   \code{\link{pbest_clinical_matrix}} when \code{method = "clinical"}.}
#'   \item{\code{family = "hyper"}}{Dispatch to \code{\link{hyper_matrix}} using
#'   \eqn{N}, \eqn{J}, and \eqn{q}.}
#'   \item{\code{family = "separable"}}{Dispatch to
#'   \code{\link{SeparableMatrix}} using \eqn{N}, \eqn{d}, optional \eqn{J},
#'   and the selected method.}
#'   \item{\code{family = "disjunct"}}{Dispatch to
#'   \code{\link{DisjunctMatrix}} using \eqn{N}, \eqn{d}, optional \eqn{J}, and
#'   the selected method.}
#'   \item{\code{family = "dorfman"}}{Dispatch to
#'   \code{\link{dorfman_matrix}} using \eqn{N} and the explicit pool-size
#'   argument \code{g}.}
#'   \item{\code{family = "array"}}{Dispatch to \code{\link{array_matrix}} or
#'   \code{\link{array_matrix_3d}} according to \code{method}.}
#'   \item{\code{family = "bibd"}}{Dispatch to \code{\link{bibd_matrix}} using
#'   \eqn{v}, \eqn{k}, and \eqn{\lambda}.}
#'   \item{\code{family = "hypercube"}}{Dispatch to
#'   \code{\link{hypercube_matrix}} using \eqn{q} and \eqn{d}.}
#'   \item{\code{family = "kirkman"}}{Dispatch to
#'   \code{\link{kirkman_matrix}} using \eqn{v}.}
#'   \item{\code{family = "pg"}}{Dispatch to \code{\link{pg_matrix}} using
#'   \eqn{q}.}
#' }
#'
#' Family-specific parameters remain explicit. In particular, the wrapper does
#' not impose \eqn{d}-only semantics on families such as HYPER, Dorfman, array,
#' BIBD, Kirkman, or projective-plane designs, where \eqn{d} is not the native
#' defining parameter. Likewise, calling \code{PoolMatrix()} should be read only
#' as "construct \eqn{\mathbf{M}} for the chosen family" and not as "run the
#' complete testing workflow" for that family.
#'
#' @examples
#' M_pp <- PoolMatrix(family = "pp", q = 4, d = 3, nl = 5, N = 64)
#' M_pp
#'
#' M_pbest <- PoolMatrix(family = "pbest")
#' M_pbest
#'
#' M_hyper <- PoolMatrix(family = "hyper", N = 12, J = 6, q = 2, reorder = FALSE)
#' M_hyper
#'
#' M_sep <- PoolMatrix(
#'   family = "separable",
#'   N = 16L,
#'   d = 2L,
#'   method = "reed-solomon",
#'   verify = TRUE
#' )
#' M_sep
#' attr(M_sep, "separable_verified")
#'
#' M_dorfman <- PoolMatrix(family = "dorfman", N = 12L, g = 3L)
#' M_dorfman
#'
#' M_array <- PoolMatrix(family = "array", method = "2d", r = 4L, c = 5L)
#' M_array
#'
#' M_bibd <- PoolMatrix(family = "bibd", N = 7L, k = 3L)
#' M_bibd
#'
#' M_cube <- PoolMatrix(family = "hypercube", q = 2L, d = 3L)
#' M_cube
#'
#' M_pg <- PoolMatrix(family = "pg", q = 2L)
#' M_pg
#'
#' # PoolMatrix() constructs the matrix only.
#' # Decoding and downstream testing decisions are handled elsewhere.
#'
#' @references
#' Tan, Y. H. I. (2020). Pooling matrix designs for group testing.
#' \emph{SIAM Undergraduate Research Online}, 13, 1-21.
#'
#' Shental, N., Levy, S., Wuvshet, V., et al. (2020). Efficient high-throughput
#' SARS-CoV-2 testing to detect asymptomatic carriers. \emph{Science Advances},
#' 6(37), eabc5961.
#'
#' Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing for
#' SARS-CoV-2: A flexible, easy-to-implement, and highly efficient method.
#' \emph{Nature Communications}, 13, Article 3626.
#'
#' D'yachkov, A. G., & Rykov, V. V. (1982). Bounds on the length of disjunctive
#' codes. *Problemy Peredachi Informatsii*, 18(3), 7-13.
#'
#' Dorfman, R. (1943). The detection of defective members of large populations.
#' \emph{The Annals of Mathematical Statistics}, 14(4), 436-440.
#'
#' @seealso
#' \code{\link{pp_matrix}}, \code{\link{pbest_clinical_matrix}},
#' \code{\link{hyper_matrix}}, \code{\link{SeparableMatrix}},
#' \code{\link{DisjunctMatrix}}, \code{\link{dorfman_matrix}},
#' \code{\link{array_matrix}}, \code{\link{bibd_matrix}},
#' \code{\link{hypercube_matrix}}, \code{\link{kirkman_matrix}},
#' \code{\link{pg_matrix}}
#'
#' @family pooling_designs
#' @export
PoolMatrix <- function(family = c(
                         "pp", "pbest", "hyper", "separable", "disjunct",
                         "dorfman", "array", "bibd", "hypercube",
                         "kirkman", "pg"
                       ),
                       d = NULL,
                       N = NULL,
                       q = NULL,
                       J = NULL,
                       method = NULL,
                       verify = TRUE,
                       seed = NULL,
                       ...,
                       verbose = FALSE) {
  family <- match.arg(family)
  dots <- list(...)

  .dot_value <- function(name, default = NULL) {
    if (name %in% names(dots)) dots[[name]] else default
  }

  .require_not_null <- function(x, name, family_name) {
    if (is.null(x)) {
      stop(
        "PoolMatrix(family = \"", family_name, "\") requires `", name, "`."
      )
    }
    x
  }

  .validate_equals <- function(actual, expected, label) {
    if (!is.null(actual) && !identical(as.integer(actual), as.integer(expected))) {
      stop(label, " must equal ", expected, " for this dispatch path.")
    }
  }

  .attach_meta <- function(M, family_name, method_name, info) {
    if (is.null(M)) {
      return(NULL)
    }
    attr(M, "poolmatrix") <- list(
      family = family_name,
      method = method_name,
      dispatch = info
    )
    M
  }

  if (isTRUE(verbose)) {
    message("PoolMatrix dispatching family = ", family)
  }

  if (family == "pp") {
    q <- as.integer(.require_not_null(q, "q", family))
    d <- as.integer(.require_not_null(d, "d", family))
    nl <- .dot_value("nl")
    if (is.null(nl)) {
      if (is.null(J)) {
        stop(
          "PoolMatrix(family = \"pp\") requires `nl`, or `J` together with `q` ",
          "so that `nl = J / q` is defined."
        )
      }
      if (J %% q != 0L) {
        stop("For family = \"pp\", `J` must be divisible by `q` to infer `nl`.")
      }
      nl <- J %/% q
    }
    out <- pp_matrix(
      q = q,
      d = d,
      nl = as.integer(nl),
      N = if (is.null(N)) NULL else as.integer(N)
    )
    return(.attach_meta(out, family, "pp_matrix", list(q = q, d = d, nl = nl, N = N)))
  }

  if (family == "pbest") {
    method_name <- if (is.null(method)) "pp" else match.arg(method, c("pp", "clinical"))
    if (method_name == "clinical") {
      out <- pbest_clinical_matrix(
        N = if (is.null(N)) 384L else as.integer(N),
        J = if (is.null(J)) 94L else as.integer(J),
        seed = if (is.null(seed)) 42L else as.integer(seed)
      )
      return(.attach_meta(
        out,
        family,
        method_name,
        list(N = N, J = J, seed = seed)
      ))
    }

    q <- if (is.null(q)) 8L else as.integer(q)
    d <- if (is.null(d)) 3L else as.integer(d)
    nl <- .dot_value("nl", 6L)
    out <- pp_matrix(
      q = q,
      d = d,
      nl = as.integer(nl),
      N = if (is.null(N)) 384L else as.integer(N)
    )
    return(.attach_meta(
      out,
      family,
      method_name,
      list(q = q, d = d, nl = nl, N = N)
    ))
  }

  if (family == "hyper") {
    N <- as.integer(.require_not_null(N, "N", family))
    J <- as.integer(.require_not_null(J, "J", family))
    q <- if (is.null(q)) 2L else as.integer(q)
    reorder <- isTRUE(.dot_value("reorder", TRUE))
    out <- hyper_matrix(n = N, m = J, q = q, reorder = reorder)
    return(.attach_meta(
      out,
      family,
      "hyper_matrix",
      list(N = N, J = J, q = q, reorder = reorder)
    ))
  }

  if (family == "separable") {
    N <- as.integer(.require_not_null(N, "N", family))
    d <- as.integer(.require_not_null(d, "d", family))
    method_name <- if (is.null(method)) "auto" else match.arg(
      method,
      c("auto", "random", "reed-solomon")
    )
    out <- SeparableMatrix(
      N = N,
      d = d,
      M = if (is.null(J)) NULL else as.integer(J),
      method = method_name,
      seed = if (is.null(seed)) 42L else as.integer(seed),
      verify = isTRUE(verify)
    )
    return(.attach_meta(
      out,
      family,
      method_name,
      list(N = N, d = d, J = J, verify = verify, seed = seed)
    ))
  }

  if (family == "disjunct") {
    N <- as.integer(.require_not_null(N, "N", family))
    d <- as.integer(.require_not_null(d, "d", family))
    method_name <- if (is.null(method)) "auto" else match.arg(
      method,
      c("auto", "random", "reed-solomon")
    )
    out <- DisjunctMatrix(
      N = N,
      d = d,
      M = if (is.null(J)) NULL else as.integer(J),
      method = method_name,
      seed = if (is.null(seed)) 42L else as.integer(seed),
      verify = isTRUE(verify)
    )
    return(.attach_meta(
      out,
      family,
      method_name,
      list(N = N, d = d, J = J, verify = verify, seed = seed)
    ))
  }

  if (family == "dorfman") {
    N <- as.integer(.require_not_null(N, "N", family))
    g <- as.integer(.require_not_null(.dot_value("g"), "g", family))
    out <- dorfman_matrix(N = N, g = g)
    return(.attach_meta(out, family, "dorfman_matrix", list(N = N, g = g)))
  }

  if (family == "array") {
    method_name <- if (is.null(method)) "2d" else match.arg(method, c("2d", "3d"))
    if (method_name == "2d") {
      r <- as.integer(.require_not_null(.dot_value("r"), "r", family))
      c <- as.integer(.require_not_null(.dot_value("c"), "c", family))
      out <- array_matrix(r = r, c = c)
      return(.attach_meta(
        out,
        family,
        method_name,
        list(r = r, c = c)
      ))
    }

    d1 <- as.integer(.require_not_null(.dot_value("d1"), "d1", family))
    d2 <- as.integer(.require_not_null(.dot_value("d2"), "d2", family))
    d3 <- as.integer(.require_not_null(.dot_value("d3"), "d3", family))
    out <- array_matrix_3d(d1 = d1, d2 = d2, d3 = d3)
    return(.attach_meta(
      out,
      family,
      method_name,
      list(d1 = d1, d2 = d2, d3 = d3)
    ))
  }

  if (family == "bibd") {
    v <- .dot_value("v")
    if (is.null(v)) {
      v <- .require_not_null(N, "v or N", family)
    }
    k <- as.integer(.require_not_null(.dot_value("k"), "k", family))
    lambda <- as.integer(.dot_value("lambda", 1L))
    out <- bibd_matrix(v = as.integer(v), k = k, lambda = lambda)
    return(.attach_meta(
      out,
      family,
      "bibd_matrix",
      list(v = v, k = k, lambda = lambda)
    ))
  }

  if (family == "hypercube") {
    q <- as.integer(.require_not_null(q, "q", family))
    d <- as.integer(.require_not_null(d, "d", family))
    if (!is.null(N)) {
      .validate_equals(N, q^d, "N")
    }
    out <- hypercube_matrix(q = q, d = d)
    return(.attach_meta(out, family, "hypercube_matrix", list(q = q, d = d)))
  }

  if (family == "kirkman") {
    v <- .dot_value("v")
    if (is.null(v)) {
      v <- .require_not_null(N, "v or N", family)
    }
    out <- kirkman_matrix(v = as.integer(v))
    return(.attach_meta(out, family, "kirkman_matrix", list(v = v)))
  }

  if (family == "pg") {
    q <- as.integer(.require_not_null(q, "q", family))
    implied_n <- q^2 + q + 1L
    if (!is.null(N)) {
      .validate_equals(N, implied_n, "N")
    }
    if (!is.null(J)) {
      .validate_equals(J, implied_n, "J")
    }
    out <- pg_matrix(q = q)
    return(.attach_meta(out, family, "pg_matrix", list(q = q)))
  }

  stop("Unsupported family: ", family)
}

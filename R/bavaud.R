#' @import magrittr
#' @import parallel
#' @importFrom stats dist lm
#' @importFrom utils combn
#' @importFrom matrixcalc vech
#' @importFrom RcppAlgos permuteGeneral comboGeneral
NULL

chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

if (nzchar(chk) && chk == "TRUE") {
  # use 2 cores in CRAN/Travis/AppVeyor
  num_workers <- 2L
} else {
  # use all cores in devtools::test()
  num_workers <- parallel::detectCores()
}

#' Compute validity matrix
#'
#' @param M a matrix containing NAs
#' @return The validity matrix containing 1 if the element is a NA 0 else.
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' validity_matrix(M)
#' @export
validity_matrix <- function(M) {
  matrix(as.numeric(!is.na(M)), nrow = nrow(M)) %>%
    `rownames<-`(rownames(M)) %>%
    `colnames<-`(colnames(M))
}

#' Compute validity weight
#'
#' @param M a matrix containing NAs
#' @return a named array containing weights
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' validity_weight(M)
#' @export
validity_weight <- function(M) {
  row_sums <- rowSums(
    validity_matrix(M)
  )
  result <- row_sums / sum(row_sums)
  names(result) <- rownames(M)
  result
}

#' Compute the L1 dissimilarity
#'
#' @param M a matrix containing NAs
#' @return a dissimilarity renormalised by disputedness
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' dissimilarity_L1(M)
#' @export
dissimilarity_L1 <- function(M) {
  n <- nrow(M)
  diss <- function(ij) {
    i <- ij[1]
    j <- ij[2]
    x <- M[i, ]
    y <- M[j, ]
    k <- which(!is.na(x + y))
    ifelse(
      is.null(k),
      NA,
      sum(abs(x[k] - y[k])) / length(k)
    )
  }

  matrix(unlist(RcppAlgos::permuteGeneral(
    c(seq(n)),2, FUN = diss, repetition = TRUE, Parallel = TRUE)),
  nrow = n)
}

#' Compute column disputedness
#'
#' @param M a matrix containing NAs
#' @param f an array containing weights
#' @return a named array containing the column disputedness
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' f <- rowSums(M, na.rm = TRUE) / sum(M, na.rm = TRUE)
#' disputedness(M,f)
#' @export
disputedness <- function(M, f = NULL) {
  n <- nrow(M)
  if (is.null(f)) f <- rep(1 / n, n)
  p <- ncol(M)
  ks <- seq(p)
  out <- lapply(ks, function(k) {
    result <- list()
    not_na <- which(!is.na(M[, k]))
    result$numk <- sum(unlist(
      RcppAlgos::comboGeneral(
        not_na, 2,
        FUN = function(d) {
          f[d[1]] * f[d[2]] * abs(M[d[1], k] - M[d[2], k])
        },
        Parallel = TRUE
      ))
    )
    result$denk <- sum(unlist(
      RcppAlgos::comboGeneral(not_na, 2, FUN = function(d) {f[d[1]] * f[d[2]]},
                              Parallel = TRUE)
    )
    )

    result
  })
  out <- do.call(rbind.data.frame, out)

  result <- out$numk / out$denk
  names(result) <- colnames(M)
  result
}


#' Compute the disputedness dissimilarity
#'
#' @param M a matrix containing NAs
#' @param f a named weight array
#' @param disp an array containing column disputedness
#' @return a dissimilarity renormalised by disputedness
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' f <- rowSums(M, na.rm = TRUE) / sum(M, na.rm = TRUE)
#' disp <- disputedness(M)
#' dissimilarity_disputedness(M, f,disp)
#' @export
dissimilarity_disputedness <- function(M, f, disp) {
  n <- nrow(M)
  if (is.null(weights)) {
    weights <- rep(1 / n, n)
  }
  diss <- function(ij) {
    i <- ij[1]
    j <- ij[2]
    x <- M[i, ]
    y <- M[j, ]
    k <- which(!is.na(x + y))
    ifelse(
      is.null(k),
      NA,
      sum(
        f[i] * f[j] * disp[k] * abs(x[k] - y[k])
      ) / sum(disp[k])
    )
  }

    dr <- matrix(unlist(RcppAlgos::permuteGeneral(
    c(seq(n)), 2, FUN = diss, repetition = TRUE, Parallel = TRUE)),
    nrow = n)

  rownames(dr) <- rownames(M)
  colnames(dr) <- rownames(M)
  dr
}

#' Compute the dissimilarity estimation
#'
#' @param D a dissimilarity
#' @param f a named weight array
#' @return a dissimilarity estimation
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' f <- rowSums(M, na.rm = TRUE) / sum(M, na.rm = TRUE)
#' disp <- disputedness(M)
#' D <- dissimilarity_disputedness(M, f, disp)
#' dissimilarity_tilde_estimation(D, f)
#' @export
dissimilarity_tilde_estimation <- function(D, f = NULL) {
  n <- nrow(D)
  if (is.null(f)) {
    f <- rep(1 / n, n)
  }
  diss <- function(ij) {
    i <- ij[1]
    j <- ij[2]
    x <- D[i, ]
    y <- D[j, ]
    k <- which(!is.na(x + y))
    ifelse(
      is.null(k),
      NA,
      sum(
        f[k] * (x[k] - y[k])^2
      ) / sum(f[k]) - (sum(f[k] * (x[k] - y[k])) / sum(f[k]))^2
    )
  }
  D_tilde <- matrix(unlist(RcppAlgos::permuteGeneral(
    c(seq(n)),2, FUN = diss, repetition = TRUE, Parallel = TRUE)),
    nrow = n)

  rownames(D_tilde) <- rownames(D)
  colnames(D_tilde) <- colnames(D)
  D_tilde
}

#' Compute the dissimilarity estimation
#'
#' @param D a dissimilarity
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' f <- rowSums(M, na.rm = TRUE) / sum(M, na.rm = TRUE)
#' d <- disputedness(M)
#' D <- dissimilarity_disputedness(M, f, d)
#' D_estim <- dissimilarity_regression_estimation(D)
#' @export
dissimilarity_regression_estimation <- function(D) {
  dm <- D
  diag(dm) <- NA
  y <- vech(sqrt(D))
  dv <- vech(dm)
  A <- stats::lm(y ~ dv - 1)
  D_estimation <- sqrt(D) / A$coefficients
  rownames(D_estimation) <- rownames(D)
  colnames(D_estimation) <- colnames(D)
  D_estimation
}

#' Compute the dissimilarity estimation
#'
#' @param D_estim an estimated dissimilarity
#' @param D a dissimilarity
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' f <- rowSums(M, na.rm = TRUE) / sum(M, na.rm = TRUE)
#' d <- disputedness(M)
#' D <- dissimilarity_disputedness(M, f, d)
#' D_estim <- dissimilarity_regression_estimation(D)
#' @export
dissimilarity_final <- function(D_estim, D) {
  n <- nrow(D)
  diss <- function(ij) {
    i <- ij[1]
    j <- ij[2]
    ifelse(
      is.na(D[i, j]),
      D_estim[i, j],
      0.5 * (D_estim[i, j] + D[i, j])
    )
  }

  D_final <- matrix(unlist(RcppAlgos::permuteGeneral(
    c(seq(n)),2, FUN = diss, repetition = TRUE, Parallel = TRUE)),
    nrow = n)


  rownames(D_final) <- rownames(D)
  colnames(D_final) <- colnames(D)
  D_final
}

#' Compute the distance estimation
#'
#' @param M a matrix
#' @param f a named weight array
#' @param disp disputedness
#' @return a dissimilarity estimation
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' f <- rowSums(M, na.rm = TRUE) / sum(M, na.rm = TRUE)
#' d <- disputedness(M)
#' D <- estimate_distance(M, f, d)
#' @export
estimate_distance <- function(M, f, disp = NULL) {
  if (is.null(disp)) disp <- disputedness(M, f)
  D_disp <- dissimilarity_disputedness(M, f, disp)
  dtilde <- dissimilarity_tilde_estimation(D_disp, f)
  D_estim <- dissimilarity_regression_estimation(dtilde)
  dissimilarity_final(D_estim, D_disp)
}

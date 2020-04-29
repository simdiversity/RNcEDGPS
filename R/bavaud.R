#' @import magrittr
#' @importFrom stats dist lm
#' @importFrom utils combn
#' @importFrom matrixcalc vech
#' @importFrom parallel mcmapply
#'
NULL
#' Compute validity matrix
#'
#' @param M a matrix containing NAs
#' @return The validity matrix containing 1 if the element is a NA 0 else.
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' validity.matrix(M)
#' @export
validity.matrix <- function(M) {
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
#' validity.weight(M)
#' @export
validity.weight <- function(M) {
    row_sums <- rowSums(
      validity.matrix(M)
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
#' f <- rowSums(M, na.rm = TRUE) / sum(M, na.rm = TRUE)
#' dissimilarity.L1(M)
#' @export
dissimilarity.L1 <- function(M){
  n <- nrow(M)
  diss <- function(i, j) {
    x <- M[i,]
    y <- M[j,]
    k <- which(!is.na(x + y))
    ifelse(
      is.null(k),
      NA,
      sum(abs(x[k] - y[k]))/length(k)
    )
  }
  opts <- options(mccores = parallel::detectCores())
  matrix(
    do.call(
      parallel::mcmapply,
      c(diss, unname(expand.grid(seq(n), seq(n))))
    ),
    nrow = n
  )
  options(opts)
}

#' Compute column disputedness
#'
#' @param M a matrix containing NAs
#' @param f an array containing weights
#' @return a named array containing the column disputedness
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' disputedness(M)
#' @export
disputedness <- function(M, f = NULL) {
  n <- nrow(M)
  if (is.null(f)) f <- rep(1/n, n)
  p <- ncol(M)
  ks <- seq(p)
  opts <- options(mccores = parallel::detectCores())
  out <- parallel::mclapply(ks, function(k) {
    result <- list()
    not_na <- which(!is.na(M[, k]))
    result$numk <- sum(
      combn(
        not_na, 2,
        function(d) {
          f[d[1]] * f[d[2]] * abs(M[d[1], k] - M[d[2], k])
        }
      )
    )
    result$denk <- sum(
      combn(not_na, 2, function(d) f[d[1]] * f[d[2]])
    )
    result
  })
  options(opts)
  out <- do.call(rbind.data.frame, out)
  result <- out$numk / out$denk
  names(result) <- colnames(M)
  result
}

#' Compute the disputedness dissimilarity
#'
#' @param M a matrix containing NAs
#' @param f a named weight array
#' @param disputedness an array containing column disputedness
#' @return a dissimilarity renormalised by disputedness
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' f <- rowSums(M, na.rm = TRUE) / sum(M, na.rm = TRUE)
#' d <- disputedness(M)
#' dissimilarity.disputedness(M, f, d)
#' @export
dissimilarity.disputedness <- function(M, f, disputedness) {
  n <- nrow(M)
  if (is.null(weights)) {
    weights <- rep(1 / n, n)
  }
  diss <- function(i, j) {
    x <- M[i, ]
    y <- M[j, ]
    k <- which(!is.na(x + y))
    ifelse(
      is.null(k),
      NA,
      sum(f[i] * f[j] * disputedness[k] * abs(x[k] - y[k])) / sum(disputedness[k])
    )
  }
  opts <- options(mccores = parallel::detectCores())
  dr <- matrix(
    do.call(
      parallel::mcmapply,
      c(diss, unname(expand.grid(seq(n), seq(n))))
    ),
    nrow = n
  )
  options(opts)
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
#' d <- disputedness(M)
#' D <- dissimilarity.disputedness(M, f, d)
#' dissimilarity.tilde_estimation(D, f)
#' @export
dissimilarity.tilde_estimation <- function(D, f = NULL) {
  n <- nrow(D)
  if (is.null(f)) {
    f <- rep(1 / n, n)
  }
  diss <- function(i, j) {
    x <- D[i, ]
    y <- D[j, ]
    k <- which(!is.na(x + y))
    ifelse(
      is.null(k),
      NA,
      sum(f[k] * (x[k] - y[k])^2) / sum(f[k]) - (sum(f[k] * (x[k] - y[k])) / sum(f[k]))^2
    )
  }
  opts <- options(mccores = parallel::detectCores())
  D_tilde <- matrix(
    do.call(
      parallel::mcmapply,
      c(diss, unname(expand.grid(seq(n), seq(n))))
    ),
    nrow = n
  )
  options(opts)
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
#' D <- dissimilarity.disputedness(M, f, d)
#' D_estim <- dissimilarity.regression_estimation(D)
#' @export
dissimilarity.regression_estimation <- function(D) {
  dm <- D
  diag(dm) <- NA
  y <- vech(sqrt(D))
  dv <- vech(dm)
  A <- lm(y ~ dv - 1)
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
#' D <- dissimilarity.disputedness(M, f, d)
#' D_estim <- dissimilarity.regression_estimation(D)
#' @export
dissimilarity.final <- function(D_estim, D) {
  n <- nrow(D)
  diss <- function(i, j) {
    ifelse(
      is.na(D[i, j]),
      D_estim[i, j],
      0.5 * (D_estim[i, j] + D[i, j])
    )
  }
  opts <- options(mccores = parallel::detectCores())
  D_final <- matrix(
    do.call(
      parallel::mcmapply,
      c(diss, unname(expand.grid(seq(n), seq(n))))
    ),
    nrow = n
  )
  options(opts)
  rownames(D_final) <- rownames(D)
  colnames(D_final) <- colnames(D)
  D_final
}

#' Compute the distance estimation
#'
#' @param M a matrix
#' @param f a named weight array
#' @return a dissimilarity estimation
#'
#' @examples
#' M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
#' f <- rowSums(M, na.rm = TRUE) / sum(M, na.rm = TRUE)
#' D <- estimate.distance(M, f)
#' @export
estimate.distance <- function(M, f) {
  disp <- disputedness(M,f)
  D_disp <- dissimilarity.disputedness(M, f, disp)
  dtilde <- dissimilarity.tilde_estimation(D_disp, f)
  D_estim <- dissimilarity.regression_estimation(dtilde)
  dissimilarity.final(D_estim, D_disp)
}


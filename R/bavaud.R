library(parallel)
library(matrixcalc)

weight_valid <- function(in_data){
  valid <- validity_matrix(in_data)
  result = rowSums(valid)/sum(valid)
  names(result) <- rownames(valid)
  result
}

validity_matrix <- function(in_data){
  n <- dim(in_data)
  p <- n[2]
  n <- n[1]
  valid <- matrix(as.numeric(!is.na(in_data)), n, p)
  rownames(valid) <- rownames(in_data)
  colnames(valid) <- colnames(in_data)
  valid
}


disputedness <- function(in_data, weights=c()) {
  n = nrow(in_data)
  p = ncol(in_data)
  numk = rep(0,p)
  denk = rep(0,p)
  out <- lapply(seq(p), function(k) {
    result = list()
    not_na <- which(!is.na(in_data[,k]), useNames = FALSE)
    result$numk <- sum(combn(not_na, 2, function(d) weights[d[1]] * weights[d[2]] * abs(in_data[d[1], k] - in_data[d[2], k])))
    result$denk <- sum(combn(not_na, 2, function(d) weights[d[1]] * weights[d[2]]))
    result
  })
  print(out)
  out = do.call(rbind.data.frame, out )
  result = out$numk/out$denk
  names(result) <- colnames(in_data)
  result
}


dissimilarity.renormalised <- function(in_data, disputedness) {
  n = nrow(in_data)
  if (is.null(weights)) {
    weights = rep(1/n, n)
  }
  diss <- function(x, y) {
    k <- which(!is.na(x + y))
    ifelse(
      is.null(k),
      NA,
      sum(disputedness[k] * abs( x[k] - y[k]))/sum(disputedness[k])
    )
  }
  dr = sapply(1:n, function(i) sapply(1:n, function(j) diss(in_data[i,], in_data[j,])))
  rownames(dr) <- rownames(in_data)
  colnames(dr) <- rownames(in_data)
  dr
}


dissimilarity.tilde_estimation = function(dissimilarity, weights) {
  n = nrow(dissimilarity)
  if (is.null(weights)) {
    weights = rep(1/n, n)
  }
  diss <- function(x, y) {
    k <- which(!is.na(x + y))
    ifelse(
      is.null(k),
      NA,
      sum(weights[k] * (x[k] - y[k])^2)/sum(weights[k]) - (sum(weights[k] * (x[k] - y[k]))/sum(weights[k]))^2
    )
  }
  dissimilarity_tilde = sapply(1:n, function(i) sapply(1:n, function(j) diss(dissimilarity[i,], dissimilarity[j,])))
  rownames(dissimilarity_tilde) <- rownames(dissimilarity)
  colnames(dissimilarity_tilde) <- colnames(dissimilarity)
  dissimilarity_tilde
}


dissimilarity.regression_estimation <- function(dissimilarity, q) {
  dm = dissimilarity
  diag(dm) <- NA
  y <- vech(sqrt(dissimilarity))
  dv <- vech(dm)
  A = lm( y ~ dv - 1)
  dissimilarity_estimation = sqrt(dissimilarity) / A$coefficients
  rownames(dissimilarity_estimation) <- rownames(dissimilarity)
  colnames(dissimilarity_estimation) <- colnames(dissimilarity)
  dissimilarity_estimation
}


dissimilarity.final <- function(d_estim, dissimilarity) {
  n = nrow(dissimilarity)
  dfinal=sapply(1:n, function(i)
    sapply(1:n, function(j)
      ifelse(
        is.na(dissimilarity[i,j]),
        d_estim[i, j],
        0.5 * (d_estim[i, j] + dissimilarity[i, j])
      )
    )
  )
  rownames(dfinal) <- rownames(dissimilarity)
  colnames(dfinal) <- colnames(dissimilarity)
  dfinal
}

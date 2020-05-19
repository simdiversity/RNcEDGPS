test_that("disputedness works", {
  M <- matrix(c(1, 2, NA, 3, NA ,5, NA, 6, 5, 4, 19, 5, 0, NA, 0, 3, 4, 5, 6, 8, 8, 5,6,6), nrow = 4)
  f <- c(.4, .25, .2, .1, .05)
  p <- ncol(M)
  n <- nrow(M)
  do <- c()
  NumK = rep(0L,p);
  DenK = rep(0L,p)

  for (k in 1:p) {
  for (i in 1:n) {
    for (j in 1:n) {
      if ((is.na(M[i,k]) == FALSE) & (is.na(M[j,k]) == FALSE)) {
        NumK[k] = NumK[k] + f[i] * f[j] * abs(M[i,k] - M[j,k])
        DenK[k] = DenK[k] + f[i] * f[j]}
      }}
  do[k] = NumK[k]/DenK[k]
  }
  dc <-  disputedness(M, f)
  expect(all.equal(do,dc), "There is a problem computing disputedness")
})

test_that("disputedness works", {
  M <- matrix(c(1, 2, NA, NA, 4, 19, 0, NA, 0), nrow = 3)
  f <- c(.5, .3, .2)
  p <- ncol(M)
  n <- nrow(M)
  do <- c()
  for (k in 1:p) {
  NumK = 0L; DenK = 0L
  for (i in 1:n) {
    for (j in 1:n) {
      if ((is.na(M[i,k]) == FALSE) & (is.na(M[j,k]) == FALSE)) {
        NumK = NumK + f[i] * f[j] * abs(M[i,k] - M[j,k])
        DenK = DenK + f[i] * f[j]}
      }}
  do[k] = NumK/DenK
  }
  expect(all.equal(do, disputedness(M, f)))
})

#' importFrom BB BBsolve
NULL
# # CONVERTED FROM:
# # Abiy Tasissa and Rongjie Lai, "Exact Reconstruction of Euclidean Distance Geometry Problem Using Low-rank Matrix Completion," in IEEE Transactions on Information Theory, 2018. doi: 10.1109/TIT.2018.2881749
# # https://github.com/abiy-tasissa/Nonconvex-Euclidean-Distance-Geometry-Problem-Solver
# #
# # TODO: needs testing!
#
# # ORIGINAL NOTE:
# # --------------------------------------------------------------------------
# # This script is a numerical solution of the Euclidean distance geometry
# # problem (EDG). Formally, consider a set of n points where partial
# # inter-point distance information is provided. The goal of EDG is to find
# # the coordinate of the points given this partial information.
# # --------------------------------------------------------------------------
# # The n points usually lie in a low dimensional space of size r << n, low
# # rank. With this, the EDG problem can be set as low-rank completion problem
# # which can be solved via nuclear norm minimization. We recover the Gram
# # matrix, the inner product matrix, and follow classical MDS to recover the
# # coordinates. For details, see the associated paper below.
# # --------------------------------------------------------------------------
# # Tasissa, Abiy, and Rongjie Lai."Exact Reconstruction of Euclidean Distance
# # Geometry Problem Using Low-rank Matrix Completion." arXiv preprint
# # arXiv:1804.04310 (2018).
# # --------------------------------------------------------------------------
# # The gram matrix is psd so nuclear norm minimization equates to trace.
# # Our algorithm uses the Augmented Lagrangian framework to find the Gram
# # matrix.
# # --------------------------------------------------------------------------
# # We assume the partial inter-distance information comes from a uniform
# # random sample. The information is assumed to be exact.
# # Dist is the full distance matrix. Weight is a binary matrix informing
# # whether a given entry of Dist is chosen or not.
# # --------------------------------------------------------------------------
# # Rongjie Lai, Abiy Tasissa
# # --------------------------------------------------------------------------
#
alternating_completion <- function(Dist,Selected,opts,lsopts) {
  norm <- function(x) sqrt(sum(x^2))
  # -------------------------------------------------------------------------
  # Constructs the operator A which captures the linear operator:
  # R_{\omega}(X) <- R_{\omega}(M)
  # -------------------------------------------------------------------------
  A_operator <- function(X, i,j) {
    X <- rowSums(X[i,] * X[j,])
      diag(X) <- rowSums(X*X)
      X_offdiag <- rowSums(X[i,]*X[j,]);
      Y <- diag(X[i,]) + diag(X[j,]) - 2 * X_offdiag
      Y
  }
  # -------------------------------------------------------------------------
  # Constructs the adjoint operator for A (see associated paper for details).
  # -------------------------------------------------------------------------
  At_operator <- function(y, edgeind) {
    n <- nrow(y)
    O <- matrix(0, n,n)
    O[edgeind,] <- -2*y[seq(n)]
    diag(O) <- -rowSums(O)
    O
  }
  # -------------------------------------------------------------------------
  # a function handle for the line search BB algorithm
  # -------------------------------------------------------------------------
  gradient <- function(P, b, D1, r, i, j, edgeind) {
    tmp <- A_operator(P, i, j) - b + D1
    F <- sum(rowSums(P*P)) + 0.5 * r * norm(tmp)^2;
    G <- 2 * P + 2 * r * At_operator(tmp, edgeind) * P
    cbind(F , G)
  }
  # aug. lagrangian penalty and estimate of the rank
  r <- opts$r;
  Rk <- opts$rank;
  # n <- number of points
  # Dist <- D^{2}(i,j) <- {d_{i,j}^{2}} matrix
  n <- nrow(D)
  # calculate the ground truth of inner-product matrix
  IPM_Truth <- Dist - rowMeans(Dist)*diag(1, n);
  IPM_Truth <- -.5 * (IPM_Truth - array(1, c(n,1)) * colMeans(IPM_Truth));
  # indices of the randomly chosen entries of D
  IJ <- which(Selected == 1);
  i <- IJ[1]
  j <- IJ[2]
  edgeind <- i + (j - 1)*n
  # the minimization problem has linear constraint R_{Omega}(X) <- R_{Omega)(M)
  # Represent this constraint as A(X) <- b where b <- R_{omega}(M)
  M <- Dist(edgeind);
  b <- M;
  # -------------------------------------------------------------------------
  # main algorithm
  # -------------------------------------------------------------------------
  # initialize P, lagrangian multiplier
  P <- runif(n,Rk)
  R <- P
  D1 <- rep(length(b),1)
  # initialize energies
  E1 <- rep(0, opts$maxit)
  E <- rep(0, opts$maxit)
  num_it <- 0
  cre <- 1
  # -------------------------------------------------------------------------
    # main iteration: BB method to solve for P
  # -------------------------------------------------------------------------
  for (i in seq(opts$maxit)) {
    num_it <- num_it + 1;
    # do line search based gradient descent for P
    P <- BB::BBsolve(P,gradient(P, b, D1, r, i, j, edgeind)$P,lsopts);
    # update multiplier D
    tmperr <- A_operator(P, i, j) - b
    D1 <- D1 + tmperr
    # total energy
    E(i) <- sum(rowSums(P*P)) + 0.5*r*norm(D1)^2
    if (opts$printenergy == 1) {
      sprintf('Iteration %d, TotalE <- %f\n',i,E(i))
    }
    # calculate the energy E1
    E1(i) <- 0.5 * r * norm(tmperr)^2
    # stopping condition
    if (i > 2) {
      cre <- abs(E[i] - E[i - 1])/E[i]
    }
    if (E1(i) < opts$tol && cre < opts$tol) {
      break
    }
  }
  IPM <- P*t(P)
  IPM_Recon <- IPM - (1/n)*matrix(rowSums(IPM),1,n) - (1/n)*matrix(colSums(IPM),n,1) + (1/(n*n))*matrix(sum(IPM),n,n)
  IPM_Recon <- IPM_Recon + t(IPM_Recon)/2
  IPM_err <- norm(IPM_Truth - IPM_Recon)/norm(IPM_Truth)
  output <- list()
  output$E1 <- E1[seq(num_it)]
  output$E <- E[num_it]
  output$ReconError <- IPM_err
  output$numit <- num_it
  output
}

#' sphere_vertices
#'
#' @param xc x cooddinate for the center of the sphere
#' @param yc y cooddinate for the center of the sphere
#' @param zc z cooddinate for the center of the sphere
#' @param r radius
#' @param lats number of latitudes
#' @param longs number of longitudes
#'
#' @export
sphere_vertices <- function(xc = 0,
                            yc = 0,
                            zc = 0,
                            r = 1,
                            lats = 16L,
                            longs = 16L
                            ) {
  # xc,yc,zc give centre of sphere, r is radius, lats/longs for resolution
  vertices <- vector(mode = "numeric", length = 12L * lats * longs)
  vi <- 1L
  for (i in 1:lats) {
    lat0 <- pi * (-0.5 + (i - 1) / lats)
    z0 <- sin(lat0) * r
    zr0 <- cos(lat0) * r
    lat1 <- pi * (-0.5 + i / lats)
    z1 <- sin(lat1) * r
    zr1 <- cos(lat1) * r
    for (j in 1:longs) {
      lng1 <- 2 * pi * (j - 1) / longs
      lng2 <- 2 * pi * (j) / longs
      x1 <- cos(lng1)
      y1 <- sin(lng1)
      x2 <- cos(lng2)
      y2 <- sin(lng2)
      vertices[vi] <- x1 * zr0 + xc
      vi <- vi + 1L
      vertices[vi] <- y1 * zr0 + yc
      vi <- vi + 1L
      vertices[vi] <- z0 + zc
      vi <- vi + 1L
      vertices[vi] <- x1 * zr1 + xc
      vi <- vi + 1L
      vertices[vi] <- y1 * zr1 + yc
      vi <- vi + 1L
      vertices[vi] <- z1 + zc
      vi <- vi + 1L
      vertices[vi] <- x2 * zr1 + xc
      vi <- vi + 1L
      vertices[vi] <- y2 * zr1 + yc
      vi <- vi + 1L
      vertices[vi] <- z1 + zc
      vi <- vi + 1L
      vertices[vi] <- x2 * zr0 + xc
      vi <- vi + 1L
      vertices[vi] <- y2 * zr0 + yc
      vi <- vi + 1L
      vertices[vi] <- z0 + zc
      vi <- vi + 1L
    }
  }
  #indices <- 1:(length(vertices) / 3)
  matrix(vertices, ncol = 3)
}

#' NA_dist_random
#' @description
#'
#' Randomly sets a percentage of elements of a symmetric matrix to NA. diag being excluded from replacements.
#'
#' @param M a symmetric matrix
#' @param percentage the percentage of elements to replace by NA
#'
#' @export

NA_dist_random <- function(M, percentage) {
  n <- sum(upper.tri(M))
  rm <- sample(n, ceiling(n*percentage))
  M[upper.tri(M)] <- replace(M[upper.tri(M)], rm, NA_real_)
  M[lower.tri(M)] <- t(M)[lower.tri(t(M))]
  M
}

evaluation <- environment()
data(list = c("eurodist"), envir = evaluation)
eurodist_m <- as.matrix(eurodist)
sphere <- sphere_vertices()
spheredist <-  as.matrix(dist(sphere))
eurodist_NA <- NA_dist_random(eurodist_m, .2)

spheredist_NA <-  NA_dist_random(spheredist, .2)

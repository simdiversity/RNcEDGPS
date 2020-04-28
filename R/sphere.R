
sphere_vertices_distance <- function() {

  sphere_vertices <- function(xc = 0,
                              yc = 0,
                              zc = 0,
                              r = 1,
                              lats = 16L,
                              longs = 16L,
                              ...
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
    indices <- 1:(length(vertices) / 3)
    matrix(vertices, ncol = 3)
  }

  sphere <- sphere_vertices()
  dist(sphere)
}

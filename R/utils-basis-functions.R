rbf.distance <- function(x, c) {
  ## x points, c centers
  x2 <- apply(x^2, 2, sum)
  c2 <- apply(c^2, 2, sum)
  outer(x2, c2, "+") - 2 * t(x) %*% c ##  was 0.5 * ( )
}

kernel.Gaussian <- function(x, c, sigma) {
  ## x points, c centers, sigma kernel sd width
  distance2 <- rbf.distance(x, c)
  ## no need to scale by constant if the same for all basis fcts:
  exp(-distance2 / sigma^2) # /(2*sigma)^n.dim,  was (2*sigma^2)
}

neighbor.kernel.Gaussian <- function(x, c, c.dists, sigma, n.dim) {
  ## turn cols with centers into rows for multiplication
  distance2 <- t(rbf.distance(x, c))
  distance2[1:5]
  ## easy to multiply now
  ## ret.t <- exp(-distance2/(2*c.dists^2*sigma^2))/(2*pi*c.dists^2 * sigma^2)^(n.dim/2)
  ret.t <- exp(-distance2 / (c.dists * sigma)^2) / c.dists^n.dim
  t(ret.t)
}

local.distances <- function(c, neigh.rank = 5) {
  c.dist <- rbf.distance(c, c)
  c.dist[c.dist < 0] <- 0
  c.dist <- c.dist^(1 / 2)
  apply(c.dist, 1, function(d) d[rank(d) >= neigh.rank][1])
}

moderate.distances <- function(d, low, high) {
  r <- range(d)
  low + (d - r[1]) * (high - low) / (r[2] - r[1])
}

lin.down <- function(x, t) {
  pmax(-as.numeric(x) + t, 0)
}

lin.up <- function(x, t) {
  pmax(as.numeric(x) - t, 0)
}

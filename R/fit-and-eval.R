#' Evaluate the density ratio estimate
#' 
#' Evaluate the density ratio estimate at \code{x.grid}
#' 
#' @param basis.fit a list with elements documented in \code{\link{clsif}}
#' @param x.grid a matrix of points to evaluate the density ratio at. Note that
#' this should have d rows (where d is the dimension of the density/ratio space)
#' and each point should be a separate column.
#' 
#' @export
eval.basis <- function(basis.fit, x.grid) {
  if (!is.matrix(x.grid)) x.grid <- matrix(x.grid, nrow = 1) # Thu 29 Nov 15:14:29 2018
  if (basis.fit$is.adaptive) {
    X.grid <- neighbor.kernel.Gaussian(
      x = x.grid, c = basis.fit$x.ce,
      c.dists = basis.fit$c.dists,
      sigma = basis.fit$sigma,
      n.dim = nrow(x.grid)
    )
  } else {
    X.grid <- kernel.Gaussian(x.grid, basis.fit$x.ce, basis.fit$sigma)
  }
  y <- as.numeric(X.grid %*% basis.fit$alpha)
  range(x.grid)
  range(y)
  as.numeric(X.grid %*% basis.fit$alpha)
}

#' Estimate the parameters for the basis functions
#' 
#' This function actually does all of the work, the \code{\link{clsif}} and 
#' \code{\link{ulsif}} functions just call this with their respective learning
#' functions.
fit.dr <- function(x.de, x.nu, lambda, sigma.chosen,
                   is.adaptive = FALSE, neigh.rank,
                   kernel.low, kernel.high,
                   b, fold, learning.fct) {
  if (!is.matrix(x.de)) x.de <- matrix(x.de, nrow = 1)
  if (!is.matrix(x.nu)) x.nu <- matrix(x.nu, nrow = 1)

  if (nrow(x.de) != nrow(x.nu)) {
    cat("x.de, x.nu unequal dimension\n")
  }
  n.nu <- ncol(x.nu)

  ## Choosing Gaussian kernel center `x_ce', at most b
  rand.index <- sample(1:n.nu)
  b <- min(b, n.nu)
  x.ce <- x.nu
  x.ce <- x.nu[, rand.index[1:b], drop = FALSE]
  dim(x.ce)

  if (sigma.chosen <= 0) {
    sigma <- 4 # very sensitive to initial sigma
    score <- -Inf
    eps <- log10(sigma) - 1
    for (eps in seq(log10(sigma) - 1, -5, -0.5)) { # feel like we should make more bigger steps
      cat("new eps", 10^eps, "\n")
      for (round in 1:9) { # why 9 ?
        sigma.new <- sigma - 10^eps
        cat("current sigma", sigma, "\n")
        cv.index <- sample(1:n.nu) # RANDOM split in CV sets
        cv.split <- floor((0:(n.nu - 1)) * fold / n.nu) + 1
        score.new <- 0

        X.de <- kernel.Gaussian(x = x.de, c = x.ce, sigma.new)
        X.nu <- kernel.Gaussian(x.nu, x.ce, sigma.new)
        for (i in 1:fold) {
          alpha.cv <- learning.fct(X.de, X.nu[cv.index[cv.split != i], ], lambda = lambda)$alpha # lambda was forgotten
          wh.cv <- X.nu[cv.index[cv.split == i], ] %*% alpha.cv
          score.new <- score.new + mean(log(wh.cv)) / fold
        }
        if (score.new <= score) break
        score <- score.new
        sigma <- sigma.new
        cat("improved score", score, " at sigma", sigma, "\n")
      }
      sigma.chosen <- sigma
      cat("final sigma", sigma.chosen, "\n")
    }
  }

  ## assuming sigma.chosen
  if (is.adaptive) {
    c.dists <- local.distances(c = x.ce, neigh.rank = neigh.rank)
    c.dists <- moderate.distances(d = c.dists, kernel.low, kernel.high)
    X.de <- neighbor.kernel.Gaussian(x = x.de, c = x.ce, c.dists, sigma.chosen, nrow(x.ce))
    X.nu <- neighbor.kernel.Gaussian(x.nu, x.ce, c.dists, sigma.chosen, nrow(x.ce))
  } else {
    c.dists <- rep(1, ncol(x.ce))
    X.de <- kernel.Gaussian(x = x.de, c = x.ce, sigma.chosen)
    X.nu <- kernel.Gaussian(x.nu, x.ce, sigma.chosen)
  }

  fit <- learning.fct(X.de, X.nu, lambda)

  list(
    alpha = fit$alpha, score = fit$score, x.ce = x.ce, sigma = sigma.chosen,
    is.adaptive = is.adaptive, c.dists = c.dists
  )
}
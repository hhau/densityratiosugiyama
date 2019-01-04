clsif.learning <- function(X.de, X.nu, lambda, max.iteration = 100,
                           eps.list = 10^seq(3, -3, -1)) {
  ## compare Sugiyama book algorithm p 68

  mean.X.de <- apply(X.de, 2, mean)
  ss.mean.X.de <- sum(mean.X.de^2) # psi_bar_de' psi_bar_de

  dim(X.de)
  H <- t(X.de) %*% X.de / nrow(X.de) # basis.fcts x basis.fcts matrix
  h <- apply(X.nu, 2, mean)
  # alpha <- rep(0.01,ncol(X.de))
  # score <- 0.5*t(alpha) %*% H %*% alpha - sum(h*alpha) + lambda*sum(alpha)
  alpha <- rep(0.01, ncol(X.de))
  score <- Inf

  i <- 1
  eps <- eps.list[1]
  for (eps in eps.list) {
    for (i in 1:max.iteration) {
      ## min 0.5*alpha H alpha - h alpha + lambda * 1vec * alpha, alpha >= 0
      ## so direction is -grad
      alpha.new <- alpha + eps * (-H %*% alpha + h - lambda)
      alpha.new <- alpha.new + ((1 - sum(mean.X.de * alpha.new)) / ss.mean.X.de) * mean.X.de
      alpha.new <- pmax(0, alpha.new)
      alpha.new <- alpha.new / sum(mean.X.de * alpha.new)
      score.new <- (0.5 * t(alpha.new) %*% H %*% alpha.new - sum(h * alpha.new)
        + lambda * sum(alpha.new))
      if (score.new >= score) break # no improvement any more
      score <- score.new
      alpha <- alpha.new
    }
  }

  list(alpha = alpha, score = score)
}

#' Constrained Least Squares Importance Fitting (CLISF)
#'
#' Sugiyama density ratio estimation method, with an L1 penalty on the 
#' parameters.
#' 
#' \code{x.de} and \code{x.nu} should be the same dimension (same number of
#' rows), but there can an uneven number of samples (number of rows)
#' 
#' 
#' @param x.de A matrix with d rows, with one sample from p(x_de) per column.
#' @param x.nu A matrix with d rows, with one sample from p(x_nu) per column.
#' @param lambda Positive real number. Regularisation parameter, see Sugiyama,
#'  Suzuki and Kanamori (2012) Section 6.2.1 for details
#' @param sigma.chosen Positive real number. Sigma for the Gaussian kernel 
#' radial basis functions. If this is set to zero, will be chosen via cross 
#' validation.
#' @param is.adaptive Boolean. Adaptively choose location of basis functions.
#' @param neigh.rank Positive integer. How many other kernels to use to compute
#' distance metrics.
#' @param kernel.low Real number. Lower bound for rescaled distances. 
#' @param kernel.high Real number. Upper bound for rescaled distances.
#' @param b Positive integer. How many kernels to use.
#' @param fold Positive integer. How many cross validation folds to use to 
#' select \code{sigma.chosen}
#'
#' @return list with the following elements:
#' \describe{
#'   \item{alpha}{basis function parameter estimates.}
#'   \item{score}{final cross validation score, used to select sigma.chosen.}
#'   \item{x.ce}{the chosen centers for the density ratio.}
#'   \item{sigma}{the value of sigma.chosen after the cross validation.} 
#'   \item{is.adaptive}{the value of is.adaptive - used to figure out which 
#'     basis function to call later.}
#'   \item{c.dists}{vector of distances between centers, used if is.adaptive 
#'     is true}
#' }
#' Note that this is list is meant to be passed to \code{\link{eval.basis}}. It
#' also serves as a small way to represent the estimated density ratio.
#' @export
clsif <- function(x.de, x.nu, lambda, sigma.chosen = 0.2, is.adaptive = FALSE,
                  neigh.rank = 5, kernel.low = 0.5, kernel.high = 2, b = 50, fold = 6) {
  fit.dr(x.de, x.nu,
    lambda = lambda, sigma.chosen = sigma.chosen,
    is.adaptive = is.adaptive, neigh.rank = neigh.rank,
    kernel.low = kernel.low, kernel.high = kernel.high,
    b = b, fold = fold, learning.fct = clsif.learning
  )
}
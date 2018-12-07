ulsif.learning <- function(X.de, X.nu, lambda, max.iteration = 100,
                           eps.list = 10^seq(3, -3, -1)) {
  ## with kernels this amounts to integrating out weight parameters for RBFs
  ## with a prior lambda I
  ## Since it has the form
  H <- t(X.de) %*% X.de / nrow(X.de) # basis.fcts x basis.fcts matrix
  h <- apply(X.nu, 2, mean)
  alpha <- solve(H + lambda * diag(nrow(H)), h)
  score <- 0.5 * t(alpha) %*% H %*% alpha - sum(h * alpha) + 0.5 * lambda * sum(alpha^2)

  list(alpha = alpha, score = score)
}
#' Unconstrained Least Square Importance Fitting (ULSIF)
#' 
#' Sugiyama, Suzuki and Kanamori (2012) density ratio estimation method
#' with L2 penalty on the basis parameters.
#' 
#' \code{x.de} and \code{x.nu} should be the same dimension (same number of
#' rows), but there can an uneven number of samples (number of rows).
#' 
#' @inheritParams clsif
#' @export
ulsif <- function(x.de, x.nu, lambda, sigma.chosen = 0.2, is.adaptive = FALSE,
                  neigh.rank = 5, kernel.low = 0.5, kernel.high = 2, b = 50, fold = 6) {
  fit.dr(x.de, x.nu,
    lambda = lambda, sigma.chosen = sigma.chosen,
    is.adaptive = is.adaptive, neigh.rank = neigh.rank,
    kernel.low = kernel.low, kernel.high = kernel.high,
    b = b, fold = fold, learning.fct = ulsif.learning
  )
}
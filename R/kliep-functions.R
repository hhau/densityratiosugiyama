kliep.projection <- function(alpha, mean.X.de, X.nu, ss.mean.X.de) {
  ## <x,a.new> = <x,a> + <x,x>*(1 - <x,a>)/<x,x> = <x,a> - 1 - <x,a> = 1
  alpha <- alpha + ((1 - sum(mean.X.de * alpha)) / ss.mean.X.de) * mean.X.de
  alpha <- pmax(0, alpha)
  alpha <- alpha / sum(mean.X.de * alpha)
  X.nu.alpha <- X.nu %*% alpha
  score <- mean(log(X.nu.alpha))
  list(alpha = alpha, X.nu.alpha = X.nu.alpha, score = score)
}

kliep.learning.proj <- function(X.de, X.nu, lambda, max.iteration = 100,
                                eps.list = 10^seq(3, -3, -1)) {
  ## compare Sugiyama book algorithm p 58
  ## mean.X.de is psi_bar_de = sum_sample_i (psi_1(x_i),psi_2(x_i),...)'
  ## X.nu = psi_nu = (psi_1(x_1), psi_2(x_2) ...; psi_1(x_2), psi_2(x_2), ... ]
  ## ie samples down rows, functions along columns

  mean.X.de <- apply(X.de, 2, mean)

  nc <- ncol(X.nu)
  dim(X.nu)
  ss.mean.X.de <- sum(mean.X.de^2) # psi_bar_de' psi_bar_de
  alpha <- rep(1, nc)

  kp <- kliep.projection(alpha, mean.X.de, X.nu, ss.mean.X.de)
  (score <- kp$score)
  alpha <- kp$alpha
  X.nu.alpha <- kp$X.nu.alpha

  i <- 1
  eps <- eps.list[1]
  for (eps in eps.list) {
    for (i in 1:max.iteration) {
      ## max 1/n.nu sum_i log (t(nu_i)*alpha), t(X.mean.de)*alpha = 1, alpha >= 1
      ## gradient d/dalpha log (t(nu_i)*alpha) = t(nu_i) * 1/t(nu_i)*alpha
      ## rows in t(X.nu) are basis fcts
      alpha.new <- alpha + eps * t(X.nu) %*% (1 / X.nu.alpha)
      kp <- kliep.projection(alpha.new, mean.X.de, X.nu, ss.mean.X.de)
      if (kp$score <= score) break
      score <- kp$score
      alpha <- kp$alpha
      X.nu.alpha <- kp$X.nu.alpha
    }
  }

  list(alpha = alpha, score = score)
}

kliep.learning <- function(X.de, X.nu, lambda, max.iteration = 100,
                           eps.list = 10^seq(3, -3, -1)) {
  ## compare Sugiyama book algorithm p 58
  ## mean.X.de is psi_bar_de = sum_sample_i (psi_1(x_i),psi_2(x_i),...)'
  ## X.nu = psi_nu = (psi_1(x_1), psi_2(x_2) ...; psi_1(x_2), psi_2(x_2), ... ]
  ## ie samples down rows, functions along columns

  mean.X.de <- apply(X.de, 2, mean)
  ss.mean.X.de <- sum(mean.X.de^2) # psi_bar_de' psi_bar_de
  alpha <- rep(10, ncol(X.nu))
  # alpha <- rep(0.1,ncol(X.nu))
  alpha <- alpha + ((1 - sum(mean.X.de * alpha)) / ss.mean.X.de) * mean.X.de
  alpha <- pmax(0, alpha)
  alpha <- alpha / sum(mean.X.de * alpha)
  score <- mean(log(X.nu %*% alpha))

  alpha <- rep(100, ncol(X.nu))
  score <- -Inf

  for (eps in eps.list) {
    for (i in 1:max.iteration) {
      ## max 1/n.nu sum_i log (t(nu_i)*alpha), t(X.mean.de)*alpha = 1, alpha >= 1
      ## gradient d/dalpha log (t(nu_i)*alpha) = t(nu_i) * 1/t(nu_i)*alpha
      ## rows in t(X.nu) are basis fcts
      alpha.new <- alpha + eps * t(X.nu) %*% (1 / (X.nu %*% alpha))
      ## <x,a.new> = <x,a> + <x,x>*(1 - <x,a>)/<x,x> = <x,a> - 1 - <x,a> = 1
      alpha.new <- alpha.new + ((1 - sum(mean.X.de * alpha.new)) / ss.mean.X.de) * mean.X.de
      alpha.new <- pmax(0, alpha.new)
      alpha.new <- alpha.new / sum(mean.X.de * alpha.new)
      score.new <- mean(log(X.nu %*% alpha.new))
      if (score.new <= score) break
      score <- score.new
      alpha <- alpha.new
    }
  }
  # alpha.new <- alpha

  list(alpha = alpha, score = score)
}

kliep <- function(x.de, x.nu, x.grid, sigma.chosen = 0, is.adaptive = FALSE,
                  neigh.rank = 5, kernel.low = 0.5, kernel.high = 2, b = 50, fold = 6) {
  ## Notes: sigma.chosen needs to be accurate close to kernel bw of x.de
  ## b must not be too large (probably ill conditioned matrices) around 50
  ## matlab code from  http://www.ms.k.u-tokyo.ac.jp/software.html#KLIEP
  ##
  ## Kullback-Leiblar importance estimation procedure (with cross validation)
  ##
  ## Estimating ratio of probability densities
  ##   \frac{ p_{nu}(x) }{ p_{de}(x) }
  ## from samples
  ##    { xde_i | xde_i\in R^{d} }_{i=1}^{n_{de}}
  ## drawn independently from p_{de}(x) and samples
  ##    { xnu_i | xnu_i\in R^{d} }_{i=1}^{n_{nu}}
  ## drawn independently from p_{nu}(x).
  ##
  ## Usage:
  ##       [wh_x_de,wh_x_re]=KLIEP(x_de,x_nu,x_re,sigma_chosen,b)
  ##
  ## Input:
  ##    x_de:         d by n_de sample matrix corresponding to `denominator' (iid from density p_de)
  ##    x_nu:         d by n_nu sample matrix corresponding to `numerator'   (iid from density p_nu)
  ##    x_re:         (OPTIONAL) d by n_re reference input matrix
  ##    sigma_chosen: (OPTIONAL) positive scalar representing Gaussian kernel width;
  ##                  if omitted (or zero), it is automatically chosen by cross validation
  ##    b:            (OPTINLAL) positive integer representing the number of kernels (default: 100);
  ##
  ## Output:
  ##    wh_x_de:      estimates of density ratio w=p_nu/p_de at x_de
  ##    wh_x_re:      estimates of density ratio w=p_nu/p_de at x_re (if x_re is provided)
  ##
  ## (c) Masashi Sugiyama, Department of Compter Science, Tokyo Institute of Technology, Japan.
  ##     sugi@cs.titech.ac.jp,     http://sugiyama-www.cs.titech.ac.jp/~sugi/software/KLIEP/

  fit.dr(x.de, x.nu, x.grid,
    lambda = 0, sigma.chosen = sigma.chosen,
    is.adaptive = is.adaptive, neigh.rank = neigh.rank,
    kernel.low = kernel.low, kernel.high = kernel.high,
    b = b, fold = fold, learning.fct = kliep.learning
  )
}

loglin.kliep.learning <- function(mean.X.nu, X.de, max.iteration = 100,
                                  eps.list = 10^seq(3, -3, -1)) {
  ## compare Sugiyama "Density ratio" book algorithm p 60, loglinear model
  ## mean.X.de is psi_bar_de = sum_sample_i (psi_1(x_i),psi_2(x_i),...)'
  ## X.nu = psi_nu = (psi_1(x_1), psi_2(x_2) ...; psi_1(x_2), psi_2(x_2), ... ]
  ## ie samples down rows, functions along columns

  n.de <- ncol(X.de)
  alpha <- rep(0.1, n.de)

  lg <- loglin.grad(alpha, mean.X.nu, X.de)
  score <- lg$score
  grad <- lg$grad

  for (eps in eps.list) {
    for (i in 1:max.iteration) {
      alpha.new <- alpha + eps * grad
      lg <- loglin.grad(alpha = alpha.new, mean.X.nu, X.de)
      if (lg$score <= score) break
      alpha <- alpha.new
      score <- lg$score
      grad <- lg$grad
    }
  }

  list(alpha = alpha, score = score)
}

loglin.kliep <- function(x.de, x.nu, x.grid, sigma.chosen = 0, b = 100) {
  if (!is.matrix(x.de)) x.de <- matrix(x.de, nrow = 1)
  if (!is.matrix(x.nu)) x.nu <- matrix(x.nu, nrow = 1)
  if (!is.matrix(x.grid)) x.grid <- matrix(x.grid, nrow = 1)

  if (nrow(x.de) != nrow(x.nu)) {
    cat("x.de, x.nu unequal dimension\n")
  }
  n.nu <- ncol(x.nu)

  ## Choosing Gaussian kernel center `x_ce', at most b
  rand.index <- sample(1:n.nu)
  b <- min(b, n.nu)
  x.ce <- x.nu
  x.ce <- x.nu[, rand.index[1:b], drop = FALSE]

  ## assuming sigma.chosen
  X.de <- kernel.Gaussian(x = x.de, c = x.ce, sigma.chosen)
  X.de <- cbind(X.de, 1, t(x.de), t(x.de^2), t(x.de^4))
  ## X.de <- cbind(X.de,lin.down(x.de,median(x.de)),lin.up(x.de,median(x.de)))
  X.nu <- kernel.Gaussian(x.nu, x.ce, sigma.chosen)
  X.nu <- cbind(X.nu, 1, t(x.nu), t(x.nu^2), t(x.nu^3))
  # X.nu <- cbind(X.nu,lin.down(x.nu,median(x.de)),lin.up(x.nu,median(x.de)))
  mean.X.nu <- apply(X.nu, 2, mean)
  alpha <- loglin.kliep.learning(mean.X.nu, X.de)$alpha

  if (!is.null(x.grid)) {
    X.grid <- kernel.Gaussian(x.grid, x.ce, sigma.chosen)
    X.grid <- cbind(X.grid, 1, t(x.grid), t(x.grid^2), t(x.grid^4))
    # X.grid <- cbind(X.grid,lin.down(x.grid,median(x.de)),lin.up(x.grid,median(x.de)))
    y.grid <- loglin.ratio(alpha, X.grid, X.de)
  }
  list(alpha = alpha, y.grid = y.grid)
}
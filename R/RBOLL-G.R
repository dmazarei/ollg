#' The Ristic-Balakrishnan  Odd log-logistic family of distributions (RBOLL-G)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the beta extended distribution due to Esmaeili et al. (2020) specified by the pdf
#' \deqn{f(x)=\frac{\alpha\,g(x)\,G(x)^{\alpha-1}\bar{G}(x)^{\alpha-1}}{\Gamma(\beta)\left[G(x)^\alpha+\bar{G}(x)^\alpha\right]^2}\,\left\{-\log\left[\frac{G(x)^\alpha}{G(x)^\alpha+\bar{G}(x)^\alpha}\right]\right\}^{\beta-1}}
#' for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf,  \eqn{\alpha > 0}, the first shape parameter, and \eqn{\beta > 0}, the second shape parameter.
#'
#' @name RBOLLG
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the value of the first shape parameter, must be positive, the default is 1.
#' @param beta the value of the second shape parameter, must be positive, the default is 1.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{prbollg} gives the distribution function,
#'  \code{drbollg} gives the density,
#'  \code{qrbollg} gives the quantile function,
#'  \code{hrbollg} gives the hazard function and
#'  \code{rrbollg} generates random variables from the The Ristic-Balakrishnan  Odd log-logistic family of
#'  distributions (RBOLL-G) for baseline cdf G.
#' @references Esmaeili, H., Lak, F., Altun, E. (2020). The Ristic-Balakrishnan odd log-logistic family of distributions: Properties and Applications. Statistics, Optimization Information Computing, 8(1), 17-35.
#' @importFrom stats numericDeriv  pnorm  rgamma qgamma pgamma uniroot  integrate
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' prbollg(x)
#' prbollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
prbollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  u <- -log(G^alpha / (G^alpha + (1 - G)^alpha))
  F0 <- pgamma(u, shape = beta) - pgamma(0, shape = beta)
  return(1 - F0)
}


#'
#' @name RBOLLG
#' @examples
#' drbollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(drbollg, -3, 3)
#' @export
drbollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * g * G^(alpha - 1) * (1 - G)^(alpha - 1) / (gamma(beta) * ((G^alpha) + (1 - G)^alpha)^2) * (-log(G^alpha / (G^alpha + (1 - G)^alpha)))^(beta - 1)
  return(df)
}


#'
#' @name RBOLLG
#' @examples
#' qrbollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qrbollg <- function(q, alpha = 1, beta = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - prbollg(t, alpha, beta, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name RBOLLG
#' @examples
#' n <- 10
#' rrbollg(n, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
rrbollg <- function(n, alpha = 1, beta = 1, G = pnorm, ...) {
  v <- runif(1e4)
  Q_G <- function(y) qrbollg(y, alpha, beta, G, ...)
  X <- Q_G(exp(-qgamma((1 - v), shape = beta) / alpha) / (exp(-qgamma(1 - v, shape = beta) / alpha) + (1 - exp(-qgamma(1 - v, shape = beta)))^(1 / alpha)))
  return(X)
}

#'
#' @name RBOLLG
#' @examples
#' hrbollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hrbollg, -3, 3)
#' @export
hrbollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  u <- -log(G^alpha / (G^alpha + (1 - G)^alpha))
  n <- length(u)
  F0 <- rep(NA, n)
  for (i in 1:n) {
    F0[i] <- integrate(function(t) 1 / gamma(beta) * t^(beta - 1) * exp(-t), 0, u[i])$value
  }
  h <- alpha * g * G^(alpha - 1) * (1 - G)^(alpha - 1) * (-log(G^alpha / (G^alpha + (1 - G)^alpha)))^(beta - 1) / (gamma(beta) * ((G^alpha) + (1 - G)^alpha)^2 * F0)
  return(h)
}

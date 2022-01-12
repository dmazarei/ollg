#' The Zografos-Balakrishnan  Odd log-logistic family of distributions (ZBOLL-G)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the beta extended distribution due to Cordeiro et al. (2016) specified by the pdf
#' \deqn{f=\frac{\alpha\,g\,G^{\alpha-1}\bar{G}^{\alpha-1}}{\Gamma(\beta)[G^\alpha+\bar{G}^\alpha]^2}\,\{-\log[1-\frac{G^\alpha}{G^\alpha+\bar{G}^\alpha}]\}^{\beta-1}}
#' for \eqn{G} any valid continuous cdf , \eqn{\bar{G}=1-G}, \eqn{g} the corresponding pdf, \eqn{\Gamma(\beta)} the Gamma funcion,  \eqn{\alpha > 0}, the first shape parameter, and \eqn{\beta > 0}, the second shape parameter.
#'
#' @name ZBOLLG
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the value of the first shape parameter, must be positive, the default is 1.
#' @param beta the value of the second shape parameter, must be positive, the default is 1.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{pzbollg} gives the distribution function,
#'  \code{dzbollg} gives the density,
#'  \code{qzbollg} gives the quantile function,
#'  \code{hzbollg} gives the hazard function and
#'  \code{rzbollg} generates random variables from the The Zografos-Balakrishnan  Odd log-logistic family of
#'  distributions (ZBOLL-G) for baseline cdf G.
#' @references Cordeiro, G. M., Alizadeh, M., Ortega, E. M., Serrano, L. H. V. (2016). The Zografos-Balakrishnan odd log-logistic family of distributions: Properties and Applications. Hacettepe Journal of Mathematics and Statistics, 45(6), 1781-1803. .
#' @importFrom stats numericDeriv  pnorm  rgamma pgamma uniroot  integrate
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' pzbollg(x)
#' pzbollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
pzbollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  u <- -log(1 - G^alpha / (G^alpha + (1 - G)^alpha))
  F0 <- pgamma(u, shape = beta) - pgamma(0, shape = beta)
  return(F0)
}


#'
#' @name ZBOLLG
#' @examples
#' dzbollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(dzbollg, -3, 3)
#' @export
dzbollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * g * G^(alpha - 1) * (1 - G)^(alpha - 1) / (gamma(beta) * ((G^alpha) + (1 - G)^alpha)^2) * (-log(1 - G^alpha / (G^alpha + (1 - G)^alpha)))^(beta - 1)
  return(df)
}


#'
#' @name ZBOLLG
#' @examples
#' qzbollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qzbollg <- function(q, alpha = 1, beta = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pzbollg(t, alpha, beta, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name ZBOLLG
#' @examples
#' n <- 10
#' rzbollg(n, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
rzbollg <- function(n, alpha = 1, beta = 1, G = pnorm, ...) {
  v <- rgamma(n, shape = beta)
  Q_G <- function(y) qzbollg(y, alpha, beta, G, ...)
  X <- Q_G((1 - exp(-v))^(1 / alpha) / ((1 - exp(-v))^(1 / alpha) + exp(-v / alpha)))
  return(X)
}


#'
#' @name ZBOLLG
#' @examples
#' hzbollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hzbollg, -3, 3)
#' @export
hzbollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
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
  h <- alpha * g * G^(alpha - 1) * (1 - G)^(alpha - 1) * (-log(1 - G^alpha / (G^alpha + (1 - G)^alpha)))^(beta - 1) / (gamma(beta) * ((G^alpha) + (1 - G)^alpha)^2 * (1 - F0))
  return(h)
}

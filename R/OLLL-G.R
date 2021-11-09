#' Odd log-logistic logarithmic family of distributions (OLLL-G)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the beta extended distribution due to Haghbin et al. (2017) specified by the pdf
#' \deqn{f(x)=\frac{\alpha\beta\,g(x)\,G(x)^{\alpha-1}\bar{G}(x)^{\alpha-1}}{-\left[G(x)^\alpha+\bar{G}(x)^\alpha\right]\left[(1-\beta)\,G(x)^\alpha+\bar{G}(x)^\alpha\right]\log(1-\beta)}}
#' for \eqn{G} any valid cdf, \eqn{g} the corresponding pdf,  \eqn{\alpha > 0}, the first shape parameter, and \eqn{0 < \beta < 1}, the second shape parameter.
#'
#' @name OLLLG
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the value of the first shape parameter, must be positive, the default is 1.
#' @param beta the value of the second shape parameter, \eqn{0 < \beta < 1}, the default is 0.1.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{polllg} gives the distribution function,
#'  \code{dolllg} gives the density,
#'  \code{qolllg} gives the quantile function,
#'  \code{holllg} gives the hazard function and
#'  \code{rolllg} generates random variables from the Odd log-logistic logarithmic family of
#'  distributions (OLLL-G) for baseline cdf G.
#' @references Alizadeh, M., MirMostafee, S. M. T. K., Ortega, E. M., Ramires, T. G., Cordeiro, G. M. (2017). The odd log-logistic logarithmic generated family of distributions with applications in different areas. Journal of Statistical Distributions and Applications, 4(1), 1-25.
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' polllg(x)
#' polllg(x, alpha = 2, beta = .2, G = pbeta, shape1 = 1, shape2 = 2)
#'
#' @export
polllg <- function(x, alpha = 1, beta = .1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  F0 <- log(1 - (beta * G^alpha / (G^alpha + (1 - G)^alpha))) / log(1 - beta)
  return(F0)
}

#' @name OLLLG
#' @examples
#' dolllg(x, alpha = 2, beta = .2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(dolllg, -3, 3)
#' @importFrom stats numericDeriv  pnorm  runif uniroot
#' @export
dolllg <- function(x, alpha = 1, beta = .1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * beta * g * G^(alpha - 1) * (1 - G)^(alpha - 1)
  return(df)
}


#'
#' @name OLLLG
#' @examples
#' qolllg(x, alpha = 2, beta = .2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qolllg <- function(q, alpha = 1, beta = .1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - polllg(t, alpha, beta, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name OLLLG
#' @examples
#' n <- 10
#' rolllg(n, alpha = 2, beta = .2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
rolllg <- function(n, alpha = 1, beta = .1, G = pnorm, ...) {
  u <- runif(n)
  Q_G <- function(y) qolllg(y, alpha, beta, G, ...)
  Q_G <- Vectorize(Q_G)
  X <- Q_G((1 - (1 - beta)^u)^(1 / alpha) / ((1 - (1 - beta)^u)^(1 / alpha) + (beta - 1 + (1 - beta)^u))^(1 / alpha))
  return(X)
}


#'
#' @name OLLLG
#' @examples
#' holllg(x, alpha = 2, G = pbeta, beta = .2, shape1 = 1, shape2 = 2)
#' curve(holllg, -3, 3)
#' @export
holllg <- function(x, alpha = 1, beta = .1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- alpha * beta * g * G^(alpha - 1) * (1 - G)^(alpha - 1) / (-1 *(G^alpha + (1 - G)^alpha) * ((1 - beta) * G^alpha + (1 - G)^alpha) * log((1 - beta) / (1 - (beta * G^alpha / (G^alpha + (1 - G)^alpha)))))
  return(h)
}

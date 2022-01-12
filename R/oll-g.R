#' Odd log-logistic family of distributions (OLL-G)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the beta extended distribution due to Gleaton et al. (2006) specified by the pdf
#' \deqn{f=\frac{\alpha\,g\,G^{\alpha-1}\bar{G}^{\alpha-1}}{[G^\alpha+\bar{G}^\alpha]^2}}
#' for \eqn{G} any valid continuous cdf , \eqn{\bar{G}=1-G}, \eqn{g} the corresponding pdf,  \eqn{\alpha > 0}, the first shape parameter.
#'
#' @name OLLG
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the value of the first shape parameter, must be positive, the default is 1.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{pollg} gives the distribution function,
#'  \code{dollg} gives the density,
#'  \code{qollg} gives the quantile function,
#'  \code{hollg} gives the hazard function and
#'  \code{rollg} generates random variables from the Odd log-logistic family of
#'  distributions (OLL-G) for baseline cdf G.
#' @references Gleaton, J. U., Lynch, J. D. (2006). Properties of generalized log-logistic families of lifetime distributions. Journal of Probability and Statistical Science, 4(1), 51-64.
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' pollg(x)
#' pollg(x, alpha = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
pollg <- function(x, alpha = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  F0 <- G^alpha / (G^alpha + (1 - G)^alpha)
  return(F0)
}

#'
#' @name OLLG
#' @examples
#' dollg(x, alpha = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(dollg, -3, 3)
#' @importFrom stats numericDeriv  pnorm  runif uniroot
#' @export
dollg <- function(x, alpha = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * g * G^(alpha - 1) * (1 - G)^(alpha - 1) / (G^alpha + (1 - G)^alpha)^2
  return(df)
}


#'
#' @name OLLG
#' @examples
#' qollg(x, alpha = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qollg <- function(q, alpha = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pollg(t, alpha, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name OLLG
#' @examples
#' n <- 10
#' rollg(n, alpha = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
rollg <- function(n, alpha = 1, G = pnorm, ...) {
  u <- runif(n)
  Q_G <- function(y) qollg(y, alpha, G, ...)
  Q_G <- Vectorize(Q_G)
  X <- Q_G(u^(1 / alpha) / (u^(1 / alpha) + (1 - u)^(1 / alpha)))
  return(X)
}


#'
#' @name OLLG
#' @examples
#' hollg(x, alpha = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hollg, -3, 3)
#' @export
hollg <- function(x, alpha = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- alpha * g * G^(alpha - 1) / ((1 - G) * (G^alpha + (1 - G)^alpha))
  return(h)
}







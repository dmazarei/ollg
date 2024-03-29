#' Marshal-Olkin Odd log-logistic family of distributions (MOOLL-G)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the beta extended distribution due to Gleaton et al. (2010) specified by the pdf
#' \deqn{f=\frac{\alpha\beta\,g\,G^{\alpha-1}\bar{G}^{\alpha-1}}{[G^\alpha+\beta\,\bar{G}^\alpha]^2}}
#' for \eqn{G} any valid continuous cdf , \eqn{\bar{G}=1-G}, \eqn{g} the corresponding pdf,  \eqn{\alpha > 0}, the first shape parameter, and \eqn{\beta > 0}, the second shape parameter.
#'
#' @name MOOLLG
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the value of the first shape parameter, must be positive, the default is 1.
#' @param beta the value of the second shape parameter, must be positive, the default is 1.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{pmoollg} gives the distribution function,
#'  \code{dmoollg} gives the density,
#'  \code{qmoollg} gives the quantile function,
#'  \code{hmoollg} gives the hazard function and
#'  \code{rmoollg} generates random variables from the Marshal-Olkin Odd log-logistic family of
#'  distributions (MOOLL-G) for baseline cdf G.
#'
#' @references Gleaton, J. U., Lynch, J. D. (2010). Extended generalized loglogistic families of lifetime distributions with an application. J. Probab. Stat.Sci, 8(1), 1-17.
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' pmoollg(x)
#' pmoollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
pmoollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  F0 <- G^alpha / (G^alpha + beta * (1 - G)^alpha)
  return(F0)
}

#'
#' @name MOOLLG
#' @examples
#' dmoollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(dmoollg, -3, 3)
#' @importFrom stats numericDeriv  pnorm  runif uniroot
#' @export
dmoollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * beta * g * G^(alpha - 1) * (1 - G)^(alpha - 1) / (G^alpha + beta * (1 - G)^alpha)^2
  return(df)
}


#'
#' @name MOOLLG
#' @examples
#' qmoollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qmoollg <- function(q, alpha = 1, beta = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pmoollg(t, alpha, beta, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name MOOLLG
#' @examples
#' n <- 10
#' rmoollg(n, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
rmoollg <- function(n, alpha = 1, beta = 1, G = pnorm, ...) {
  u <- runif(n)
  Q_G <- function(y) qmoollg(y, alpha, beta, G, ...)
  Q_G <- Vectorize(Q_G)
  X <- Q_G((beta * u)^(1 / alpha) / ((beta * u)^(1 / alpha) + (1 - u)^(1 / alpha)))
  return(X)
}


#'
#' @name MOOLLG
#' @examples
#' hmoollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hmoollg, -3, 3)
#' @export
hmoollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- alpha * g * G^(alpha - 1) / ((1 - G) * (G^alpha + beta * (1 - G)^alpha))
  return(h)
}

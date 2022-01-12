#' Exponentiated Odd log-logistic family of distributions (EOLL-G)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers
#' of the beta extended distribution due to Alizadeh et al. (2020) specified by the pdf
#' \deqn{f=\frac{\alpha\beta\,g\,G^{\alpha\beta-1}\bar{G}^{\alpha-1}}{[G^\alpha+\bar{G}^\alpha]^{\beta+1}}}
#' for \eqn{G} any valid continuous cdf , \eqn{\bar{G}=1-G}, \eqn{g} the corresponding pdf,  \eqn{\alpha > 0}, the first shape parameter, and \eqn{\beta > 0}, the second shape parameter.
#'
#' @name EOLLG
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the value of the first shape parameter, must be positive, the default is 1.
#' @param beta the value of the second shape parameter, must be positive, the default is 1.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{peollg} gives the distribution function,
#'  \code{deollg} gives the density,
#'  \code{qeollg} gives the quantile function,
#'  \code{heollg} gives the hazard function and
#'  \code{reollg} generates random variables from the Exponentiated Odd log-logistic family of
#'  distributions (EOLL-G) for baseline cdf G.
#' @references ALIZADEH, Morad; TAHMASEBI, Saeid; HAGHBIN, Hossein. The exponentiated odd log-logistic family of distributions: Properties and applications. Journal of Statistical Modelling: Theory and Applications, 2020, 1. Jg., Nr. 1, S. 29-52.
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' peollg(x)
#' peollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
peollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  F0 <- G^(alpha * beta) / (G^alpha + (1 - G)^alpha)^beta
  return(F0)
}

#'
#' @name EOLLG
#' @examples
#' deollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(deollg, -3, 3)
#' @importFrom stats numericDeriv  pnorm  runif uniroot
#' @export
deollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * beta * g * G^(alpha * beta - 1) * (1 - G)^(alpha - 1) / ((G^alpha) + (1 - G)^alpha)^(beta + 1)
  return(df)
}


#'
#' @name EOLLG
#' @examples
#' qeollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qeollg <- function(q, alpha = 1, beta = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - peollg(t, alpha, beta, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name EOLLG
#' @examples
#' n <- 10
#' reollg(n, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
reollg <- function(n, alpha = 1, beta = 1, G = pnorm, ...) {
  u <- runif(n)
  Q_G <- function(y) qeollg(y, alpha, beta, G, ...)
  X <- Q_G(u^(1 / (alpha * beta)) / (u^(1 / (alpha * beta)) + (1 - u^((1 / beta)))^(1 / alpha)))
  return(X)
}


#'
#' @name EOLLG
#' @examples
#' heollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(heollg, -3, 3)
#' @export
heollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- alpha * beta * g * G^(alpha * beta - 1) * (1 - G)^(alpha - 1) / (((G^alpha) + (1 - G)^alpha) * (((G^alpha) + (1 - G)^alpha)^beta - G^(alpha * beta)))
  return(h)
}

#' Exponentiated Odd log-logistic family of distributions (EOLL-G)
#'
#'  Distribution function, density, quantile function, hazard
#'  function and random generation for Exponentiated Odd log-logistic family
#'  of distributions (EOLL-G) with baseline cdf G.
#'
#' @name EOLLG
#' @param x,q A numeric/quantiles	vector.
#' @param n number of observations. If \code{length(n) > 1},
#'  the length is taken to be the number required.
#' @param alpha non-negative parameters.
#' @param beta non-negative parameters.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{peollg} gives the distribution function,
#'  \code{deollg} gives the density,
#'  \code{qeollg} gives the quantile function,
#'  \code{heollg} gives the hazard function and
#'  \code{reollg} generates random variables from the Exponentiated Odd log-logistic family of
#'  distributions (EOLL-G) for baseline cdf G.
#'
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
#' heollg(x, alpha = 2, beta = 1, G = pbeta, shape1 = 1, shape2 = 2)
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

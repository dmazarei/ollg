#' The beta Odd log-logistic family of distributions (BOLL-G)
#'
#'  Distribution function, density, quantile function, hazard
#'  function and random generation for The beta Odd log-logistic family
#'  of distributions (BOLL-G) with baseline cdf G.
#'
#' @name BOLLG
#' @param x,q A numeric/quantiles	vector.
#' @param n number of observations. If \code{length(n) > 1},
#'  the length is taken to be the number required.
#' @param alpha non-negative parameter.
#' @param a non-negative parameter.
#' @param b non-negative parameter.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{pbollg} gives the distribution function,
#'  \code{dbollg} gives the density,
#'  \code{qbollg} gives the quantile function,
#'  \code{hbollg} gives the hazard function and
#'  \code{rbollg} generates random variables from the The beta Odd log-logistic family of
#'  distributions (BOLL-G) for baseline cdf G.
#' @references Cordeiro, G. M., Alizadeh, M., Tahir, M. H., Mansoor, M., Bourguignon, M., Hamedani, G. G. (2016). The beta odd log-logistic generalized family of distributions. Hacettepe Journal of Mathematics and Statistics, 45(4), 1175-1202.
#' @importFrom stats numericDeriv  pnorm  rbeta pbeta  uniroot  integrate
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' pbollg(x)
#' pbollg(x, alpha = 2, a = 1, b = 1, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
pbollg <- function(x, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  u <- G^alpha / (G^alpha + (1 - G)^alpha)
  F0 <- pbeta(u,a, b) - pbeta(0,a, b)
  return(F0)
}


#'
#' @name BOLLG
#' @examples
#' dbollg(x, alpha = 2, a = 1, b = 1, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(dbollg, -3, 3)
#' @export
dbollg <- function(x, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * g * G^(a * alpha - 1) * (1 - G)^(b * alpha - 1) / (beta(a, b) * ((G^alpha) + (1 - G)^alpha)^(a + b))
  return(df)
}


#'
#' @name BOLLG
#' @examples
#' qbollg(x, alpha = 2, a = 1, b = 1, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qbollg <- function(q, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pbollg(t, alpha, a, b, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name BOLLG
#' @examples
#' n <- 10
#' rbollg(n, alpha = 2, a = 1, b = 1, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
rbollg <- function(n, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  v <- rbeta(n, a, b)
  Q_G <- function(y) qbollg(y, alpha, a, b, G, ...)
  X <- Q_G(v^(1 / alpha) / (v^(1 / alpha) + (1 - v)^(1 / alpha)))
  return(X)
}


#'
#' @name BOLLG
#' @examples
#' hbollg(x, alpha = 2, a = 1, b = 1, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hbollg, -3, 3)
#' @export
hbollg <- function(x, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- alpha * g * G^(a * alpha - 1) * (1 - G)^(b * alpha - 1) / (beta(a, b) * ((G^alpha) + (1 - G)^alpha)^(a + b) * (1 - alpha * g * G^(a * alpha - 1) * (1 - G)^(b * alpha - 1) / (beta(a, b) * ((G^alpha) + (1 - G)^alpha)^(a + b))))
  return(h)
}

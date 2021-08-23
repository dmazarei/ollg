#'  New Odd log-logistic family of distributions (NOLL-G)
#'
#'  Distribution function, density, quantile function, hazard
#'  function and random generation for New Odd log-logistic family
#'  of distributions (NOLL-G) with baseline cdf G.
#'
#' @name NOLLG
#' @param x,q A numeric/quantiles	vector.
#' @param n number of observations. If \code{length(n) > 1},
#'  the length is taken to be the number required.
#' @param alpha non-negative parameters.
#' @param beta non-negative parameters.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{pnollg} gives the distribution function,
#'  \code{dnollg} gives the density,
#'  \code{qnollg} gives the quantile function,
#'  \code{hnollg} gives the hazard function and
#'  \code{rnollg} generates random variables from the New Odd log-logistic family of
#'  distributions (NOLL-G) for baseline cdf G.
#' @references Alizadeh, M., Altun, E., Ozel, G., Afshari, M., Eftekharian, A. (2019). A new odd log-logistic lindley distribution with properties and applications. Sankhya A, 81(2), 323-346.
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' pnollg(x)
#' pnollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
pnollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  F0 <- G^alpha / (G^alpha + (1 - G)^beta)
  return(F0)
}

#'
#' @name NOLLG
#' @examples
#' dnollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(dnollg, -3, 3)
#' @importFrom stats numericDeriv  pnorm  runif uniroot
#' @export
dnollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- g * G^(alpha - 1) * (1 - G)^(beta - 1) * (alpha + (beta - alpha) * G) / (G^alpha + (1 - G)^beta)^2
  return(df)
}


#'
#' @name NOLLG
#' @examples
#' qnollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qnollg <- function(q, alpha = 1, beta = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pnollg(t, alpha, beta, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name NOLLG
#' @examples
#' n <- 10
#' rnollg(n, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
rnollg <- function(n, alpha = 1, beta = 1, G = pnorm, ...) {
  u <- runif(n)
  Q_G <- function(y) qnollg(y, alpha, beta, G, ...)
  X <- Q_G(u)
  return(X)
}

#'
#' @name NOLLG
#' @examples
#' hnollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hnollg, -3, 3)
#' @export
hnollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- g * G^(alpha - 1) * (alpha + (beta - alpha) * G) / ((1 - G) * (G^alpha + (1 - G)^beta))
  return(h)
}

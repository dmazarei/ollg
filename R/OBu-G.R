#' Odd Burr generated family of distributions (OBu-G)
#'
#'  Distribution function, density, quantile function, hazard
#'  function and random generation for Odd Burr generated family
#'  of distributions (OBu-G) with baseline cdf G.
#'
#' @name OBuG
#' @param x,q A numeric/quantiles	vector.
#' @param n number of observations. If \code{length(n) > 1},
#'  the length is taken to be the number required.
#' @param alpha non-negative parameter.
#' @param beta non-negative parameter.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{pobug} gives the distribution function,
#'  \code{dobug} gives the density,
#'  \code{qobug} gives the quantile function,
#'  \code{hobug} gives the hazard function and
#'  \code{robug} generates random variables from the Odd Burr generated family of
#'  distributions (OBu-G) for baseline cdf G.
#' @references Alizadeh, M., Cordeiro, G. M., Nascimento, A. D., Lima, M. D. C. S., Ortega, E. M. (2017). Odd-Burr generalized family of distributions with some applications. Journal of statistical computation and simulation, 87(2), 367-389.
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' pobug(x)
#' pobug(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#'
#' @export
pobug <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  F0 <- 1 - (1 - (G^alpha / (G^alpha + (1 - G)^alpha)))^beta
  return(F0)
}

#'
#' @name OBuG
#' @examples
#' dobug(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(dobug, -3, 3)
#' @importFrom stats numericDeriv  pnorm  runif uniroot
#' @export
dobug <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * beta * g * G^(alpha - 1) * (1 - G)^(alpha * beta - 1) / (G^alpha + (1 - G)^alpha)^(beta + 1)
  return(df)
}


#'
#' @name OBuG
#' @examples
#' qobug(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qobug <- function(q, alpha = 1, beta = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pobug(t, alpha, beta, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name OBuG
#' @examples
#' n <- 10
#' robug(n, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
robug <- function(n, alpha = 1, beta = 1, G = pnorm, ...) {
  u <- runif(n)
  Q_G <- function(y) qobug(y, alpha, beta, G, ...)
  Q_G <- Vectorize(Q_G)
  X <- Q_G((1 - ( 1 - u)^(1 / beta) )^(1 / alpha) / ((1 - ( 1 - u)^(1 / beta) )^(1 / alpha) + (1 - u)^(1 / (alpha * beta))))
  return(X)
}


#'
#' @name OBuG
#' @examples
#' hobug(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hobug, -3, 3)
#' @export
hobug <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- alpha * beta * g * G^(alpha - 1) / ((1 - G) * (G^alpha + (1 - G)^alpha))
  return(h)
}

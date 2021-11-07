#' A New Odd log-logistic family of distributions (ANOLL-G)
#'
#'
#'
#' @name ANOLLG
#' @param x,q A numeric/quantiles	vector.
#' @param n number of random numbers to be generated.
#' @param alpha the value of the first shape parameter, must be positive, the default is 1
#' @param beta the value of the first shape parameter, must be positive, the default is 1
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{panollg} gives the distribution function,
#'  \code{danollg} gives the density,
#'  \code{qanollg} gives the quantile function,
#'  \code{hanollg} gives the hazard function and
#'  \code{ranollg} generates random variables from the A New Odd log-logistic family of
#'  distributions (ANOLL-G) for baseline cdf G.
#' @references Haghbin, Hossein, et al. "A new generalized odd log-logistic family of distributions." Communications in Statistics-Theory and Methods 46.20(2017): 9897-9920.
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' panollg(x)
#' panollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
panollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  F0 <- (1 - (1 - G)^alpha)^beta / ((1 - (1 - G)^alpha)^beta + (1 - G)^(alpha * beta))
  return(F0)
}

#'
#' @name ANOLLG
#' @examples
#' danollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(danollg, -3, 3)
#' @importFrom stats numericDeriv  pnorm  runif uniroot
#' @export
danollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * beta * g * (1 - G)^(alpha * beta - 1) * (1 - (1 - G)^(alpha))^(beta - 1) / ((1 - (1 - G)^(alpha))^(beta - 1) + (1 - G)^(alpha * beta))^2
  return(df)
}


#'
#' @name ANOLLG
#' @examples
#' qanollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qanollg <- function(q, alpha = 1, beta = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - panollg(t, alpha, beta, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name ANOLLG
#' @examples
#' n <- 10
#' ranollg(n, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
ranollg <- function(n, alpha = 1, beta = 1, G = pnorm, ...) {
  u <- runif(n)
  Q_G <- function(y) qanollg(y, alpha, beta, G, ...)
  X <- Q_G(1 - ((1 - u)^(1 / (alpha * beta)) / (u^(1 / beta) + (1 - u)^(1 / beta)))^(1 / alpha))
  return(X)
}


#'
#' @name ANOLLG
#' @examples
#' hanollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hanollg, -3, 3)
#' @export
hanollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- alpha * beta * g * (1 - G)^(alpha * beta - 1) / ((1 - G) * ((1 - (1 - G)^(alpha))^(beta) + (1 - G)^(alpha * beta)))
  return(h)
}

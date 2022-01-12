#' Kumaraswamy Odd log-logistic family of distributions (KwOLL-G)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers
#' of the beta extended distribution due to Alizadeh et al. (2017) specified by the pdf
#' \deqn{f=\frac{a\,b\,\alpha\,g\,G^{a\,\alpha-1}\bar{G}^{\alpha-1}}{[G^\alpha+\bar{G}^\alpha]^{a+1}}\times \{1-[\frac{G^\alpha}{G^\alpha+\bar{G}^\alpha}]^a\}^{b-1}}
#' for \eqn{G} any valid continuous cdf , \eqn{\bar{G}=1-G}, \eqn{g} the corresponding pdf, \eqn{a, b > 0}, the shape parameter, \eqn{\alpha > 0}, the first shape parameter.
#'
#' @name KwOLLG
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the value of the first shape parameter, must be positive, the default is 1.
#' @param a the value of the shape parameter, must be positive, the default is 1.
#' @param b the value of the shape parameter, must be positive, the default is 1.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{pkwollg} gives the distribution function,
#'  \code{dkwollg} gives the density,
#'  \code{qkwollg} gives the quantile function,
#'  \code{hkwollg} gives the hazard function and
#'  \code{rkwollg} generates random variables from the Kumaraswamy Odd log-logistic family of
#'  distributions (KwOLL-G) for baseline cdf G.
#' @references Alizadeh, M., Emadi, M., Doostparast, M., Cordeiro, G. M., Ortega, E. M., Pescim, R. R. (2015). A new family of distributions: the Kumaraswamy odd log-logistic, properties and applications. Hacettepe Journal of Mathematics and Statistics, 44(6), 1491-1512.
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' pkwollg(x)
#' pkwollg(x, alpha = 2, a = 2, b = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
pkwollg <- function(x, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  F0 <- 1 - (1 - (G^alpha / (G^alpha + (1 - G)^alpha))^a)^b
  return(F0)
}

#'
#' @name KwOLLG
#' @examples
#' dkwollg(x, alpha = 2, a = 2, b = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(dkwollg, -3, 3)
#' @importFrom stats numericDeriv  pnorm  runif uniroot
#' @export
dkwollg <- function(x, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- (a * b * alpha * g * G^(a * alpha - 1) * (1 - G)^(alpha - 1) / ((G^alpha) + (1 - G)^alpha)^(a + 1)) * (1 - (G^alpha / (G^alpha + (1 - G)^alpha))^a)^(b - 1)
  return(df)
}


#'
#' @name KwOLLG
#' @examples
#' qkwollg(x, alpha = 2, a = 2, b = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qkwollg <- function(q, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pkwollg(t, alpha, a, b, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name KwOLLG
#' @examples
#' n <- 10
#' rkwollg(n, alpha = 2, a = 2, b = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
rkwollg <- function(n, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  u <- runif(n)
  Q_G <- function(y) qkwollg(y, alpha, a, b, G, ...)
  X <- Q_G((1 - (1 - u)^(1 / (b)))^(1 / (a * alpha)) / (((1 - (1 - u)^(1 / (b)))^(1 / (a * alpha))) + (1 - (1 - (1 - u)^(1 / b))^(1 / a))^(1 / alpha)))
  return(X)
}


#'
#' @name KwOLLG
#' @examples
#' hkwollg(x, alpha = 2, a = 2, b = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hkwollg, -3, 3)
#' @export
hkwollg <- function(x, alpha = 1, a = 1, b = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- a * b * alpha * g * G^(a * alpha - 1) * (1 - G)^(alpha - 1) / (((G^alpha) + (1 - G)^alpha)^(a + 1) * (1 - (G^alpha / (G^alpha + (1 - G)^alpha))^a))
  return(h)
}

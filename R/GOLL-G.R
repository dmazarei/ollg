#' Generalized Odd log-logistic family of distributions (GOLL-G)
#'
#' Computes the pdf, cdf, hdf, quantile and random numbers of the beta extended distribution due to Cordeiro et al. (2017) specified by the pdf
#' \deqn{f=\frac{\alpha\beta\,g\,G^{\alpha\beta-1}[1-G^\alpha]^{\beta-1}}{[G^{\alpha\beta}+[1-G^\alpha]^\beta]^2}}
#' for \eqn{G} any valid continuous cdf , \eqn{\bar{G}=1-G}, \eqn{g} the corresponding pdf,  \eqn{\alpha > 0}, the first shape parameter, and \eqn{\beta > 0}, the second shape parameter.
#'
#' @name GOLLG
#' @param x scaler or vector of values at which the pdf or cdf needs to be computed.
#' @param q scaler or vector of probabilities at which the quantile needs to be computed.
#' @param n number of random numbers to be generated.
#' @param alpha the value of the first shape parameter, must be positive, the default is 1.
#' @param beta the value of the second shape parameter, must be positive, the default is 1.
#' @param G A baseline continuous cdf.
#' @param ... The baseline cdf parameters.
#' @return  \code{pgollg} gives the distribution function,
#'  \code{dgollg} gives the density,
#'  \code{qgollg} gives the quantile function,
#'  \code{hgollg} gives the hazard function and
#'  \code{rgollg} generates random variables from the Generalized Odd log-logistic family of
#'  distributions (GOLL-G) for baseline cdf G.
#' @references Cordeiro, G.M., Alizadeh, M., Ozel, G., Hosseini, B., Ortega, E.M.M., Altun, E. (2017). The generalized odd log-logistic family of distributions : properties, regression models and applications. Journal of Statistical Computation and Simulation ,87(5),908-932.
#' @examples
#' x <- seq(0, 1, length.out = 21)
#' pgollg(x)
#' pgollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
pgollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G <- sapply(x, G, ...)
  F0 <- G^(alpha * beta) / (G^(alpha * beta) + ((1 - G)^alpha)^beta)
  return(F0)
}

#'
#' @name GOLLG
#' @examples
#' dgollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(dgollg, -3, 3)
#' @importFrom stats numericDeriv  pnorm  runif uniroot
#' @export
dgollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  df <- alpha * beta * g * G^(alpha * beta - 1) * ((1 - G)^alpha)^beta / (G^(alpha * beta) + ((1 - G)^alpha)^beta)^2
  return(df)
}


#'
#' @name GOLLG
#' @examples
#' qgollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
qgollg <- function(q, alpha = 1, beta = 1, G = pnorm, ...) {
  q0 <- function(x0) {
    if (x0 < 0 || x0 > 1) stop(message = "[Warning] 0 < x < 1.")
    F0 <- function(t) x0 - pgollg(t, alpha, beta, G, ...)
    F0 <- Vectorize(F0)
    x0 <- uniroot(F0, interval = c(-1e+15, 1e+15))$root
    return(x0)
  }
  return(sapply(q, q0))
}


#'
#' @name GOLLG
#' @examples
#' n <- 10
#' rgollg(n, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' @export
rgollg <- function(n, alpha = 1, beta = 1, G = pnorm, ...) {
  u <- runif(n)
  Q_G <- function(y) qgollg(y, alpha, beta, G, ...)
  X <- Q_G((u^(1 / beta) / (u^(1 / beta) + ((1 - u)^(1 / beta))))^(1 / alpha))
  return(X)
}


#'
#' @name GOLLG
#' @examples
#' hgollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
#' curve(hgollg, -3, 3)
#' @export
hgollg <- function(x, alpha = 1, beta = 1, G = pnorm, ...) {
  G0 <- function(y) G(y, ...)
  myenv <- new.env()
  myenv$par <- list(...)
  myenv$x <- as.numeric(x)
  g0 <- numericDeriv(quote(G0(x)), "x", myenv)
  g <- diag(attr(g0, "gradient"))
  G <- sapply(x, G0)
  h <- alpha * beta * g * G^(alpha * beta - 1) / (((G^(alpha * beta) + ((1 - G)^alpha))^beta) * (1 - G)^alpha)
  return(h)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NOLL-G.R
\name{NOLLG}
\alias{NOLLG}
\alias{pnollg}
\alias{dnollg}
\alias{qnollg}
\alias{rnollg}
\alias{hnollg}
\title{New Odd log-logistic family of distributions (NOLL-G)}
\usage{
pnollg(x, alpha = 1, beta = 1, G = pnorm, ...)

dnollg(x, alpha = 1, beta = 1, G = pnorm, ...)

qnollg(q, alpha = 1, beta = 1, G = pnorm, ...)

rnollg(n, alpha = 1, beta = 1, G = pnorm, ...)

hnollg(x, alpha = 1, beta = 1, G = pnorm, ...)
}
\arguments{
\item{x}{scaler or vector of values at which the pdf or cdf needs to be computed.}

\item{alpha}{the value of the first shape parameter, must be positive, the default is 1.}

\item{beta}{the value of the second shape parameter, must be positive, the default is 1.}

\item{G}{A baseline continuous cdf.}

\item{...}{The baseline cdf parameters.}

\item{q}{scaler or vector of probabilities at which the quantile needs to be computed.}

\item{n}{number of random numbers to be generated.}
}
\value{
\code{pnollg} gives the distribution function,
 \code{dnollg} gives the density,
 \code{qnollg} gives the quantile function,
 \code{hnollg} gives the hazard function and
 \code{rnollg} generates random variables from the New Odd log-logistic family of
 distributions (NOLL-G) for baseline cdf G.
}
\description{
Computes the pdf, cdf, hdf, quantile and random numbers of the beta extended distribution due to Alizadeh et al. (2019) specified by the pdf
\deqn{f=\frac{g\,G^{\alpha-1}\bar{G}^{\beta-1}[\alpha+(\beta-\alpha)G]}{[G^\alpha+\bar{G}^\beta]^2}}
for \eqn{G} any valid continuous cdf , \eqn{\bar{G}=1-G}, \eqn{g} the corresponding pdf,  \eqn{\alpha > 0}, the first shape parameter, and \eqn{\beta > 0}, the second shape parameter.
}
\examples{
x <- seq(0, 1, length.out = 21)
pnollg(x)
pnollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
dnollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
curve(dnollg, -3, 3)
qnollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
n <- 10
rnollg(n, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
hnollg(x, alpha = 2, beta = 2, G = pbeta, shape1 = 1, shape2 = 2)
curve(hnollg, -3, 3)
}
\references{
Alizadeh, M., Altun, E., Ozel, G., Afshari, M., Eftekharian, A. (2019). A new odd log-logistic lindley distribution with properties and applications. Sankhya A, 81(2), 323-346.
}

\name{TiCIR}
\alias{TiCIR}
\docType{data}
\title{
Simulated Bivariate Time In-Homogeneous CIR process.
}
\description{
A simulated diffusion of the form 
\bold{\deqn{dX_t = X_t(1+sin(2*pi*t)-X_t^2)dt +dW_t^1,}}
\bold{\deqn{dY_t = Y_t(1+sin(2*pi*t)-Y_t^2)dt +dW_t^2,}}
with \bold{\eqn{X_0 = 1}} and \bold{\eqn{Y_0 = 1}}.}

\usage{data(TiCIR)}
\format{
  A data frame with 101 observations on the following 3 variables.
  \describe{
    \item{\code{X_t}}{A numeric vector of X-coordinates for the simulated data.}
    \item{\code{Y_t}}{A numeric vector of X-coordinates for the simulated data.}
    \item{\code{t}}{A numeric vector of time nodes at which the simulated data were observed.}
  }
}
\examples{
\donttest{
data(TiCIR)
}
}
\keyword{datasets}

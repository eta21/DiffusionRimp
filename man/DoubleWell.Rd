\name{DoubleWell}
\alias{DoubleWell}
\docType{data}
\title{
Simulted Double Well with Shifting Modality
}
\description{
A simulated diffusion of the form 
\bold{\deqn{dX_t = X_t(1+sin(2*pi*t)-X^2)dt +dW_t,}}
with \bold{\eqn{X_0 = 1}}.}
\usage{data(DoubleWell)}
\format{
  A data frame with 201 observations on the following 2 variables.
  \describe{
    \item{\code{X_t}}{Vector containing the simulated trajectory.}
    \item{\code{t}}{Time points at which the diffusion was observed.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(DoubleWell) 
plot(rev(DoubleWell),type='l') 
}
\keyword{datasets}

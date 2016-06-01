\name{DiffusionRimp-package}
\alias{DiffusionRimp-package}
\alias{DiffusionRimp}
\docType{package}
\title{
Data-imputation and density approximations for diffusion processes.
}
\description{
A package for performing data imputation on discretely observed diffusion processes as well as calculating numerical approximations to transition and first passage time densities. 
}
\details{
\tabular{ll}{
Package: \tab DiffusionRimp\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2015-12-01\cr
License: \tab GPL (>= 2)\cr
}
Functions included in the package:
\tabular{lcl}{
\code{\link{RS.impute}}    \tab:\tab Perform inference on a diffusion model using the random walk Metropolis-Hastings algorithm using the data-imputation algorithm. \cr
\code{\link{BiRS.impute}}  \tab:\tab Perform inference on a bivariate diffusion model using the random walk Metropolis-Hastings algorithm using the data-imputation algorithm. \cr
\code{\link{MOL.density}}  \tab:\tab Calculate the transitional density of a diffusion model using the method of lines.\cr
\code{\link{BiMOL.density}}\tab:\tab  Calculate the transitional density of a bivariate diffusion model using the method of lines.\cr
\code{\link{MOL.passage}}\tab:\tab Calculate the first passage time density of a time-homogeneous diffusion model with fixed barriers (i.e., a two-barrier first passage time problem). \cr
\code{\link{BiMOL.passage}}\tab:\tab Calculate the first passage time density of a time-homogeneous bivariate diffusion model with fixed barriers (i.e., a four-barrier problem in two dimensions). \cr
\code{\link{MOL.aic}}*\tab:\tab Calculate a pseudo-AIC value for a diffusion model using the method of lines.  \cr
\code{\link{BiMOL.aic}}*\tab:\tab Calculate pseudo-AIC value for a bivariate diffusion model using the method of lines.\cr
}
* Functions use C++.
}
\author{
Etienne A.D. Pienaar \email{etiennead@gmail.com}
}


\keyword{ package }
\keyword{ C++}

\examples{
# example(RS.impute)
# example(MOL.density)
# example(MOL.passage)
}

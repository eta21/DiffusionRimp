\name{MOL.aic}
\alias{MOL.aic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate a Pseudo AIC Value for a Diffusion Model via the Method of Lines
}
\description{
\code{MOL.aic()} approximates the likelihood function for a diffusion model under a given dataset and parameter vector.
}
\usage{
MOL.aic(X, time, delt, xlims, N, theta, plt = TRUE, wrt = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Vector of coordinates at which the diffusion process was observed, where N is the number of observations.}  
  \item{time}{A vector of time nodes at which the process was observed.}
  \item{xlims}{Limits for the spatial nodes of the lattice. These limits should be wide enough for the transition density to be negligibly small at the endpoints (see note \bold{[2]}).}
  \item{N}{The number of nodes in the spatial domain at which to evaluate the transitional density. Increase \code{N} for more accurate approximations (see note \bold{[3]} and warning \bold{[2]}).}
  \item{delt}{Step size for the time domain of the lattice (see note \bold{[4]}).}
  \item{theta}{Parameter vector at which the AIC should be evaluated. Typically the parameter vector is calculated using \code{RS.impute()}}
  \item{plt}{Draw a plot of the calculation as it takes place.}
  \item{wrt}{Write a .cpp file. Useful for inspection purposes.}
}
\value{
\item{AIC}{An approximate AIC value.}
\item{likelihood}{The approximate likelihood value used in the calculation.}
\item{p}{The dimesnion of the parameter vector.}
}

\author{
Etienne A.D. Pienaar \email{etiennead@gmail.com}
}
\section{Interface}{
\code{MOL.aic} uses a function-input interface whereby the drift and diffusion components of the stochastic differential equation (SDE)
\bold{\deqn{dX_t = mu(X_t, t)dt+sigma(X_t, t)dW_t,}}
are defined as R-functions in the current workspace. That is by defining the drift and diffusion components
 
 \code{mu=function(X,t){some expression i.t.o. X and t}}
 
 \code{sig=function(X,t){some expression i.t.o. X and t}}

further analysis may be conducted by calling the function \code{MOL.aic()}.
}

\note{
\bold{Note [1]}: Although the spatial limits of the lattice is defined by the user using \code{xlims}, if the initial value \code{Xs} does not fall on one of the lattice nodes, then the lattice is shifted accordingly. 

\bold{Note [2]}: The approximation assumes that the entire density of the process falls on a finite interval. Thus the algorithm may breakdown for certain problems. Depending on the parameters of the process, the limits may be very far apart or near. Some experimentation may be required. Otherwise, set \code{autofind = TRUE} to find useful limits. Note also that larger \code{N} may be required for wider limits.

\bold{Note [3]}: Increasing \code{N} will likely require smaller \code{delt}, thus increasing computation time. For some problems, slight increases in \code{N} will require much smaller \code{delt}.

\bold{Note [4]}: \code{delt} is used as the step size for a 10(8)-th order Runge-Kutta method for solving the resulting system of ODEs. Note again the inter-play between \code{N} and \code{delt} (see note \bold{[3]}).
}

\section{Warning }{
\bold{Warning [1]}: Note that if the lattice is shifted, degeneracies may occur for certain drift/volatility specifications if the shifted lattice moves into non-nonsensical values of the drift/volatility functions' domains. 

\bold{Warning [2]}: Although increasing the spatial resolution of the lattice by increasing \code{N} improves approximations, instabilities will occur if \code{delt} is not sufficiently small. This tends to manifest as jagged/spiked solutions that oscillate between positive and negative values.
}

\seealso{\code{\link{RS.impute}}, \code{\link{BiMOL.aic}}.}

\examples{
#===============================================================================
}

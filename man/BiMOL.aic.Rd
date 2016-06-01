\name{BiMOL.aic}
\alias{BiMOL.aic}

\title{
Calculate a Pseudo AIC Value for a Bivariate Diffusion Model via the Method of Lines
}
\description{
\code{MOL.aic()} approximates the likelihood function for a bivariate diffusion model under a given dataset and parameter vector.
}
\usage{
BiMOL.aic(X, time, delt, xlims, ylims, N, theta, plt = TRUE, wrt = FALSE)
}

\arguments{
  \item{X}{N x 2 matrix of coordinates at which the diffusion process was observed, where N is the number of observations.}  
  \item{time}{A vector of time nodes at which the process was observed.}
  \item{delt}{Step size for the time domain of the lattice (see note \bold{[4]}). }
  \item{xlims}{X-limits for the spatial nodes of the lattice. These limits should be wide enough for the transition density to be negligibly small at the endpoints (see note \bold{[2]}).}
  \item{ylims}{Y-limits for the spatial nodes of the lattice. These limits should be wide enough for the transition density to be negligibly small at the endpoints (see note \bold{[2]}).}
  \item{N}{The number of nodes in each spatial domain at which to evaluate the transitional density. Increase \code{N} for more accurate approximations (see note \bold{[3]} and warning \bold{[2]}).}
 \item{theta}{Parameter vector at which the AIC should be evaluated. Typically the parameter vector is calculated using \code{RS.impute()}}
 \item{plt}{Draw a plot of the calculation as it takes place.}
 \item{wrt}{Write a .cpp file. Useful for inspection purposes.}
}

\value{
\item{AIC}{An approximate AIC value.}
\item{likelihood}{The approximate likelihood value used in the calculation.}
\item{p}{The dimesnion of the parameter vector.}
}

\references{
Hamdi, S., Schiesser, W. E. and Griffiths, G. W. 2007 Method of lines. \emph{Scholarpedia}, \bold{2(7)}:2859.
Revision #124335.

Feagin, T. 2007 A tenth-order Runge-Kutta method with error estimate. \emph{In Proceedings of the IAENG
Conf. on Scientifc Computing}.
}

\author{
Etienne A.D. Pienaar \email{etiennead@gmail.com}
}

\section{Interface}{
\code{BiMOL.aic} uses a function-input interface whereby the drift and diffusion components of the stochastic differential equation (SDE)
\bold{\deqn{dX_t = mu_1(X_t, Y_t, t)dt+sigma_1(X_t, Y_t, t)dW_t^1}}
\bold{\deqn{dY_t = mu_2(X_t, Y_t, t)dt+sigma_2(X_t, Y_t, t)dW_t^2,}}
are defined as R-functions in the current workspace. That is by defining the drift and diffusion components
 
 \code{mu1=function(X,Y,t){some expression i.t.o. X, Y and t}}
 
 \code{mu2=function(X,Y,t){some expression i.t.o. X, Y and t}}

 \code{sig11=function(X,Y,t){some expression i.t.o. X, Y and t}}

 \code{sig22=function(X,Y,t){some expression i.t.o. X, Y and t}}

further analysis may be conducted by calling the function \code{BiMOL.aic()}.
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

\seealso{\code{\link{MOL.aic}}, \code{\link{BiRS.impute}}.}

\examples{
#===============================================================================
}

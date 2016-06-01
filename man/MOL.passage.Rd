\name{MOL.passage}
\alias{MOL.passage}
\title{Approximate the First Passage Time Density of a Two-Barrier Problem for Time-Homogeneous Scalar Diffusions.}

\description{
For scalar diffusions with drift \code{mu=function(X){}} and diffusion \code{sig=function(X){}}, moving in relation to lower and upper bounds \code{barriers[1]} and \code{barriers[2]} respectively, \code{MOL.passage()} approximates a solution to the partial differential equation (PDE) that governs the evolution of the survaval distribution of the first passage time density via the method of lines (MOL). 
}

\usage{
MOL.passage(Xs, t, barriers, N, delt, mu, sig, desc = 1)
}

\arguments{
  \item{Xs}{Starting/Initial value for the diffusion process (see note \bold{[1]}).}
  \item{t}{Value (>0) giving the terminal point for the transition horizon (the final time at which to evaluate the transition density).}
  \item{barriers}{Limits for the spatial nodes of the lattice. These limits now represent the barriers in the spatial domain (see note \bold{[2]}).}
  \item{N}{The number of nodes in the spatial domain at which to evaluate the transitional density. Increase \code{N} for more accurate approximations (see note \bold{[3]} and warning \bold{[2]}).}
  \item{delt}{Step size for the time domain of the lattice (see note \bold{[4]}).}
   \item{mu}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
   \item{sig}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
  \item{desc}{The type of discretization used (see note \bold{[5]}). }
}

\value{
\item{surface}{A matrix giving the approximate survival probability over time for all starting values contained in the discretization of the interval enclosed by \code{barriers}.}
\item{density}{A vector containing the approximate first passage time density for trajectories starting at \code{Xs} (see note \bold{[i]}).}
\item{time}{A vector of time nodes at which the approximation was evaluated.}
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
\code{MOL.passage} uses a function-input interface whereby the drift and diffusion components of the time-homogeneous stochastic differential equation (SDE):
\bold{\deqn{dX_t = mu(X_t)dt+sigma(X_t)dW_t,}}
are defined as R-functions in the current workspace. That is by defining the drift and diffusion components
 
 \code{mu=function(X){some expression i.t.o. X}}
 
 \code{sig=function(X){some expression i.t.o. X}}

further analysis may be conducted by calling the function \code{MOL.passage()}.
}

\note{
\bold{Note [1]}: If the initial value \code{Xs} does not fall on one of the lattice nodes, then the first passage time density is calculated by linearly interpolating between approximations at the two nearest lattice nodes. 

\bold{Note [2]}:  Note that that enough nodes, \code{N}, are needed in order to generate a sufficiently accurate approximation, especially when \code{barriers[1]} and \code{barriers[2]} are far apart.

\bold{Note [3]}: Increasing \code{N} will likely require smaller \code{delt}, thus increasing computation time. For some problems, slight increases in \code{N} will require much smaller \code{delt}.

\bold{Note [4]}: \code{delt} is used as the step size for a 10(8)-th order Runge-Kutta method for solving the resulting system of ODEs. Note again the inter-play between \code{N} and \code{delt} (see note \bold{[3]}).

\bold{Note [5]}: When one of the barriers is sufficiently far away to not be accessible within the provided time-horizon, instabilities may occur and an alternative discretization may be required in order to ensure smooth operation. Possible values are \code{desc = 1} (close barriers), \code{desc = 2} (\code{barriers[1]} is inaccessible) and \code{desc = 3}  (\code{barriers[2]} is inaccessible).

}
\section{Warnings}{
\bold{Warning [1]:}
}
\seealso{\code{\link{MOL.density}}, \code{\link{BiMOL.density}}.}

\examples{
\donttest{
#===============================================================================
# Calculate an approximate FPT density for a cubic diffusion confined to
# a finite strip.
# Given dXt = 0.5Xt(1-Xt^2)dt+0.5dBt with X0 = 0.5, approximate the first passage
# time density to the boundaries [-0.5,2].
#===============================================================================

  mu  <- function(X){(0.5*X*(1-X^2))}
  sig <- function(X) {0.5}
  res.PDE <- MOL.passage(Xs=0.5,t=10,barriers=c(-0.5,2),N=51,delt=1/250)
  MOL.plot(res.PDE)
  
  plot(res.PDE$density~res.PDE$time, col ='#222299',type='l',
       main = 'First passage time density')

#-------------------------------------------------------------------------------
}
}

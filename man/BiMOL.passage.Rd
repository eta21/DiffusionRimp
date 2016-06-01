\name{BiMOL.passage}
\alias{BiMOL.passage}
\title{Approximate the First Passage Time Density of a Four-Barrier Problem for Time-Homogeneous Bivariate Diffusions.}
\description{
\code{BiMOL.passage()} approximates a solution to partial differential equation (PDE) that governs the evolution of the survaval distribution of the first passage time density of a bivariate diffusion passing through fixed thresholds \code{barriers[1]} or \code{barriers[2]} in the X-dimension and \code{barriers[3]} or \code{barriers[4]} in the Y-dimension.
}

\usage{
BiMOL.passage(Xs, Ys, t, barriers, N, delt, mu1, mu2, sig11, sig12, 
              sig21, sig22, desc = 1)
}

\arguments{
  \item{Xs}{Starting/Initial X-coordinate for the diffusion process (see note \bold{[1]}).}
  \item{Ys}{Starting/Initial Y-coordinate for the diffusion process (see note \bold{[1]}).}
  \item{t}{Value (>0) giving the terminal point for the transition horizon (the final time at which to evaluate the transition density).}
  \item{barriers}{Limits for the spatial nodes of the lattice. These limits now represent the barriers in the spatial domain (see note \bold{[2]}).}
  \item{N}{The number of nodes in the spatial domain at which to evaluate the transitional density. Increase \code{N} for more accurate approximations (see note \bold{[3]} and warning \bold{[2]}).}
  \item{delt}{Step size for the time domain of the lattice (see note \bold{[4]}).}
 \item{mu1}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
 \item{mu2}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
  \item{sig11}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
 \item{sig12}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
  \item{sig21}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
 \item{sig22}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
  \item{desc}{The type of discretization used (see note \bold{[5]}). }
}


\value{
\item{surface}{An array giving the approximate survival probability volume for all starting values contained in the discretization of the polygon enclosed by the barriers.}
\item{density}{A vector containing the approximate first passage time density for trajectories starting at (\code{Xs},\code{Ys}) (see note \bold{[i]}).}
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
\code{BiMOL.passage} uses a function-input interface whereby the drift and diffusion components of the time-homogeneous bivariate stochastic differential equation (SDE):
\bold{\deqn{dX_t = mu_1(X_t,Y_t)dt+sigma_1(X_t,Y_t)dW_t^1}}
\bold{\deqn{dY_t = mu_2(X_t,Y_t)dt+sigma_2(X_t,Y_t)dW_t^2,}}
are defined as R-functions in the current workspace. That is by defining the drift and diffusion components
 
 \code{mu1=function(X,Y){some expression i.t.o. X and Y}}
 
 \code{mu2=function(X,Y){some expression i.t.o. X and Y}}

 \code{sig11=function(X,Y){some expression i.t.o. X and Y}}

 \code{sig22=function(X,Y){some expression i.t.o. X and Y}}

further analysis may be conducted by calling the function \code{BiMOL.passage()}.
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

\seealso{\code{\link{MOL.passage}}, \code{\link{BiMOL.density}}.}

\examples{
 \donttest{
#===============================================================================
# Bi-cubic diffusion with concentration in the even quadrants.
#===============================================================================
 
 # Define drift and diffusion terms.
 mu1   <- function(X,Y){0.5*(1-X^2)*X-1*Y}
 mu2   <- function(X,Y){0.5*(1-Y^2)*Y-1*X}
 sig11 <- function(X,Y){1.0}
 sig22 <- function(X,Y){1.0}

 # Parameters of the problem.
 Xs    <- 0.5      # Starting X-coordinate
 Ys    <- 0.5      # Starting Y-coordinate
 t     <- 10       # Final horizon time
 Xbar  <- c(-2,2)  # Barriers in X-dim
 Ybar  <- c(-2,2)  # Barriers in Y-dim
 Nodes <- 51       # How many nodes x nodes (incl. ends)
 delt  <- 1/250    # Time step size

 res.PDE <- BiMOL.passage(Xs, Ys, t, c(Xbar, Ybar), Nodes, delt)
 MOL.plot(res.PDE)

#-------------------------------------------------------------------------------
}
}


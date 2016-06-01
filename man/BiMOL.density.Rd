\name{BiMOL.density}
\alias{BiMOL.density}
\title{Approximate the Transition Density of a Bivariate Diffusion with Arbitrary Drift and Volatility Specification}

\description{
\code{BiMOL.density} approximates the transition density of a bivariate diffusion on a lattice \code{[xlims[1],xlims[2]] x [ylims[1],ylims[2]] x [s,t]} with \code{N x N} spatial nodes and time discretization \code{delt}, via the method of lines. The method of lines approximates the solution of the Fokker-Planck equation by an \code{N x N}-dimensional system of ordinary differential equations (ODEs) evaluated on [s,t]. 
}

\usage{
BiMOL.density(Xs, Ys, s , t , xlims, ylims, N, delt, mu1, mu2, 
             sig11, sig12, sig21, sig22, final.only = FALSE, 
             show.lattice = FALSE)
}

\arguments{
  \item{Xs}{Starting/Initial X-coordinate for the diffusion process (see note \bold{[1]} and warning \bold{[1]}).}
  \item{Ys}{Starting/Initial Y-coordinate for the diffusion process (see note \bold{[1]} and warning \bold{[1]}).}  
  \item{s}{Starting time for the diffusion process.}
  \item{t}{Value (>s) giving the terminal point for the transition horizon (the final time at which to evaluate the transition density).}
  \item{xlims}{X-limits for the spatial nodes of the lattice. These limits should be wide enough for the transition density to be negligibly small at the endpoints (see note \bold{[2]}).}
  \item{ylims}{Y-limits for the spatial nodes of the lattice. These limits should be wide enough for the transition density to be negligibly small at the endpoints (see note \bold{[2]}).}
  \item{N}{The number of nodes in each spatial domain at which to evaluate the transitional density. Increase \code{N} for more accurate approximations (see note \bold{[3]} and warning \bold{[2]}).}
 \item{mu1}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
 \item{mu2}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
  \item{sig11}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
 \item{sig12}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
  \item{sig21}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
 \item{sig22}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
  \item{delt}{Step size for the time domain of the lattice (see note \bold{[4]}). }
  \item{final.only}{Should the transition density on the entire lattice be returned (\code{FALSE}) or only the at the terminal point of the transition horizon \code{t} (\code{TRUE}). Default = \code{FALSE}.}
  \item{show.lattice}{If \code{=TRUE} (default) then the X-Y lattice is drawn and the initial value is indicated in red.}
}



\value{
  \item{density}{3D array containing the density approximation (matrix if \code{final.only = TRUE}).}
  \item{Xt}{
  Vector of X-coordinates at which approximation was carried out.
  }
  \item{Yt}{Vector of Y-coordinates at which approximation was carried out.}
  \item{time}{Vector of time nodes at which the approximation was evaluated.}
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
\code{BiMOL.density} uses a function-input interface whereby the drift and diffusion components of the stochastic differential equation (SDE)
\bold{\deqn{dX_t = mu_1(X_t, Y_t, t)dt+sigma_1(X_t, Y_t, t)dW_t^1}}
\bold{\deqn{dY_t = mu_2(X_t, Y_t, t)dt+sigma_2(X_t, Y_t, t)dW_t^2,}}
are defined as R-functions in the current workspace. That is by defining the drift and diffusion components
 
 \code{mu1=function(X,Y,t){some expression i.t.o. X, Y and t}}
 
 \code{mu2=function(X,Y,t){some expression i.t.o. X, Y and t}}

 \code{sig11=function(X,Y,t){some expression i.t.o. X, Y and t}}

 \code{sig22=function(X,Y,t){some expression i.t.o. X, Y and t}}

further analysis may be conducted by calling the function \code{BiMOL.density()}.
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

\seealso{\code{\link{MOL.passage}}, \code{\link{MOL.density}}.}

\examples{
\donttest{
#===============================================================================
# Bi-cubic diffusion with oscillating spin.
#===============================================================================

 # Define drift and diffusion terms:
 mu1   <- function(X,Y,t){X*(1-X^2)+sin(2*pi*t/4)*Y}
 mu2   <- function(X,Y,t){Y*(1-Y^2)-sin(2*pi*t/4)*X}
 sig11 <- function(X,Y,t){0.5}
 sig22 <- function(X,Y,t){0.5}

 # Parameters of the problem:
 Xs     <- 1            # Starting X-coordinate
 Ys     <- 1            # Starting Y-coordinate
 s      <- 0            # Starting time
 t      <- 15           # Final horizon time
 Xlim   <- c(-2.2, 2.2) # Lattice endpoints in X dim
 Ylim   <- c(-2.2, 2.2) # Lattice endpoints in Y dim
 Nodes  <- 51           # How many nodes (incl. ends)
 delt   <- 1/100        # Time stepsize
 
 # Run the Method of Lines:
 res <- BiMOL.density(Xs, Ys, s, t, Xlim, Ylim, Nodes, delt)
 
 # Plot the density
 MOL.plot(res)
 
#-------------------------------------------------------------------------------
}
}

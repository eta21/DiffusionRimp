\name{MOL.density}
\alias{MOL.density}
\title{Approximate the Transition Density of a Scalar Diffusion with Arbitrary Drift and Volatility Specification}

\description{
For scalar diffusions with drift \code{mu=function(X,t){}} and diffusion \code{sig=function(X,t){}}, \code{MOL.density} approximates the transition density of a scalar diffusion on a lattice \code{[xlims[1],xlims[2]] x [s,t]} with \code{N} spatial nodes and time discretization \code{delt}, via the method of lines. The method of lines approximates the solution of the Fokker-Planck equation by an \code{N}-dimensional system of ordinary differential equations (ODEs) evaluated on [s,t].
}

\usage{
MOL.density(Xs, s, t, xlims, N = 31, delt, mu, sig, final.only = FALSE)
}

\arguments{
  \item{Xs}{Starting/Initial value for the diffusion process (see note \bold{[1]} and warning \bold{[1]}).}
  \item{s}{Starting time for the diffusion process.}
  \item{t}{Value (>s) giving the terminal point for the transition horizon (the final time at which to evaluate the transition density).}
  \item{xlims}{Limits for the spatial nodes of the lattice. These limits should be wide enough for the transition density to be negligibly small at the endpoints (see note \bold{[2]}).}
  \item{N}{The number of nodes in the spatial domain at which to evaluate the transitional density. Increase \code{N} for more accurate approximations (see note \bold{[3]} and warning \bold{[2]}).}
  \item{delt}{Step size for the time domain of the lattice (see note \bold{[4]}).}
   \item{mu}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
   \item{sig}{Optional (if drift and diffusion coefficients are already defined) arguments giving the drift and diffusion coefficients as text expressions (See \bold{Interface} below).}
  \item{final.only}{Should the transition density on the entire lattice be returned (\code{FALSE}) or only the at the terminal point of the transition horizon \code{t} (\code{TRUE}). Default = \code{FALSE}.}
}



\value{
  \item{density}{Matrix containing the density approximation (vector if \code{final.only = TRUE}).}
  \item{Xt}{Vector of spatial values at which approximation was carried out.}
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
\code{MOL.density} uses a function-input interface whereby the drift and diffusion components of the stochastic differential equation (SDE)
\bold{\deqn{dX_t = mu(X_t, t)dt+sigma(X_t, t)dW_t,}}
are defined as R-functions in the current workspace. That is by defining the drift and diffusion components
 
 \code{mu=function(X,t){some expression i.t.o. X and t}}
 
 \code{sig=function(X,t){some expression i.t.o. X and t}}

further analysis may be conducted by calling the function \code{MOL.density()}.
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

\seealso{\code{\link{MOL.passage}}, \code{\link{BiMOL.density}}.}

\examples{
\donttest{
#===============================================================================
# Cubic diffusion with time-varying coefficients.
#===============================================================================

# Define drift and diffusion terms.
mu  <- function(X,t){2*X*(1+0.5*sin(2*pi*t/4)-X^2)}
sig <- function(X,t){1}
res <- MOL.density(c(-1),0,15,c(-3,3),81,1/250)

MOL.plot(res)

#-------------------------------------------------------------------------------
}
}

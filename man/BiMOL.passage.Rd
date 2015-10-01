\name{BiMOL.passage}
\alias{BiMOL.passage}
\title{Approximate the First Passage Time Density of a Four-Barrier Problem for Time-Homogeneous Bivariate Diffusions.}
\description{
\code{BiMOL.passage()} approximates a solution to partial differential equation (PDE) that governs the evolution of the survaval distribution of the first passage time density of a bivariate diffusion passing through fixed thresholds \code{barriers[1]} or \code{barriers[2]} in the X-dimension and \code{barriers[3]} or \code{barriers[4]} in the Y-dimension.
}

\usage{
BiMOL.passage(Xs, Ys, t = 1, barriers, N = 31, delt = 1/100, desc = 1)
}

\arguments{
  \item{Xs}{Starting/Initial X-coordinate for the diffusion process (see note \bold{[1]}).}
  \item{Ys}{Starting/Initial Y-coordinate for the diffusion process (see note \bold{[1]}).}
  \item{t}{Value (>0) giving the terminal point for the transition horizon (the final time at which to evaluate the transition density).}
  \item{barriers}{Limits for the spatial nodes of the lattice. These limits now represent the barriers in the spatial domain (see note \bold{[2]}).}
  \item{N}{The number of nodes in the spatial domain at which to evaluate the transitional density. Increase \code{N} for more accurate approximations (see note \bold{[3]} and warning \bold{[2]}).}
  \item{delt}{Step size for the time domain of the lattice (see note \bold{[4]}).}
  \item{desc}{The type of discretization used (see note \bold{[5]}). }
}

\details{
DETAILS ABOUT THE METHODOLOGY
}
\value{
\item{volume}{An array giving the approximate survival probability volume for all starting values contained in the discretization of the polygon enclosed by the barriers.}
\item{density}{A vector containing the approximate first passage time density for trajectories starting at (\code{Xs},\code{Ys}) (see note \bold{[i]}).}
\item{time}{A vector of time nodes at which the approximation was evaluated.}
}
\references{
DANIELS!!!!!!

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

 \code{sig1=function(X,Y){some expression i.t.o. X and Y}}

 \code{sig2=function(X,Y){some expression i.t.o. X and Y}}

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
#===============================================================================
# Bivariate CIR process with interaction.
#===============================================================================
 rm(list=ls(all=TRUE))
  
 # Define drift and diffusion terms.
 mu1=function(X,Y){0.5*(6-X)+0.1*Y}
 mu2=function(X,Y){0.5*(6-Y)}
 sig1=function(X,Y){0.3*sqrt(X)}
 sig2=function(X,Y){0.3*sqrt(Y)}
  
 Xstart = 5.5213        # Starting X-coordinate
 Ystart = 5.5130        # Starting Y-coordinate
 tmax   = 20            # Final horizon time
 X_barriers =  c(4,8)   # Barriers in X-dim
 Y_barriers =  c(4,8)   # Barriers in Y-dim
 N.points = 51          # How many nodes x nodes (incl. ends)
 delt =1/250            # Time step size
  
 simulate=function(Xs,Ys,xLt,xUt,yLt,yUt,N=10000,delt=1/1000)
 {
    X=rep(Xs,N)
    Y=rep(Ys,N)
    t=0
    k=1
    times=X*0
    while(N>=1)
    {
      X=X+mu1(X,Y)*delt+sig1(X,Y)*rnorm(N,sd=sqrt(delt))
      Y=Y+mu2(X,Y)*delt+sig2(X,Y)*rnorm(N,sd=sqrt(delt))
      
      t=t+delt
      if(any((X<xLt)|(X>xUt)|(Y<xLt)|(Y>xUt)))
      {
        wh=which((X>xLt)&(X<xUt)&(Y>yLt)&(Y<yUt))
        X=X[wh]
        Y=Y[wh]
        Nnew=length(X)
        d=N-Nnew
        times[k:(k+d)]=t
        k=k+d+1
        N=Nnew
      }
    }
    return(list(times=times))
 }
  
 res.SIM=simulate(Xstart,Ystart,X_barriers[1],X_barriers[2],Y_barriers[1],Y_barriers[2],10000,1/1000)
  
  
 res=BiMOL.passage(Xstart,Ystart,tmax,c(X_barriers, Y_barriers),N.points ,delt)

 # Draw a nice plot
 x11()
 par(ask=F)
 
 hist(res.SIM$times,breaks=50,col='darkblue',freq=F,ylim=range(res$density))
 lines(res$density~res$time,col='red')

 ss=round(seq(1,dim(res$volume)[3],length=100))
 cols=colorRampPalette(c('blue','green','yellow','red'))
 for(i in ss)
 {
     nlevs=51
     levs=seq(0,max(res$volume[,,i]),length=nlevs)
     filled.contour(res$X,res$Y,res$volume[,,i],levels=levs,
                col=cols(nlevs+1),ylab='Y_s',xlab='X_s',
                plot.axes={abline(h=0,v=0);axis(1);axis(2);},
                main=paste0('Contour Plot of Survival Probability (Time=',
                round((i-1)*diff(res$time)[1],3),')'))
 }

#===============================================================================
# Bi-cubic diffusion with concentration in the even quadrants.
#===============================================================================
 # Define drift and diffusion terms.
 mu1=function(X,Y){0.5*(1-X^2)*X-1*Y}
 mu2=function(X,Y){0.5*(1-Y^2)*Y-1*X}
 sig1=function(X,Y){1.0+0*X}
 sig2=function(X,Y){1.0+0*Y}

 # Parameters of the problem.
 Xstart = 0.5           # Starting X-coordinate
 Ystart = 0.5           # Starting Y-coordinate
 tmax   = 10            # Final horizon time
 X_barriers =  c(-2,2)  # Barriers in X-dim
 Y_barriers =  c(-2,2)  # Barriers in Y-dim
 N.points = 51          # How many nodes x nodes (incl. ends)
 delt =1/500            # Time step size

 res=BiMOL.passage(Xstart,Ystart,tmax,c(X_barriers, Y_barriers),N.points ,delt)

 # Make a nice plot again
 ss=round(seq(1,dim(res$volume)[3],length=100))
 cols=colorRampPalette(c('blue','green','yellow','red'))
 for(i in ss)
 {
     nlevs=51
     levs=seq(0,max(res$volume[,,i]),length=nlevs)
     filled.contour(res$X,res$Y,res$volume[,,i],levels=levs,
                col=cols(nlevs+1),ylab='Y_s',xlab='X_s',
                plot.axes={abline(h=0,v=0);axis(1);axis(2);},
                main=paste0('Contour Plot of Survival Probability (Time=',
                round((i-1)*diff(res$time)[1],3),')'))
 }
#-------------------------------------------------------------------------------
}


\name{BiRS.impute}
\alias{BiRS.impute}
\title{Brownian Bridge Data Imputation for Bivariate Diffusion Processes.}

\description{\code{BiRS.impute} performs inference on bivariate diffusion processes with quite arbitrary drift functionals by imputing missing sample paths with Brownian bridges. The procedure was developed by Roberts and Stramer (2001) and subsequently extended in later papers (). Currently, the diffusion is assumed to take on the form 
\bold{\deqn{dX_t = mux(X_t, Y_t, t, theta)dt+sigma[1]f(X_t)dW_t^1}}
\bold{\deqn{dY_t = muy(X_t, Y_t, t, theta)dt+sigma[2]g(Y_t)dW_t^2,}}
where \bold{\eqn{mux(X_t, Y_t, t, theta)}} and \bold{\eqn{muy(X_t, Y_t, t, theta)}} may be defined with near-complete freedom as R-functions in the current workspace. \bold{\eqn{f(.)}} and \bold{\eqn{g(.)}} may take on predefined forms (see details).
}

\usage{
BiRS.impute(X, T.seq, M = 10, theta1, theta2, prop.sds1, prop.sds2, transform, 
             N.updates, exclude = NULL, bridge.plot = FALSE, plot.chain = TRUE)
}

\arguments{
  \item{X}{N x 2 matrix of coordinates at which the diffusion process was observed, where N is the number of observations.}  
  \item{T.seq}{A vector of time nodes at which the process was observed.}
  \item{M}{The number of points to impute between successive observations. Note that as the transition horizon increases, more points may be required in order to get desirable acceptance rates for Brownian bridges. Note that some transition horizons may be removed from likelihood calculations - see \code{exclude}.}
  \item{theta1}{Starting parameters for the drift parameters. }
  \item{theta2}{Starting parameters for \eqn{sigma[1]} and \eqn{sigma[2]}.}
  \item{prop.sds1}{Proposal standard deviations for the drift parameter chain where \code{theta1[i]_new ~ N(\code{theta1[i]_old},\code{prop.sds1[i]}) }}
  \item{prop.sds2}{Proposal standard deviations for the volatility parameter chain where \code{theta2[i]_new ~ N(\code{theta2[i]_old},\code{prop.sds2[i]}) }}
  \item{transform}{2-Component vector of indicators, each assuming values 1, 2 or 3, for which of the predefined volatility structures to impose - see note \bold{[3]}.}
  \item{N.updates}{The number of updates to perform.}
  \item{exclude}{ Vector indicating which transitions to exclude from the analysis. Default = \code{NULL}.}
  \item{bridge.plot}{Display imputed paths for successive updates - see note \bold{[4]}. Default = \code{FALSE}.}
  \item{plot.chain}{Display the resulting MCMC chain - see notes \bold{[5]}, \bold{[6]} and \bold{[7]}.}
}

\details{
DETAILS ABOUT THE METHODOLOGY
}
\value{

  \item{per.matrix}{Matrix containing the MCMC chain updates.}
  \item{acceptence.rate}{Acceptance rates for the two parameter chains.}
  \item{bridge.rate}{Average Acceptance rate for Brownian bridge updates.}
  \item{run.time}{The total run time of the algorithm.}
}
\references{
Dellaportas, P. 2006 Bayesian model selection for partially observed diffusion models. \emph{Biometrika}, \bold{93(4)}:
809.

Kalogeropoulos, K., Dellaportas, P., and Roberts, G. O. 2011 Likelihood based inference for correlated
diffusions. \emph{Canadian Journal of Statistics}, \bold{39(1)}:52--72.

Roberts, G. and Stramer, O. 2001 On inference for partially observed nonlinear diffusion models using the
Metropolis-Hastings algorithm. \emph{Biometrika}, \bold{88}:603--621.
}
\author{Etienne A.D. Pienaar \email{etiennead@gmail.com}}
\note{
\bold{Note [1]:} The functional form of the drift components may be defined by the user as R-functions. Although the function body may take on arbitrary forms both the name of the drift functions and the input structure must assume the form 

\code{mux = function(variable1, variable2, time, theta){...}} and

\code{muy = function(variable1, variable2, time, theta){...}}.

Note that \code{theta} is a reserved variable for the parameters. It is left to the user to ensure that the functional forms are not nonsensical. See the examples below.

\bold{Note [2]:} Both drift functions have to be defined in order for the algorithm to execute.

\bold{Note [3]:} The functional form for the volatility of each dimension may take on one of the following: \eqn{f(X_t)=1}, corresponding to constant volatility
, \eqn{f(X_t)=\sqrt{X_t}} corresponding to CIR type models and \eqn{f(X_t)=X_t} such as for geometric Brownian motion. Corresponding indicators are simply given buy 1, 2, and 3.

\bold{Note [4]:} Paths are imputed on a unit -volatility process. For display purposes the back-transformed imputed trajectories may plotted along with vertical lines indicating which bridges remained unchanged over successive updates. 

\bold{Note [5]:} By default the MCMC chain is displayed in a panel plot. Standard MCMC diagnostics may be performed on this chain which is returned in the value list as \code{per.matrix}.

\bold{Note [6]:} In addition to the MCMC chain, acceptance rates are given for both the drift vector parameter updates as well as the volatility parameter updates. A, target region for the acceptance rate trajectories is displayed in blue.

\bold{Note [7]:} Average acceptance rates are indicated for bridges per the transition number by blue bars. A target region for the bridge acceptance rates is given in light blue. Average acceptance rates that are lower than 60\% are indicated along with their respective transition number. Vertical grey lines indicate \code{exclude}-d transitions.

}

\seealso{\code{\link{RS.impute}}, \code{\link{BiMOL.density}}.}

\examples{ 
%\dontrun{
rm(list=ls(all=TRUE))
#===============================================================================
# Lets perform inference on a time in-homogeneous 2D CIR process.
#===============================================================================

 # Read in the simulated data.
 data(TiCIR)
 x=TiCIR
 plot(x$X_t~x$t,col='blue',type='l')
 lines(x$Y_t~x$t,col=2)

 # Define the drift functions.
 mux=function(X,Y,t,theta){theta[1]*(theta[2]-X)+theta[3]*Y}
 muy=function(X,Y,t,theta){theta[4]*(theta[5]+theta[6]*sin(2*pi*t/12)-Y)}

 # Impute 25 points per transition and ignore the 50th and 60th transition.
 res=BiRS.impute(cbind(x$X_t,x$Y_t),x$t,25,c(5,15,2,5,15,5),c(2,2),c(1,1,1,1,1,1)/10,c(1,1)/100,N.updates=50000,transform=c(2,2),exclude=c(50,60))
 comp=cbind(c(1.5,5,0.1,2.5,5,1,0.5,0.5),apply(res$per.matrix[10000:50000,],2,mean))
 colnames(comp) = c('True','Ests.')
 rownames(comp)=c(paste0('theta[',1:6,']'),paste0('sigma[',1:2,']'))
 print(data.frame(format(comp,digits=3)))
 
#-------------------------------------------------------------------------------
%}
}

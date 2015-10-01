\name{RS.impute}
\alias{RS.impute}

\title{Brownian Bridge Data Imputation for Scalar Diffusion Processes.}

\description{
\code{RS.impute} performs inference on bivariate diffusion processes with quite arbitrary drift functionals by imputing missing sample paths with Brownian bridges. The procedure was developed by Roberts and Stramer (2001) and subsequently extended in later papers (Dellaportas et al, 2006; Kalogeropoulos et al., 2011). Currently, the diffusion is assumed to take on the form 
\bold{\deqn{dX_t = mu(X_t, t, theta)dt+sigma[1]f(X_t)dW_t,}}
where \bold{\eqn{mu(X_t, theta)}} may be defined with near-complete freedom as R-functions in the current workspace. \bold{\eqn{f(.)}}  may take on predefined forms (see details).
}

\usage{
RS.impute(X, T.seq, M = 10, theta1, theta2, prop.sds1, prop.sds2, transform, 
           N.updates, exclude = NULL, bridge.plot = FALSE, plot.chain = TRUE) 
}

\arguments{
  \item{X}{Vector of coordinates at which the diffusion process was observed, where N is the number of observations.}  
  \item{T.seq}{A vector of time nodes at which the process was observed.}
  \item{M}{The number of points to impute between successive observations. Note that as the transition horizon increases, more points may be required in order to get desirable acceptance rates for Brownian bridges. Note that some transition horizons may be removed from likelihood calculations - see \code{exclude}.}
  \item{theta1}{Starting parameters for the drift parameters. }
  \item{theta2}{Starting parameters for \eqn{sigma}.}
  \item{prop.sds1}{Proposal standard deviations for the drift parameter chain where \code{theta1[i]_new ~ N(\code{theta1[i]_old},\code{prop.sds1[i]}) }}
  \item{prop.sds2}{Proposal standard deviations for the volatility parameter chain where \code{theta2_new ~ N(\code{theta2_old},\code{prop.sds2}) }}
  \item{transform}{An indicator assuming values 1, 2 or 3, for which of the predefined volatility structures to impose - see note \bold{[3]}.}
  \item{N.updates}{The number of updates to perform.}
  \item{exclude}{ Vector indicating which transitions to exclude from the analysis. Default = \code{NULL}.}
  \item{bridge.plot}{Display imputed paths for successive updates - see note \bold{[4]}. Default = \code{FALSE}.}
  \item{plot.chain}{Display the resulting MCMC chain - see notes \bold{[5]}, \bold{[6]} and \bold{[7]}.}
}

\details{
}

\value{
  \item{per.matrix}{Matrix containing the MCMC chain updates.}
  \item{acceptence.rate}{Vector of acceptance rates for the two parameter chains.}
  \item{bridge.rate}{Vector of average acceptance rates for Brownian bridge updates.}
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
\author{
Etienne A.D. Pienaar \email{etiennead@gmail.com}
}


\note{
\bold{Note [1]:} The functional form of the drift components may be defined by the user as R-functions. Although the function body may take on arbitrary forms both the name of the drift functions and the input structure must assume the form 

\code{mu = function(variable, time, theta){...}}.

Note that \code{theta} is a reserved variable for the parameters. It is left to the user to ensure that the functional forms do not degenerate. See the examples below.

\bold{Note [2]:} The drift function has to be defined in order for the algorithm to execute.

\bold{Note [3]:} The functional form for the volatility may take on one of the following: \eqn{f(X_t)=1}, corresponding to constant volatility
, \eqn{f(X_t)=\sqrt{X_t}} corresponding to CIR type models and \eqn{f(X_t)=X_t} such as for geometric Brownian motion. Corresponding indicators are simply given buy 1, 2, and 3.

\bold{Note [4]:} Paths are imputed on a unit -volatility process. For display purposes the back-transformed imputed trajectories may plotted along with vertical lines indicating which bridges remained unchanged over successive updates. 

\bold{Note [5]:} By default the MCMC chain is displayed in a panel plot. Standard MCMC diagnostics may be performed on this chain which is returned in the value list as \code{per.matrix}.

\bold{Note [6]:} In addition to the MCMC chain, acceptance rates are given for both the drift vector parameter updates as well as the volatility parameter updates. A, target region for the acceptance rate trajectories is displayed in blue.

\bold{Note [7]:} Average acceptance rates are indicated for bridges per the transition number by blue bars. A target region for the bridge acceptance rates is given in light blue. Average acceptance rates that are lower than 60\% are indicated along with their respective transition number. Vertical grey lines indicate \code{exclude}-d transitions.
}

\seealso{\code{\link{BiRS.impute}}, \code{\link{MOL.density}}.}

\examples{
%\dontrun{
rm(list=ls(all=TRUE))
#===============================================================================
# Lets perform inference on a cubic diffusion with time varying modality.
#===============================================================================

 # Load a simulated dataset.
 data(DoubleWell)
 x=DoubleWell
 plot(x$X_t~x$t,type='l') 

 # Define the drift function.
 mu=function(X,t,theta){theta[1]*X*(theta[2]+theta[3]*sin(2*pi*t/4)-X^2)}

 # Impute 50 points per transition.
 res=RS.impute(x$X_t,x$t,50,c(4,4,1),c(2),c(1,1,1)/10,c(0.01),N.updates=50000)
 comp=cbind(c(1,1,1,0.5),apply(res$per.matrix[10000:50000,],2,mean))
 colnames(comp) = c('True','Ests.')
 rownames(comp)=c(paste0('theta[',1:3,']'),paste0('sigma[',1,']'))
 print(data.frame(format(comp,digits=3)))
 
#-------------------------------------------------------------------------------
%}
}


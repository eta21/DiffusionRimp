\name{STA.density}
\alias{STA.density}
\title{
Approximate the stationary distribution of arbitrary scalar diffusion process
}
\description{
BLAH
}
\usage{
STA.density(xmin, xmax, dx)
}
\arguments{
  \item{xmin}{Lower spatial bound.}
  \item{xmax}{Upper spatial bound.}
  \item{dx}{Step size for solver.}
}
\details{
DETAILS ABOUT THE METHODOLOGY
}
\value{
  \item{density }{Vector containing the approximate stationary distribution.}
  \item{Xt }{Vector of points in the spatial domain at which the approxximation is calculated.}
}
\references{
REFERENCES
}
\author{
Etienne A.D. Pienaar \email{etiennead@gmail.com}
}
\note{\bold{Note [1]:}
}
\seealso{
}
\examples{
  #=============================================================================
  #
  #=============================================================================
  
    mu=function(x){0.5*(5-x)*(sin(2*pi*x)^2*(1+cos(2*pi*x/4)))}
    sig=function(x){0.5}
    res=STA.density(1,12,1/50)
    plot(res$density~res$Xt,type='l',col='blue',lwd=2)

    delt=1/100
    X=rep(5,25000)
   
    for(i in 1:10000)
    {
      X=(X+mu(X)*delt+sig(X)*rnorm(length(X),sd=sqrt(delt)))
    }

    hist(X,breaks=100,col='blue',freq=F,xlim=range(X))
    lines(res$density~res$Xt,type='l',col='red',lwd=2)
   
    res2=outer(res$density,res$density)
    library(rgl)
    persp3d(res2,col='green')
    
#-------------------------------------------------------------------------------
}


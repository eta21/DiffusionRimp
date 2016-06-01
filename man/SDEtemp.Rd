\name{SDEtemp}
\alias{SDEtemp}
\docType{data}
\title{
An unspecified temperature dataset.
}
\description{
Some data. 
}
\usage{data(SDEtemp)}
\format{
  A data frame with 216 observations on the following 4 variables.
  \describe{
    \item{\code{Date}}{The dates to which obsevations pertain (year-month-day).}
    \item{\code{Density}}{Density values.}
    \item{\code{Temp}}{Temperature in degress celsius.}
    \item{\code{Time.Nodes}}{A numeric encoding of the dates at a monthly resolution using Date[1] as time 0.}
  }
}

\examples{
\donttest{
data(SDEtemp)
attach(SDEtemp)
par(mfrow=c(2,1))
plot(log10(Density)~as.Date(Date),type='l',col='#222299',lwd=2)
plot(Temp~as.Date(Date),col='#BBCCEE',type='l',lwd=2)
}}
\keyword{datasets}

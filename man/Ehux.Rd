\name{Ehux}
\alias{Ehux}
\docType{data}
\title{
Abundance Data for Emiliania huxleyi and Corresponding Sea Surface Temperatures (2002 - 2009).
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(Ehux)}
\format{
  A data frame with 216 observations on the following 4 variables.
  \describe{
    \item{\code{Date}}{The dates to which obsevations pertain (year-month-day).}
    \item{\code{Cell.Density}}{# Cells per litre.}
    \item{\code{Surface.Temp}}{Surface temperature in degress celsius.}
    \item{\code{Time.Nodes}}{A numeric encoding of the dates at a monthly resolution using Date[1] as time 0.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Ehux)
attach(Ehux)
par(mfrow=c(2,1))
plot(log10(Cell.Density)~as.Date(Date),type='l',col='blue')
plot(log(Surface.Temp)~as.Date(Date),col='red',type='l')
}
\keyword{datasets}

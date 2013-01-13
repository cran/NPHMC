\name{e1684szdata}
\alias{e1684szdata}
\docType{data}
\title{
Eastern Cooperative Oncology Group (ECOG) data with treatment as only covariate
}
\description{
Sample Size Calculation Example Data
}
\usage{data(e1684szdata)}
\format{
  A data frame with 285 observations on the following 3 variables.
  \describe{
    \item{\code{Time}}{observed relapse-free time}
    \item{\code{Status}}{censoring indicator, 1=event of interest happens, and 0=censoring}
    \item{\code{X}}{arm indicator, 0=control and 1=treatment}
  }
}

\examples{
data(e1684szdata)
}
\keyword{datasets}

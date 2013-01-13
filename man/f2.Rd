\name{f2}
\alias{f2}
\title{
function two
}
\description{
The second integrate function
}
\usage{
f2(t, accrualtime, followuptime, accrualdist, survdist, k, lambda0)
}
\arguments{
  \item{t}{
time variable
}
  \item{accrualtime}{
accrual time. The default is 3
}
  \item{followuptime}{
follow up time. The default is 4
}
  \item{accrualdist}{
accrual rate distribution. It can be "uniform", "increasing" or "decreasing".
}
  \item{survdist}{
survival distribution of uncured patients, it can be "exp" or "weib"
}
  \item{k}{
the shape parameter of weibull distribution. By default, it is set to be 1 for exponential distribution
}
  \item{lambda0}{
the scale parameter of exponential distribution or weibull distribution
}
}
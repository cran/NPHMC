\name{f4}
\alias{f4}
\title{
function four
}
\description{
the fourth integrate function
}
\usage{
f4(t, accrualtime, followuptime, hazard, accrualdist, beta0, gamma0, pi0, survdist, k)
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
  \item{hazard}{
the scale parameter of exponential distribution or weibull distribution
}
  \item{accrualdist}{
accrual rate distribution. It can be uniform, increasing or decreasing.
}
  \item{beta0}{
log hazard ratio of uncured patients
}
  \item{gamma0}{
log odds ratio of cure rates between the two groups
}
  \item{pi0}{
cure rate for control group. The defaultis 0.1
}
  \item{survdist}{
survival distribution of uncured patients, it can be exponential or weibull
}
  \item{k}{
the shape parameter of weibull distribution. By default, it is set to be 1 for exponential distribution
}
}

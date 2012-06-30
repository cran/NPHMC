\name{f3}
\alias{f3}
\title{
function three
}
\description{
The third integrate function
}
\usage{
f3(t, hazard, beta0, gamma0, pi0, survdist, k)
}
\arguments{
  \item{t}{
time variable
}
  \item{hazard}{
the scale parameter of exponential distribution or weibull distribution
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

\name{f3}
\alias{f3}
\title{
function three
}
\description{
The third integrate function
}
\usage{
f3(t, beta0, gamma0, pi0, survdist, k, lambda0)
}
\arguments{
  \item{t}{
time variable
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
survival distribution of uncured patients, it can be "exp" or "weib"
}
  \item{k}{
the shape parameter of weibull distribution. By default, it is set to be 1 for exponential distribution
}
  \item{lambda0}{
the scale parameter of exponential distribution or weibull distribution
}
}

\name{f1}
\alias{f1}
\title{
function one
}
\description{
The first integrate function
}
\usage{
f1(t, hazard, survdist, k)
}
\arguments{
  \item{t}{
time variable
}
  \item{hazard}{
the scale parameter of exponential distribution or weibull distribution
}
  \item{survdist}{
survival distribution of uncured patients, it can be exponential or weibull
}
  \item{k}{
the shape parameter of weibull distribution. By default, it is set to be 1 for exponential distribution
}
}
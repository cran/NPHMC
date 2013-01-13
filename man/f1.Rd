\name{f1}
\alias{f1}
\title{
function one
}
\description{
The first integrate function
}
\usage{
f1(t, survdist, k, lambda0)
}
\arguments{
  \item{t}{
time variable
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
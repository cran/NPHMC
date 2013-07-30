\name{S0}
\alias{S0}
\title{
S0 Function
}
\description{
Baseline survival function for mixture cure model
}
\usage{
S0(t, pi0, survdist, k, lambda0)
}
\arguments{
  \item{t}{
time variable
}
  \item{pi0}{
cure rate for the control arm, which is between 0 and 1.
}
  \item{survdist}{
survival distribution of uncured patients. It can be "\code{exp}" or "\code{weib}".
}
  \item{k}{
if \code{survdist} = "\code{weib}", the shape parameter \code{k} needs to be specified. By default \code{k} = 1, which refers to the exponential distribution.
}
  \item{lambda0}{
scale parameter of exponential distribution or Weibull distribution for survival times of uncured patients in the control arm. 

 The density function of Weibull distribution with shape parameter \code{k} and scale parameter \eqn{\lambda_0} is given by
\deqn{f(t)=\lambda_{0}k(\lambda_{0}t)^{k-1}\exp(-(\lambda_{0}t)^k),} for \eqn{t > 0},  
 and the corresponding survival distribution is 
\deqn{S(t)=\exp(-(\lambda_0 t)^k).
 }
}
}

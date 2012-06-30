\name{S0}
\alias{S0}
\title{
S0 Function
}
\description{
Baseline Survival Function for Mixture Cure Model
}
\usage{
S0(t, hazard, pi0, survdist, k)
}
\arguments{
  \item{t}{
time variable
}
  \item{hazard}{
the scale parameter of exponential distribution or weibull distribution
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

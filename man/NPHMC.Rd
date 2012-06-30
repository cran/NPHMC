\name{NPHMC}
\alias{NPHMC}
\title{
An R-package for Estimating Sample Size of Proportional Hazard Cure Models
}
\description{
An R-package for Estimating Sample Size of Proportional Hazard Cure (PHMC) Model and Standard Proportional Hazard (PH) Model
}
\usage{
NPHMC(power = 0.9, alpha = 0.05, accrualtime = 3, followuptime = 4, p = 0.5, hazard = 1, accrualdist = c("uniform", "increasing", "decreasing"), hazardratio = 0.5, oddsratio = 2.15, pi0 = 0.1, survdist = c("exp", "weib"), k = 1, data = NULL)
}
\arguments{
  \item{power}{
required power. The default power is 0.9
}
  \item{alpha}{
required type I error. By default alpha is 0.05
}
  \item{accrualtime}{
accrual time. The default is 3
}
  \item{followuptime}{
follow up time. The default is 4
}
  \item{p}{
the proportion of subjects in each arm. The default is 0.5
}
  \item{hazard}{
the scale parameter of exponential distribution or weibull distribution
}
  \item{accrualdist}{
accrual rate distribution. It can be uniform, increasing or decreasing.
}
  \item{hazardratio}{
hazard ratio of uncured patients
}
  \item{oddsratio}{
odds ratio of cure rates between control and treatment groups
}
  \item{pi0}{
cure rate for control group. The default is 0.1
}
  \item{survdist}{
survival distribution of uncured patients, it can be exponential or weibull
}
  \item{k}{
the shape parameter of weibull distribution. By default, it is set to be 1 for exponential distribution
}
  \item{data}{
historical data with only one covarite. The data has four columns named "Time","Status","X" and "Z".
}
}
\details{
None
}

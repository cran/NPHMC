\name{NPHMC}
\alias{NPHMC}
\title{
An R-package for Estimating Sample Size of Proportional Hazard Cure Models
}
\description{
An R-package for Estimating Sample Size of Proportional Hazard Cure (PHMC) Model and Standard Proportional Hazard (PH) Model
}
\usage{
NPHMC(power = 0.8, alpha = 0.05, accrualtime = 3, followuptime = 4, p = 0.5, accrualdist = c("uniform", "increasing", "decreasing"), hazardratio = 0.5, oddsratio = 2.15, pi0 = 0.1, survdist = c("exp", "weib"), k = 1, lambda0 = 1, data = NULL)
}
\arguments{
  \item{power}{
specifies the required power. The default power is 0.8
}
  \item{alpha}{
specifies the level significance of the statistical test. The default is 0.05.
}
  \item{accrualtime}{
accrual time. The default is 3
}
  \item{followuptime}{
follow up time. The default is 4
}
  \item{p}{
specifies the proportion of subjects in each arm. The default p is 0.5
}
  \item{accrualdist}{
specifies the accrual rate distribution. It can be "uniform", "increasing" or "decreasing".
}
  \item{hazardratio}{
specifies the hazard ratio of uncured patients between two arms, which is equivalent to \eqn{e^{\beta_0}=\lambda_1(t) / \lambda_0(t)}. The value must be greater than 0.
}
  \item{oddsratio}{
specifies the odds ratio of cure rates between control and treatment arms, which is calculated by \eqn{e^{\gamma_0}=\frac{\pi_1}{1-\pi_1}/\frac{\pi_0}{1-\pi_0}}. The value should be greater than 0 if cure rates exist. It can be 0 if there's no cure rate.
}
  \item{pi0}{
specifies the cure rate for control arm which is between 0 and 1. The default is 0.1.
}
  \item{survdist}{
specifies the survival distribution of uncured patients, it can be "exp" or "weib".
}
  \item{k}{
if survdist = "weib", the shape parameter k needs to be specified. By default, k = 1, which refers to the exponential distribution.
}
  \item{lambda0}{
specifies the scale parameter of exponential distribution or Weibull distribution for survival times of uncured patients in the control arm. The default lambda0 = 1.\
 The density function of Weibull distribution with shape parameter k and scale parameter lambda0 is given by
\deqn{f(t)=\lambda_{0}k(\lambda_{0}t)^{k-1}\exp(-(\lambda_{0}t)^k)
 } for t > 0
}
  \item{data}{
if observed/historical data is avaialble, the sample size can be calculated based on the nonparametric estimators from the proportional hazards mixture model by 'smcure' package in R. The data must contain three columns with order of "Time","Status" and "X" where "Time" refers to time to event of interest, "Status" refers to censoring indicator, 1=event of interest happens, and 0=censoring and "X" refers to arm indicator, 0=control and 1=treatment.
}
}
\details{
None
}

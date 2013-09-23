\name{NPHMC}
\alias{NPHMC}
\title{
An R-package for Estimating Sample Size and Power of Proportional Hazards Mixture Cure Model
}
\description{
Estimating sample size and power of survival trial based on PH mixture cure model and standard PH model
}
\usage{
NPHMC(n=NULL, power = 0.8, alpha = 0.05, accrualtime = NULL, 
followuptime = NULL, p = 0.5, accrualdist = c("uniform", "increasing",
"decreasing"), hazardratio = NULL, oddsratio = NULL, pi0 = NULL, 
survdist = c("exp", "weib"), k = 1, lambda0 = NULL, data = NULL)
}
\arguments{
  \item{n}{
the sample size needed for the power calculation. 
}
  \item{power}{
the power needed for sample size calculation. The default \code{power} is \eqn{80\%}.
}
  \item{alpha}{
the level of significance of the statistical test. The default \code{alpha} is 0.05.
}
  \item{accrualtime}{
the length of accrual period.
}
  \item{followuptime}{
the length of follow-up time.
}
  \item{p}{
the proportion of subjects in the treatment arm. The default \code{p} is 0.5.
}
  \item{accrualdist}{
the accrual pattern. It can be "\code{uniform}", "\code{increasing}" or "\code{decreasing}".
}
  \item{hazardratio}{
the hazard ratio of uncured patients between two arms, which is defined as \eqn{e^{\beta_0}=\lambda_1(t) / \lambda_0(t)}. The value must be greater than 0 but cannot be 1 because \eqn{\beta_0} is the denominator of the sample size formula.
}
  \item{oddsratio}{
the odds ratio of cure rates between two arms, which is equivalent to \eqn{e^{\gamma_0}=\frac{\pi_1}{1-\pi_1}/\frac{\pi_0}{1-\pi_0}}. The value should be greater than 0 if there is cured fraction. When it is 0, the model is reduced to the standard proportional hazards model, which means there is no cure rate.
}
  \item{pi0}{
the cure rate for the control arm, which is between 0 and 1.
}
  \item{survdist}{
the survival distribution of uncured patients. It can be "\code{exp}" or "\code{weib}".
}
  \item{k}{
if \code{survdist} = "\code{weib}", the shape parameter \code{k} needs to be specified. By default \code{k} = 1, which refers to the exponential distribution.
}
  \item{lambda0}{
the scale parameter of exponential distribution or Weibull distribution for survival times of uncured patients in the control arm. 

 The density function of Weibull distribution with shape parameter \code{k} and scale parameter \eqn{\lambda_0} is given by
\deqn{f(t)=\lambda_{0}k(\lambda_{0}t)^{k-1}\exp(-(\lambda_{0}t)^k),} for \eqn{t > 0},  
 and the corresponding survival distribution is 
\deqn{S(t)=\exp(-(\lambda_0 t)^k).
 }
}
  \item{data}{
if observed/historical \code{data} is avaialble, the sample size can be calculated based on the nonparametric estimators from the proportional hazards mixture model by '\code{smcure}' package in R. The data must contain three columns with order of "\code{Time}","\code{Status}" and "\code{X}" where "\code{Time}" refers to time to event of interest, "\code{Status}" refers to censoring indicator (1=event of interest happens, and 0=censoring) and "\code{X}" refers to arm indicator (1=treatment and 0=control). By default, \code{data}=NULL.
}
}
\value{
when \code{data} is not specified, the package can return the following estimated sample size (power) values:
    \item{nsize }{estimated sample size based on the PH mixture cure model}
    \item{nsizeph }{estimated sample size based on the standard PH model}
    \item{pw }{estimated power based on the PH mixture cure model}
    \item{pwph }{estimated power based on the standard PH model}

when \code{data} is specified, the package will display the estimators from the "\code{smcure}" package in R and the returned values list as follows:
    \item{f }{a list of estimators from the smcure package}
    \item{nonpar }{estimated nonparametric sample size estimation based on the PH mixture cure model and observed data}
    \item{nonparPH }{estimated nonparametric sample size estimation based on the standard PH model and observed data}
    \item{nonparpw }{estimated nonparametric power estimation based on the PH mixture cure model and observed data}
    \item{nonparpwph }{estimated nonparametric power estimation based on the standard PH model and observed data}    
}
\examples{

# parametric sample size calculation
NPHMC(power=0.90,alpha=0.05,accrualtime=3,followuptime=4,p=0.5,accrualdist="uniform",
hazardratio=2/2.5,oddsratio=2.25,pi0=0.1,survdist="exp",k=1,lambda0=0.5)

# nonparametric sample size calculation
data(e1684szdata)
NPHMC(power=0.80,alpha=0.05,accrualtime=4,followuptime=3,p=0.5,accrualdist="uniform",
data=e1684szdata)

# parametric power calculation
n=seq(100, 500, by=50)
NPHMC(n=n, alpha=0.05,accrualtime=3,followuptime=4,p=0.5,
accrualdist="uniform", hazardratio=2/2.5,oddsratio=2.25,pi0=0.1,survdist="exp",
k=1,lambda0=0.5)

# nonparametric power calculation
n=seq(100, 500, by=50)
NPHMC(n=n,alpha=0.05,accrualtime=4,followuptime=3,p=0.5, 
accrualdist="uniform",data=e1684szdata)
}
\keyword{PH mixture cure model }
\keyword{sample size }
\keyword{power }



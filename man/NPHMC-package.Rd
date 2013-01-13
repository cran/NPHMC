\name{NPHMC-package}
\alias{NPHMC-package}
\docType{package}
\title{
An R-package for Estimating Sample Size of Proportional Hazards Mixture Cure Model
}
\description{
An R-package for Estimating Sample Size of Proportional Hazards Mixture Cure (PHMC) Model and Standard Proportional Hazards (PH) Model
}
\details{
\tabular{ll}{
Package: \tab NPHMC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-13\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Chao Cai, Songfeng Wang, Wenbin Lu, Jiajia Zhang

Maintainer: Chao Cai <caic@email.sc.edu>
}
\references{
S. Wang, J. Zhang, and W. Lu. Sample size calculation for the proportional hazards cure model. Statistics in medicine, 31:3959-3971, 2012

C. Cai, et al., smcure: An R-Package for estimating semiparametric mixture cure models, Computer Methods and Programs in Biomedicine. (2012), http://dx.doi.org/10.1016/j.cmpb.2012.08.013.

}
\examples{
# parametric
NPHMC(power=0.90,alpha=0.05,accrualtime=3,followuptime=4,p=0.5,accrualdist="uniform",
hazardratio=2/2.5,oddsratio=2.25,pi0=0.1,survdist="exp",k=1,lambda0=0.5)

# nonparametric
data(e1684szdata)
NPHMC(power=0.80,alpha=0.05,accrualtime=4,followuptime=3,p=0.5,accrualdist="uniform",
data=e1684szdata)
}

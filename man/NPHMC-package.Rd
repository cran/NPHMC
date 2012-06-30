\name{NPHMC-package}
\alias{NPHMC-package}
\docType{package}
\title{
An R-package for Estimating Sample Size of Proportional Hazard Models
}
\description{
An R-package for Estimating Sample Size of Proportional Hazard Cure (PHMC) Model and Standard Proportional Hazard (PH) Model
}
\details{
\tabular{ll}{
Package: \tab NPHMC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-07-01\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
}
\author{
Chao Cai, Songfeng Wang, Wenbin Lu, Jiajia Zhang

Maintainer: Chao Cai <caic@email.sc.edu>
}
\references{
None
}
\examples{
# parametric
NPHMC(power=0.9,alpha=0.05,accrualtime=3,followuptime=4,hazard=0.5,accrualdist="uniform",hazardratio=2/2.5,oddsratio=2.25,p=0.5,pi=0.1,survdist="exp",k=1)

# nonparametric
data(e1684szdata)
NPHMC(power=0.80,alpha=0.05,accrualtime=4,followuptime=3,p=0.5,accrualdist="uniform",data=e1684szdata)
}

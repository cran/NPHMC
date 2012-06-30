f4 <-
function(t,accrualtime,followuptime,hazard,accrualdist,beta0,gamma0,pi0,survdist,k){
 m(t,hazard,beta0,gamma0,pi0,survdist,k)*f2(t,accrualtime,followuptime,hazard,accrualdist,survdist,k)
}


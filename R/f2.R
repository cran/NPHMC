f2 <-
function(t,accrualtime,followuptime,hazard,accrualdist,survdist,k){
 Sc(t,accrualtime,followuptime,accrualdist)*f1(t,hazard,survdist,k)
}


m <-
function(t,hazard,beta0,gamma0,pi0,survdist,k){
 (gamma0/beta0+H0(t,survdist,hazard,k))*pi0/S0(t,hazard,pi0,survdist,k)-1}


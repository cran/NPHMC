S0 <-
function(t,hazard,pi0,survdist,k){
 if(survdist=="exp") {k=1;return(pi0+(1-pi0)*exp(-hazard*t^k))}
 if(survdist=="weib") {return(pi0+(1-pi0)*exp(-hazard*t^k))}
}


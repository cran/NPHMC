f1 <-
function(t,hazard,survdist,k){
 if(survdist=="exp") {k=1; return(hazard*k*(hazard*t)^(k-1)*exp(-(hazard*t)^k))}
 if(survdist=="weib") {return(hazard*k*(hazard*t)^(k-1)*exp(-(hazard*t)^k))}
}


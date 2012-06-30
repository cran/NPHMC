H0 <-
function(t,survdist,hazard,k){
 if(survdist=="exp") {return(hazard*t)}
 if(survdist=="weib") {return((hazard*t)^k)}
}


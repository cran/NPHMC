NPHMC <-
function(power=0.8,alpha=0.05,accrualtime=3,followuptime=4,p=0.5,accrualdist=c("uniform","increasing","decreasing"),
		hazardratio=0.5,oddsratio=2.15,pi0=0.1,survdist=c("exp","weib"),k=1,lambda0=1,data=NULL){
	n<-list()
	class(n) <- c("NPHMC")  
if (is.null(data)){
      if (hazardratio<=0) stop("Hazardratio must be greater than 0")
	if (oddsratio<0) stop("Oddsratio cannot be less than 0")
        if (pi0==0 | oddsratio==0) {
            i1 <- integrate(f1,0,followuptime,survdist,k,lambda0)$value
	      i2 <- integrate(f2,followuptime,(accrualtime+followuptime),accrualtime,followuptime,accrualdist,survdist,k,lambda0)$value	  
	      beta0 <- log(hazardratio)
         	pdeath <- i1+i2
	      nsizeph <- ceiling((qnorm(power)-qnorm(alpha/2))^2/(p*(1-p)*beta0^2*pdeath))
            cat("====================================================================== \n")
	      cat("SAMPLE SIZE CALCULATION BASED ON STANDARD PH MODEL (NO CURE FRACTION) \n")
	      cat("====================================================================== \n")
            cat("PH Standard Model: n =",nsizeph,"\n")}
        else {
	  i1 <- integrate(f1,0,followuptime,survdist,k,lambda0)$value
	  i2 <- integrate(f2,followuptime,(accrualtime+followuptime),accrualtime,followuptime,accrualdist,survdist,k,lambda0)$value	  
	  beta0 <- log(hazardratio)
	  gamma0 <- log(oddsratio)
	  i3 <- integrate(f3,0,followuptime,beta0,gamma0,pi0,survdist,k,lambda0)$value
	  i4 <- integrate(f4,followuptime,(accrualtime+followuptime),accrualtime,followuptime,accrualdist,beta0,gamma0,pi0,survdist,k,lambda0)$value
	nsize <- ceiling((qnorm(power)-qnorm(alpha/2))^2*(i1+i2)/((i3+i4)^2*p*(1-p)*(1-pi0)*beta0^2))
	pdeath <- i1+i2
	nsizeph <- ceiling((qnorm(power)-qnorm(alpha/2))^2/(p*(1-p)*beta0^2*pdeath))
	cat("\n")
	cat("======================================================================== \n")
	cat("SAMPLE SIZE CALCULATION FOR PH MIXTURE CURE MODEL AND STANDARD PH MODEL \n")
	cat("======================================================================== \n")
	cat("PH Mixture Cure Model: n =",nsize,"\n")
	cat("PH Standard Model: n =",nsizeph,"\n")
      n$nsize <- nsize }	
     }
  if (!is.null(data)){
      
	ta=accrualtime
	tf=followuptime
	ttot=ta+tf
	t<-data[,1]
      colnames(data)<-c("Time","Status","X")
Time=data[,1]
Status=data[,2]
X=data[,3]
cat("test here1..","\n")
      f=smcure(Surv(Time, Status)~X,~X,data=data,model="ph",Var=FALSE)
	time<-sort(t[Status==1])
	#status<-data[,2]
	
	beta0nocure <- coxph(Surv(Time, Status)~X,method="breslow", data=data)$coef
  
	death_point <- sort(unique(subset(Time, Status==1)))

	coxexp <- exp(beta0nocure*X)

	lambda <- numeric()
    	event <- numeric()
      for(i in 1: length(death_point)){
       event[i] <- sum(Status*as.numeric(Time==death_point[i]))
        temp <- sum(as.numeric(Time>=death_point[i])*Status*drop(coxexp))
       	temp1 <- event[i]
       lambda[i] <- temp1/temp
        }
    HHazard <- numeric()
    for(i in 1:length(Time)){
        HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
        if(Time[i]>max(death_point))HHazard[i] <- Inf
        if(Time[i]<min(death_point))HHazard[i] <- 0
        }
 	  snocure <- exp(-HHazard)

	beta0 <- f$beta
	print(beta0)
	gamma0 <- -f$b[2]
	pi0=1-exp(f$b[1])/(1 + exp(f$b[1]))
    s=sort(f$s[Status==1],decreasing = TRUE)
	snocure <- sort(snocure[Status==1],decreasing = TRUE)
     f0<-diff(s)
    
      f0nocure <- diff(snocure)

     s1 <- sum(-f0*as.numeric(time<=tf)[-1])
        s1nocure <- sum(-f0nocure*as.numeric(time<=tf)[-1])
	sc=(ta+tf-time)/ta
	  s2 <- sum(-diff(s)*sc[-length(sc)]*as.numeric(time>tf)[-1])
        s2nocure <- sum(-diff(snocure)*sc[-length(sc)]*as.numeric(time>tf)[-1])
	Spop=pi0+(1-pi0)*s
	m=(gamma0/beta0-log(s))*pi0/Spop-1
	  s3 <- sum(-diff(s)*m[-length(m)]*as.numeric(time<=tf)[-1])
	Spop4=pi0+(1-pi0)*s
	m4=(gamma0/beta0-log(s))*pi0/Spop4-1
	  s4 <- sum(-diff(s)*m4[-length(m4)]*sc[-length(sc)]*as.numeric((time>tf) & (time<=ttot))[-1])
	nonpar=ceiling((qnorm(power)-qnorm(alpha/2))^2*(s1+s2)/((s3+s4)^2*p*(1-p)*(1-pi0)*beta0^2))
	n$nonpar<- nonpar
	n$HR <- exp(beta0) 
	n$OR <- exp(gamma0) 
	n$pi0<- pi0
	cat("\n")
	cat("======================================================================== \n")
	cat("SAMPLE SIZE CALCULATION FOR PH MIXTURE CURE MODEL AND STANDARD PH MODEL \n")
	cat("======================================================================== \n")
	cat("PH Mixture Cure Model with KM estimators: n =",nonpar,"\n")
	pdeathNonpar <- s1+s2
	 pdeathNonpar <- s1nocure+s2nocure
	#cat("Probability of Death: p =",pdeathNonpar,"\n")
	nonparPH<- ceiling((qnorm(power)-qnorm(alpha/2))^2/(p*(1-p)*beta0nocure^2*pdeathNonpar))
	cat("PH Standard Model with KM estimators: n =",nonparPH,"\n")
	n$nonparPH<- nonparPH	
  }
n
 }


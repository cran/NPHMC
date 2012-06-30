NPHMC <-
function(power=0.9,alpha=0.05,accrualtime=3,followuptime=4,p=0.5,hazard=1,accrualdist=c("uniform","increasing","decreasing"),
		hazardratio=0.5,oddsratio=2.15,pi0=0.1,survdist=c("exp","weib"),k=1, data=NULL){
	n<-list()
	class(n) <- c("NPHMC")  
if (is.null(data)){
	beta0 <- log(hazardratio)
	gamma0 <- log(oddsratio)
	  i1 <- integrate(f1,0,followuptime,hazard,survdist,k)$value
	  i2 <- integrate(f2,followuptime,(accrualtime+followuptime),accrualtime,followuptime,hazard,accrualdist,survdist,k)$value	  
	  i3 <- integrate(f3,0,followuptime,hazard,beta0,gamma0,pi0,survdist,k)$value
	  i4 <- integrate(f4,followuptime,(accrualtime+followuptime),accrualtime,followuptime,hazard,accrualdist,beta0,gamma0,pi0,survdist,k)$value
	nsize <- ceiling((qnorm(power)-qnorm(alpha/2))^2*(i1+i2)/((i3+i4)^2*p*(1-p)*(1-pi0)*beta0^2))
	cat("\n")
	cat("======================================================================== \n")
	cat("SAMPLE SIZE CALCULATION FOR PH MIXTURE CURE MODEL AND STANDARD PH MODEL \n")
	cat("======================================================================== \n")
	cat("PH Mixture Cure Model: n =",nsize,"\n")
	pdeath <- i1+i2
	cat("Probability of Death: p =",pdeath,"\n")
	nsizeph <- ceiling((qnorm(power)-qnorm(alpha/2))^2/(p*(1-p)*beta0^2*pdeath))
	cat("PH Standard Model: n =",nsizeph,"\n")
	n$nsize <- nsize 
     }
  if (!is.null(data)){
	ta=accrualtime
	tf=followuptime
	ttot=ta+tf
	t<-data[,1]
      colnames(data)<-c("Time","Status","X","Z")
      f=smcure(Surv(Time, Status)~X,~Z,data=data,model="ph",Var=FALSE)
	time<-sort(t[f$Status==1])
	status<-data[,2]
	beta0 <- f$beta
	gamma0 <- -f$b[2]
	pi0=1-exp(f$b[1])/(1 + exp(f$b[1]))
    s=sort(f$s[f$Status==1],decreasing = TRUE)
     f0<-diff(s)
     s1 <- sum(-f0*as.numeric(time<=tf)[-1])
	sc=(ta+tf-time)/ta
	  s2 <- sum(-diff(s)*sc[-length(sc)]*as.numeric(time>tf)[-1])
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
	cat("Probability of Death: p =",pdeathNonpar,"\n")
	nonparPH<- ceiling((qnorm(power)-qnorm(alpha/2))^2/(p*(1-p)*beta0^2*pdeathNonpar))
	cat("PH Standard Model with KM estimators: n =",nonparPH,"\n")
	n$nonparPH<- nonparPH	
  }
 }


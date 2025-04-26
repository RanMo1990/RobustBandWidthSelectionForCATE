## -----------------------------------------------------------------------------
## Title: R-function EstTauWithBwSelection()
## -----------------------------------------------------------------------------
## Authors: Ran Mo, Honglang Wang
## -----------------------------------------------------------------------------
## Description:  
## This function estimates treatment effects (tau) using different loss functions 
## (L2, density power, Huber, and Tukey) with bandwidth selection through 
## cross-validation. It computes the treatment effect (tau) along with its 
## variance, bias, and confidence intervals based on bootstrapped samples.
## -----------------------------------------------------------------------------
## Reference:  
## BiomJ article: "Outlier Resistant Inference for Conditional
##                 Average Treatment Effect"
## -----------------------------------------------------------------------------
## Required Packages: stats, MASS, KernSmooth, robustbase
## -----------------------------------------------------------------------------
## Usage:  
## EstTauWithCV(
##   X0 = 0,
##   X,
##   D,
##   Y,
##   Hlist=seq(0.01,0.1,length.out =10),
##   HlistL2=seq(0.01,0.5,length.out =10),
##   rep = 1,
##   rhotype = "Huber",
##   bw_type="Plug-in"
## )
## -----------------------------------------------------------------------------
## Arguments:
##
##  X0             A vector of values for x1 at which tau is estimated.
##
##  X              A matrix of dimension n x p, where n is the number of 
##                 observations and p is the number of covariates. The first 
##                 column must be the conditioned covariate X1.
##
##  D              A vector of treatment indicators.
##
##  Y              A vector of responses.
##
##  Hlist          The list of bandwidth used in the cross validation except
##                 the L2 loss (Since the bandwidths chosen for L2 is usually 
##                 greater under contaminated models)
##
##  HlistL2        The list of bandwidth used in the cross validation for
##                 the L2 loss
##
##  rep            Number of bootstrap replications for confidence interval 
##                 estimation.
##
##  rhotype        Specifies the loss function for cross-validation. Options 
##                 include "L1", "L2", "Huber", and "Tukey".
##
##  bw_type        Specifies the choice of bandwidth seletion method. Options
##                 include "CrossValidation", "Plug-in"
##
## -----------------------------------------------------------------------------
## Value:  
## A matrix containing the foL2owing results:
##
## - `X0`: Values of the conditioned covariate X1.
## - `tauL2`, `tauDP`, `tauHuber`, `tauTukey`: Estimated treatment effects using 
##   L2, density power, Huber, and Tukey loss functions, respectively.
## - `mu1` and `mu0`: Estimated means for treated and control groups, respectively.
## - `EstVarL2`, `EstVarDP`, `EstVarHuber`, `EstVarTukey`: Estimated variances for 
##   treatment effects under L2, density power, Huber, and Tukey loss functions.
## - Bootstrap Confidence intervals:
##   - `lbL295`, `ubL295`, `lbDP95`, `ubDP95`, `lbM95`, `ubM95`, `lbMM95`, `ubMM95`: 
##     95% confidence intervals for tau using different loss functions.
##   - `lbL290`, `ubL290`, `lbDP90`, `ubDP90`, `lbM90`, `ubM90`, `lbMM90`, `ubMM90`: 
##     90% confidence intervals for tau using different loss functions.
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
EstTauWithBwSelection<-function(X0=0,X,D,Y,Hlist,HlistL2,rep=1,rhotype="L1",bw_type="Plug-in"){
  SimOneLoop1<-function(X,D,Y,Hlist,method="DP",rhotype="L1"){
    CrossValidation1<-function(X,D,Y,Hlist,OutlierType="N",method="L2",rhotype="L1"){
      rho<-function(u,rhotype="L1"){
        if(rhotype=="L1"){
          Out=abs(u)  
        }
        if(rhotype=="L2"){
          Out=u^2  
        }
        if(rhotype=="Huber"){
          c=1.547
          Out=u^2/2*(abs(u)<c)+(c*abs(u)-c^2/2)*(abs(u)>=c)
        }
        if(rhotype=="Tukey"){
          c=4.685
          Out=1/6*(1-(1-u^2)^3)*(abs(u)<=c)+1/6*(abs(u)>=c)    
        }
        return(Out)
      }
      OneTimeFit<-function(Xminus1,X,D,Y,gma=0.5,bw=0.04,bwNon=0.04,bwM=0.04,bwMM=0.04,method="L2"){
        Step1RIPW<-function(Xpred,X,D){
          X=data.frame(X)
          g=glm(D ~ .,data=X,family = "binomial")
          coe=g$coefficients
          piHat=predict(g,newdata = Xpred,type = "response")
          return(c(coe,piHat))
        }
        
        if(is.vector(Xminus1)==TRUE){
          X0=Xminus1[1]  
        }else{
          X0=Xminus1[,1]
        }
        
        X=data.frame(X)
        Step=Step1RIPW(X,X,D)
        dimX=dim(X)[2]
        StepAlpha=Step[(1:(dimX+1))]
        StepPiHat=Step[-(1:(dimX+1))]
        
        #truncate extreme pihats
        alpha=0.00
        D=subset(D,StepPiHat<=(1-alpha) & StepPiHat>=(alpha))
        X=subset(X,StepPiHat<=(1-alpha) & StepPiHat>=(alpha))
        Y=subset(Y,StepPiHat<=(1-alpha) & StepPiHat>=(alpha))
        StepPiHat=subset(StepPiHat,StepPiHat<=(1-alpha) & StepPiHat>=(alpha))
        
        X1=subset(X[,1],D==1)
        Y1=subset(Y,D==1)
        piHat1=subset(StepPiHat,D==1)
        
        X2=subset(X[,1],D==0)
        Y2=subset(Y,D==0)
        piHat2=subset(StepPiHat,D==0)
        
        PdIpw<-function(x0=X0,gma=0.5,bw=0.04,bwNon=0.04,bwM=0.04,bwMM=0.04){
          ##Non          
            h1=bwNon
            K1=dnorm((X1-x0)/h1)
            l1=locCteSmootherC(X1, Y1, x0, h1, kernel=gaussK, weig = 1/piHat1)
            mu1=l1$beta0
            mu1LocalLinear=l1$beta0
            mu1LL=mu1LocalLinear
            mu1Out=mu1LL  
          
          
          
          ##PD
          if(method=="DP"){
            h1=bw
            K1=dnorm((X1-x0)/h1)
            #Function:calcualte Power Density
            PowerDensity<-function(Y,mu,sd,gma){
              return(dnorm(Y,mean=mu,sd=sd)^gma)
            }
            mu1=mu1LL
            for(j in 1:1000){
              sd1=1/0.675*median(abs(Y1-mu1))
              H1=PowerDensity(Y1,mu1,sd1,gma)
              mu1new=sum(H1*Y1*K1/piHat1)/sum(H1*K1/piHat1)
              
              convi<- abs(mu1new-mu1)
              done <- (convi <= 1e-7)
              
              mu1=mu1new
              if(done) break
            }
            mu1pdIPW=mu1
            mu1Out=mu1
          }        
          
          ##M
          if(method=="Huber"){
            h1=bwM
            K1=dnorm((X1-x0)/h1)
            mu1=mu1LL
            for(j in 1:1000){
              e=Y1-mu1
              sigmahat=mad(e)/0.6745
              sd1=sigmahat
              c=1.547*sigmahat
              u=e
              
              Psi=(u*(abs(u)<=c)-c*(u<(-c))+c*(u>c))
              dPsidu=(abs(u)<=c)
              Wx=1:length(e)
              for (i in 1:length(e)) {
                if(u[i]==0){
                  Wx[i]=dPsidu[i]
                }else{
                  Wx[i]=Psi[i]/u[i]    
                }
              }
              
              
              
              
              lm1=lm(Y1~1,weights=K1/piHat1*Wx)
              mu1new=lm1$coefficients
              
              convi<- abs(mu1new-mu1)
              done <- (convi <= 1e-7)
              
              mu1=mu1new
              if(done) break
            }
            mu1M=lm1$coefficients
            
            mu1Out=mu1M
          }
          
          #Tukey
          if(method=="Tukey"){
            h1=bwMM
            K1=dnorm((X1-x0)/h1)
            
            Mn1=rlm(Y1~1,weights = K1/piHat1,psi = psi.bisquare,method="MM",scale.est="MAD",wt.method = "case")
            sd1=Mn1$s
            
            mu1MM=Mn1$coefficients[1]
            mu1Out=mu1MM
          }    
          
          
          
          
          return(mu1Out)
        }
        l=lapply(X=X0,PdIpw,gma=gma,bw=bw,bwNon=bwNon,bwM=bwM,bwMM=bwMM)
        
        Xminus1=data.frame(t(Xminus1))
        colnames(Xminus1) <- paste0("V", seq_len(ncol(Xminus1)))
        OutIPW=Step1RIPW(Xminus1,X,D)[-(1:5)]
        
        return(cbind(l,OutIPW))
      }
      N=length(D)
      optH=999999999999999999999999
      minRCV=999999999999999999999999
      for (j in 1:length(Hlist)) {
        bw=Hlist[j]
        rholist=NULL
        for (i in 1:N) {
          Xt=X[-i,]
          Dt=D[-i]
          Yt=Y[-i]
          result=OneTimeFit(Xminus1=X[i,],X=Xt,D=Dt,Y=Yt,gma=0.5,bw=bw,bwNon=bw,bwM=bw,bwMM=bw,method=method)
          mu1hat=unlist(result[,1])
          pihat=unlist(result[,2])
          rholist=c(rholist,rho(u=Y[i]-mu1hat,rhotype=rhotype)*D[i]/pihat)  
        }
        RCV=sum(rholist)/N
        if(is.na(RCV)){
          RCV=9999999999999999999999999
        }
        if(RCV<minRCV){
          minRCV=RCV
          optH=bw
        }
      }
      return(c(optH,minRCV))
    }
    #CrossValidation
    CVresult=CrossValidation1(X,D,Y,Hlist,method=method,rhotype=rhotype)
    opth1order=match(CVresult[1],Hlist)
    if(opth1order==1){
      HlistNew=seq(Hlist[opth1order],(Hlist[opth1order+1]+Hlist[opth1order])/2,length.out =6)
    }
    if(opth1order==length(Hlist)){
      HlistNew=seq((Hlist[opth1order]+Hlist[opth1order-1])/2,Hlist[opth1order],length.out =6)
    }
    if((opth1order<length(Hlist))*(opth1order>1)){
      HlistNew=seq((Hlist[opth1order-1]+Hlist[opth1order])/2,(Hlist[opth1order+1]+Hlist[opth1order])/2,length.out =11)
    }
    
    
    CVresult=CrossValidation1(X,D,Y,HlistNew,method=method,rhotype=rhotype)
    
    h=CVresult[1]
    print(paste("method:",method,"bandwidth:",h))
    #
    return(h)
  }
  SimOneLoop0<-function(X,D,Y,Hlist,method="DP",rhotype="L1"){
    CrossValidation0<-function(X,D,Y,Hlist,method="L2",rhotype="L1"){
      rho<-function(u,rhotype="L1"){
        if(rhotype=="L1"){
          Out=abs(u)  
        }
        if(rhotype=="L2"){
          Out=u^2  
        }
        if(rhotype=="Huber"){
          c=1.547
          Out=u^2/2*(abs(u)<c)+(c*abs(u)-c^2/2)*(abs(u)>=c)    
        }
        if(rhotype=="Tukey"){
          c=4.685
          Out=1/6*(1-(1-u^2)^3)*(abs(u)<=c)+1/6*(abs(u)>=c)    
        }
        return(Out)
      }
      OneTimeFit0<-function(X0=24,X,D,Y,gma=0.5,bw=0.04,bwNon=0.04,bwM=0.04,bwMM=0.04,method="L2"){
        Step1RIPW<-function(X,D){
          X=data.frame(X)
          g=glm(D ~ .,data=X,family = "binomial")
          coe=g$coefficients
          piHat=predict(g,newdata = X,type = "response")
          return(c(coe,piHat))
        }
        Step=Step1RIPW(X,D)
        dimX=dim(X)[2]
        StepAlpha=Step[(1:(dimX+1))]
        StepPiHat=Step[-(1:(dimX+1))]
        
        alpha=0.00
        D=subset(D,StepPiHat<=(1-alpha) & StepPiHat>=(alpha))
        X=subset(X,StepPiHat<=(1-alpha) & StepPiHat>=(alpha))
        Y=subset(Y,StepPiHat<=(1-alpha) & StepPiHat>=(alpha))
        StepPiHat=subset(StepPiHat,StepPiHat<=(1-alpha) & StepPiHat>=(alpha))
        
        X1=subset(X[,1],D==0)
        Y1=subset(Y,D==0)
        piHat1=1-subset(StepPiHat,D==0)
        
        
        PdIpw<-function(x0=X0,gma=0.5,bw=0.04,bwNon=0.04,bwM=0.04,bwMM=0.04){
          ##Non
            h1=bwNon
            K1=dnorm((X1-x0)/h1)
            l1=locCteSmootherC(X1, Y1, x0, h1, kernel=gaussK, weig = 1/(piHat1))
            mu1=l1$beta0
            mu1LocalLinear=l1$beta0
            mu1LL=mu1LocalLinear
            mu1Out=mu1LL  

          
          ##PD
          if(method=="DP"){
            h1=bw
            K1=dnorm((X1-x0)/h1)
            PowerDensity<-function(Y,mu,sd,gma){
              return(dnorm(Y,mean=mu,sd=sd)^gma)
            }
            mu1=mu1LL
            for(j in 1:1000){
              sd1=1/0.675*median(abs(Y1-mu1))
              H1=PowerDensity(Y1,mu1,sd1,gma)
              mu1new=sum(H1*Y1*K1/piHat1)/sum(H1*K1/piHat1)
              
              convi<- abs(mu1new-mu1)
              done <- (convi <= 1e-7)
              
              mu1=mu1new
              if(done) break
            }
            mu1pdIPW=mu1
            mu1Out=mu1
          }        
          
          ##M
          if(method=="Huber"){
            h1=bwM
            K1=dnorm((X1-x0)/h1)
            mu1=mu1LL
            for(j in 1:1000){
              e=Y1-mu1
              sigmahat=mad(e)/0.6745
              sd1=sigmahat
              c=1.547*sigmahat
              u=e
              
              Psi=(u*(abs(u)<=c)-c*(u<(-c))+c*(u>c))
              dPsidu=(abs(u)<=c)
              Wx=1:length(e)
              for (i in 1:length(e)) {
                if(u[i]==0){
                  Wx[i]=dPsidu[i]
                }else{
                  Wx[i]=Psi[i]/u[i]    
                }
              }
              
              
              
              
              lm1=lm(Y1~1,weights=K1/piHat1*Wx)
              mu1new=lm1$coefficients
              
              convi<- abs(mu1new-mu1)
              done <- (convi <= 1e-7)
              
              mu1=mu1new
              if(done) break
            }
            mu1M=lm1$coefficients
            
            mu1Out=mu1M
          }
          
          #MM
          if(method=="Tukey"){
            #MM
            h1=bwMM
            K1=dnorm((X1-x0)/h1)
            
            Mn1=rlm(Y1~1,weights = K1/piHat1,psi = psi.bisquare,method="MM",scale.est="MAD",wt.method = "case")
            sd1=Mn1$s
            
            mu1MM=Mn1$coefficients[1]
            mu1Out=mu1MM
          }    
          
          
          return(mu1Out)
        }
        l=lapply(X=X0,PdIpw,gma=gma,bw=bw,bwNon=bwNon,bwM=bwM,bwMM=bwMM)
        
        OutIPW=Step1RIPW(X,D)[-(1:5)]
        
        return(cbind(l,OutIPW))
      }
      N=length(D)
      optH=9999999999999999999999999999
      minRCV=9999999999999999999999999999
      for (j in 1:length(Hlist)) {
        bw=Hlist[j]
        rholist=NULL
        for (i in 1:N) {
          X0=X[i,1]
          Xt=X[-i,]
          Dt=D[-i]
          Yt=Y[-i]
          result=OneTimeFit0(X0=X0,X=Xt,D=Dt,Y=Yt,gma=0.5,bw=bw,bwNon=bw,bwM=bw,bwMM=bw,method=method)
          mu1hat=unlist(result[,1])
          pihat=unlist(result[,2])
          rholist=c(rholist,rho(u=Y[i]-mu1hat,rhotype=rhotype)*(1-D[i])/(1-pihat))  
        }
        RCV=sum(rholist)/N
        if(is.na(RCV)){
          RCV=99999999999999999999999999999999
        }
        if(RCV<minRCV){
          minRCV=RCV
          optH=bw
        }
      }
      return(c(optH,minRCV))
    }
    #CrossValidation
    CVresult=CrossValidation0(X,D,Y,Hlist,method=method,rhotype=rhotype)
    opth1order=match(CVresult[1],Hlist)
    if(opth1order==1){
      HlistNew=seq(Hlist[opth1order],(Hlist[opth1order+1]+Hlist[opth1order])/2,length.out =6)
    }
    if(opth1order==length(Hlist)){
      HlistNew=seq((Hlist[opth1order]+Hlist[opth1order-1])/2,Hlist[opth1order],length.out =6)
    }
    if((opth1order<length(Hlist))*(opth1order>1)){
      HlistNew=seq((Hlist[opth1order-1]+Hlist[opth1order])/2,(Hlist[opth1order+1]+Hlist[opth1order])/2,length.out =11)
    }
    
    
    CVresult=CrossValidation0(X,D,Y,HlistNew,method=method,rhotype=rhotype)
    
    h=CVresult[1]
    print(paste("method:",method,"bandwidth:",h))
    return(h)
  }
  ##Bandwidth Selection
  #Hlist <- seq(0.04, 1, length.out = 10)
  #HlistL2 <- seq(0.02, 4, length.out = 10)
  #bw1=SimOneLoop1(X,D,Y,Hlist,method="DP",rhotype=rhotype)
  
  
  Step1RIPW<-function(X,D){
    X=data.frame(X)
    g=glm(D ~ .,data=X,family = "binomial")
    coe=g$coefficients
    piHat=predict(g,newdata = X,type = "response")
    return(c(coe,piHat))
  }
  Step=Step1RIPW(X,D)
  length(Y)
  length(D)
  length(Step)
  
  PlugInBWForL2<-function(x=X,t=D,y=Y,boundarytrim=0.05,numDesignPoint=100){
    #Generate a sequence of equal distant points from the lowest x to the highest x.
    n=length(x[,1])
    x0=seq(from = x[n*boundarytrim,1], to = x[n-n*boundarytrim+1,1], length.out=numDesignPoint)
    x0width=x0[2]-x0[1]
    Step1RIPW<-function(X,D){
      X=data.frame(X)
      g=glm(D ~ .,data=X,family = "binomial")
      coe=g$coefficients
      piHat=predict(g,newdata = X,type = "response")
      return(c(coe,piHat))
    }
    phiHat=Step1RIPW(x,t)[-(1:5)]
    #Step1: 
    #1.1, estiamte the second derivative of mu1(x_1) on a sequence of equal distant points.
    x1=subset(x,t==1)
    y1=subset(y,t==1)
    phiHat1=subset(phiHat,t==1)
    bw1=dpill(x1[,1],y1,trim = 0.05)
    
    regresult=locCuadSmootherC(x1[,1],y1,x0,bw=bw1,weig = 1/phiHat1,kernel=gaussK)
    Mu1dd=regresult$beta2
    #bw=dpill(x[,1],y*t/phiHat)
    #locpoly(x[,1], y*t/phiHat, bandwidth = bw,drv =2,gridsize = length(x0), range.x = range(x0))
    denx1=density(x[,1], from = min(x0), to = max(x0), n = length(x0))$y
    intMu1dd=sum(Mu1dd^2*denx1*x0width)
    #Step2:
    #2.1, estiamte mu1
    #library(np)
    mu1=locCteSmootherC(x1[,1],y1,x0,bw=bw1,weig = 1/phiHat1,kernel=gaussK)$beta0
    #plot(x1[,1],y1)
    #points(x0,mu1$beta0,type="l")
    #2.2, estiamte sigma^2(x_1)
    bw=bw1#dpill(x[,1],(y-mu1)^2*t/phiHat^2)
    regresult=locCteSmootherC(x[,1], (y-mu1)^2*t/phiHat^2,x0,bw=bw,kernel=gaussK)[,2]
    #locpoly(x[,1], (y-mu1)^2*t/phiHat^2, bandwidth = bw,
    #       gridsize = length(x0), range.x = range(x0))$y
    intSigma2=sum(regresult*x0width)
    #Step3:          
    hopt=(intSigma2/(2*sqrt(pi)*intMu1dd*n))^(1/5)
    return(hopt)
  }
  

PlugInBWForRobust<-function(x=X,t=D,y=Y,boundarytrim=0.1,numDesignPoint=100,psiType="L2"){
  #Generate a sequence of equal distant points from the lowest x to the highest x.
  n=length(x[,1])
  x1sorted=x[order(x[,1]), ]
  x0=seq(from = x1sorted[n*boundarytrim,1], to = x1sorted[n-n*boundarytrim+1,1], length.out=numDesignPoint)
  x0width=x0[2]-x0[1]
  Step1RIPW<-function(X,D){
    X=data.frame(X)
    g=glm(D ~ .,data=X,family = "binomial")
    coe=g$coefficients
    piHat=predict(g,newdata = X,type = "response")
    return(c(coe,piHat))
  }
  piHat=Step1RIPW(x,t)[-(1:5)]
  #Step1: 
  #1.1, estiamte the second derivative of mu1(x_1) on a sequence of equal distant points.
  x1=subset(x,t==1)
  y1=subset(y,t==1)
  piHat1=subset(piHat,t==1)

  
  divisor=20
  blockmax=3
  b=max(x1[,1])
  a=min(x1[,1])
  n1=length(y1)
  Nmax <- max(min(floor(n1/divisor), blockmax), 1)
  Nval <- max(KernSmooth:::cpblock(x1[,1], y1, Nmax, 4),2)
  print(paste("Nval=",Nval))
  out <- KernSmooth:::blkest(x1[,1], y1, Nval, 4)
  sigsqQ <- out$sigsqe
  th24Q <- out$th24e
  gamseh <- (sigsqQ * (b - a)/(abs(th24Q) * n1))
  if (th24Q < 0) 
    gamseh <- (3 * gamseh/(8 * sqrt(pi)))^(1/7)
  if (th24Q > 0) 
    gamseh <- (15 * gamseh/(16 * sqrt(pi)))^(1/7)
  bw1=gamseh
  #print(bw1)
  print(psiType)
  if(psiType=="L2"){
    #regresult=locPolSmootherC(x1[,1], y1, x0, bw=bw1, deg = 3, kernel = gaussK, weig = 1/piHat1)#
    #Mu1dd=regresult$beta2
    Mu1dd=1:numDesignPoint
    for (j in 1:numDesignPoint) {
      K1=dnorm((x1[,1]-x0[j])/bw1)
      Mu1dd[j]=lm(y1~poly(x1[,1], degree = 3, raw = TRUE),weights = K1/piHat1)$coefficients[3] 
    }
    #plot(x0,Mu1dd,type = "l",col="red")
    #points(x0,-25*cos(5*x0),type = "l",col="blue")
  }else if(psiType=="Huber"){
    Mu1dd=1:numDesignPoint
    for (j in 1:numDesignPoint) {
      K1=dnorm((x1[,1]-x0[j])/bw1)
      Mu1dd[j]=rlm(y1~poly(x1[,1], degree = 3, raw = TRUE),weights = K1/piHat1,psi = psi.bisquare,method="M",scale.est="MAD",wt.method = "case")$coefficients[3] 
    }
    #plot(x0,Mu1dd,type = "l",col="red")
    #points(x0,-25*cos(5*x0),type = "l",col="blue")
  }else if(psiType=="Tukey"){
    Mu1dd=1:numDesignPoint
    for (j in 1:numDesignPoint) {
      K1=dnorm((x1[,1]-x0[j])/bw1)
      Mu1dd[j]=rlm(y1~poly(x1[,1], degree = 3, raw = TRUE),weights = K1/piHat1,psi = psi.bisquare,method="MM",scale.est="MAD",wt.method = "case")$coefficients[3]
    }
  }
  
  denx1=density(x[,1], from = min(x0), to = max(x0), n = length(x0))$y
  intMu1dd=sum(Mu1dd^2*denx1*x0width)
  
  #estiamte mu1
  bw=dpill(x1[,1],y1,trim = 0.01)
  if(psiType=="L2"){
    mu1=locCteSmootherC(x=x1[,1],y=y1,xeval=x[,1],bw=bw,kernel=gaussK,weig = 1/piHat1)$beta0

  }else if(psiType=="Huber"){
    mu1=1:n
    for (j in 1:n) {
      K1=dnorm((x1[,1]-x[j,1])/bw)
      mu1[j]=rlm(y1~1,weights = K1/piHat1,psi = psi.bisquare,method="M",scale.est="MAD",wt.method = "case")$coefficients  
    }

  }else if(psiType=="Tukey"){
    mu1=1:n
    for (j in 1:n) {
      K1=dnorm((x1[,1]-x[j,1])/bw)
      mu1[j]=rlm(y1~1,weights = K1/piHat1,psi = psi.bisquare,method="MM",scale.est="MAD",wt.method = "case")
    }
 
  }
  mu1=as.numeric(mu1)
  
  Psi<-function(psitype="L2",u,sd1){
    if(psitype=="L2"){
      Out=u  
    }
    if(psitype=="Huber"){
      c=1.547*sd1
      Out=u*(abs(u)<c)+c*(u>=c)-c*(u<=-c)    
    }
    if(psitype=="Tukey"){
      c=4.685*sd1
      Out=u*(1-(u/c)^2)^2*(abs(u)<=c)
    }
    return(Out)
  }
  dPsi<-function(psitype="L2",u,sd1){

    if(psitype=="L2"){
      
      Out=rep(1,length(u)) 
    }
    if(psitype=="Huber"){
      c=1.547*sd1
      Out=1*(abs(u)<c) 
    }
    if(psitype=="Tukey"){
      c=4.685*sd1
      Out=(5*u^4-6*c^2*u^2+c^4)/c^4*(abs(u)<=c)
    }
    return(Out)
  }
  
  sd1=mad(subset(y-mu1,t==1))/0.6745
  bw=dpill(x[,1],Psi(u=y-mu1,psitype=psiType,sd1)^2*t/piHat^2,trim = 0.01)
  Numer=locCteSmootherC(x[,1], Psi(u=y-mu1,psitype=psiType,sd1)^2*t/piHat^2,x0,bw=bw,kernel=gaussK)[,2]
  bw=dpill(x[,1],dPsi(u=y-mu1,psitype=psiType,sd1)*t/piHat,trim = 0.01)
  Deno=locCteSmootherC(x[,1], dPsi(u=y-mu1,psitype=psiType,sd1)*t/piHat,x0,bw=bw,kernel=gaussK)[,2]^2
  regresult=Numer/Deno
  
  intSigma2=sum(regresult*x0width)
  intSigma2
  
  hopt=(intSigma2/(2*sqrt(pi)*intMu1dd*n))^(1/5)
  return(hopt)
}
if(bw_type=="Plug-in"){
  bwNon1=PlugInBWForRobust(x=X,t=D,y=Y,psiType = "L2")
  bwM1=PlugInBWForRobust(x=X,t=D,y=Y,psiType = "Huber")
  bwMM1=PlugInBWForRobust(x=X,t=D,y=Y,psiType = "Tukey")
}else if(bw_type=="CrossValidation"){
  bwNon1=SimOneLoop1(X,D,Y,Hlist,method="Huber",rhotype=rhotype)
  bwM1=SimOneLoop1(X,D,Y,Hlist,method="Huber",rhotype=rhotype)
  bwMM1=SimOneLoop1(X,D,Y,Hlist,method="Tukey",rhotype=rhotype)
}
  print(paste("finish bw selection1",",bwL21",bwNon1,",bwHuber1",bwM1,"bwTukey1",bwMM1))
  
  bwNon0=bwNon1
  bwM0=bwM1
  bwMM0=bwMM1
  print(paste("finish bw selection0",",bwL20",bwNon0,",bwHuber0",bwM0,"bwTukey0",bwMM0))
  ##Fit the mean
  gma=0.5
  Tau=EstTau(x0=X0,X,D,Y,gma=gma,bw1=10,bwNon1=bwNon1,bwM1=bwM1,bwMM1=bwMM1,bw0=10,bwNon0=bwNon0,bwM0=bwM0,bwMM0=bwMM0)

  TT=cbind(Tau)
  colnames(TT)=c("X0","tauLL","tauDP","tauM","tauMM","tau1LL","tau1DP","tau1M","tau1MM","tau0LL","tau0DP","tau0M","tau0MM","EstVarLL","EstVarDP","EstVarM","EstVarMM")
  rownames(TT)=X0
  
  return(TT)
  
}


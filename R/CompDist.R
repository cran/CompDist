rm(list=ls())


#install.packages("MASS", repos='http://cran.us.r-project.org')
#library(MASS)
#install.packages("numDeriv", repos='http://cran.us.r-project.org')
#library(numDeriv)
#install.packages("fExtremes", repos='http://cran.us.r-project.org')
#library(fExtremes)
#install.packages("actuar", repos='http://cran.us.r-project.org')
#library(actuar)
#install.packages("flexsurv", repos='http://cran.us.r-project.org')
#library(flexsurv)
#install.packages("Renext", repos='http://cran.us.r-project.org')
#library(Renext)
#install.packages("rmutil", repos='http://cran.us.r-project.org')
#library(rmutil)
#install.packages("VGAM", repos='http://cran.us.r-project.org')
#library(VGAM)
#install.packages("PearsonDS", repos='http://cran.us.r-project.org')
#library(PearsonDS)



#------------------  Density Function  -----------------


dcomp<-function(xx,dists,par,borders,par.pos,buffer=c(0,0)){

f<-function(y,string,par){
  if (string=="lnorm"||string=="f"||string=="burr"||string=="gompertz"||string=="levy"||string=="pareto"||string=="gpareto"||string=="lomax"||string=="gamma"||string=="exp"||string=="weibull")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2]))}
  if (string=="chisq"||string=="rayleigh")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1]))}
  if (string=="gev"||string=="loglap")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2], par[3]))}
  if (string=="beta"||string=="pearson")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2], par[3], par[4]))}
  yy=eee(y)
  return(yy)
}


F<-function(y,string,par){
  if (string=="lnorm"||string=="f"||string=="burr"||string=="gompertz"||string=="levy"||string=="pareto"||string=="gpareto"||string=="lomax"||string=="gamma"||string=="exp"||string=="weibull")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1], par[2]))}
  if (string=="chisq"||string=="rayleigh")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1]))}
  if (string=="gev"||string=="loglap")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1], par[2], par[3]))}
  if (string=="beta"||string=="pearson")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1], par[2], par[3], par[4]))}
  yy=eee(y)
  return(yy)
}


if(missing(par.pos)){par.pos<-rep(1,length(dists)-1)}

if(!missing(borders)){
  n<-length(dists)
  theta<-par[[1]]
  phi<-append(1,par[[2]])

  if(is.unsorted(theta)){
    stop()
    print("Theta not increasing")
  }

  if(any(phi<0) | any(phi==0)){
    stop()
    print("Non-positive Phi")
  }

  par[[1]]<-NULL
  par[[1]]<-NULL

  for(i in 1:(n-1)){
      rt<-function(p){
        par[[i+1]][par.pos[i]]<-p
        fi<-function(x){f(x,dists[i],par[[i]])}

        fi1<-function(x){f(x,dists[i+1],par[[i+1]])}
        return(grad(fi,theta[i+1])*fi1(theta[i+1])/fi(theta[i+1])-grad(fi1,theta[i+1]))
      }
      par[[i+1]][par.pos[i]]<-uniroot(rt,lower=theta[i]+buffer[1],upper=theta[i+1]-buffer[2])$root

      c1<-F(theta[i+1],dists[i],par[[i]])-F(theta[i],dists[i],par[[i]])
      c2<-F(theta[i+2],dists[i+1],par[[i+1]])-F(theta[i+1],dists[i+1],par[[i+1]])

      phi[i+1]<-phi[i]*f(theta[i+1],dists[i],par[[i]])/f(theta[i+1],dists[i+1],par[[i+1]])*c2/c1
#    cat(theta[i]," ",theta[i+1]," ",theta[i+2],"    ","\n")
#    cat(dists[i]," ",dists[i+1],"\n")
#    cat(par[[i]]," ",par[[i+1]],"\n")
#    cat(phi[i]," ",phi[i+1],"\n")
#    cat(c1," ",c2,"\n")
   }

  }

  if(missing(borders)){
    theta<-par[[1]]
    phi<-append(1,par[[2]])

    if(is.unsorted(theta)){
      stop()
      print("Theta not increasing")
    }

    if(any(phi<0) | any(phi==0)){
      stop()
      print("Non-positive Phi")
    }

    par[[1]]<-NULL
    par[[1]]<-NULL
  }

  pdf<-function(y){

    eval<-function(x){

      eta<-1/(sum(phi))

      sec<-which.max(theta>x)-1

      res=phi[sec]*(f(x,dists[sec],par[[sec]]))/(F(theta[sec+1],dists[sec],par[[sec]])-F(theta[sec],dists[sec],par[[sec]]))

      res=res*eta
      return(res)
    }
    return(as.numeric(unlist(lapply(y,eval))))
  }

  return(pdf(xx))
}



# ================================================================
#   
#   EXAMPLE 1
# 
# ================================================================
#   
# par<-list()
# distvec<-c("lnorm","gamma")
# par[[1]]<-c(0,1,Inf)
# par[[2]]<-c(1)
# par[[3]]<-c(0,1)
# par[[4]]<-c(1,1)
# 
# 
# 
# x<-seq(0,3,0.01)
# # non-continuous case
# y1<-dcomp(x,distvec,par)
# # continuous case
# y2<-dcomp(x,distvec,par,borders=list(c(0.00001,10)),buffer=c(10e-5,0))
# 
# 
# 
# par(mfrow=c(1,2),oma=rep(0,4))
# xrange<-range(x)
# yrange<-range(y1,y2)
# plot(x,y1,type="l",xlab="x",ylab="Density function",xlim=xrange,ylim=yrange)
# abline(v=1)
# plot(x,y2,type="l",xlab="x",ylab="Density function",xlim=xrange,ylim=yrange)
# abline(v=1)



#------------------  Cumulative Distribution Function  -----------------


pcomp<-function(xx,dists,par,borders,par.pos,buffer=c(0,0)){



f<-function(y,string,par){
  if (string=="lnorm"||string=="f"||string=="burr"||string=="gompertz"||string=="levy"||string=="pareto"||string=="gpareto"||string=="lomax"||string=="gamma"||string=="exp"||string=="weibull")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2]))}
  if (string=="chisq"||string=="rayleigh")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1]))}
  if (string=="gev"||string=="loglap")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2], par[3]))}
  if (string=="beta"||string=="pearson")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2], par[3], par[4]))}
  yy=eee(y)
  return(yy)
}


F<-function(y,string,par){
  if (string=="lnorm"||string=="f"||string=="burr"||string=="gompertz"||string=="levy"||string=="pareto"||string=="gpareto"||string=="lomax"||string=="gamma"||string=="exp"||string=="weibull")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1], par[2]))}
  if (string=="chisq"||string=="rayleigh")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1]))}
  if (string=="gev"||string=="loglap")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1], par[2], par[3]))}
  if (string=="beta"||string=="pearson")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1], par[2], par[3], par[4]))}
  yy=eee(y)
  return(yy)
}


if(missing(par.pos)){par.pos<-rep(1,length(dists)-1)}

  if(!missing(borders)){
    n<-length(dists)
    theta<-par[[1]]
    phi<-append(1,par[[2]])

    if(is.unsorted(theta)){
      stop()
      print("Theta not increasing")
    }

    if(any(phi<0) | any(phi==0)){
      stop()
      print("Non-positive Phi")
    }

    par[[1]]<-NULL
    par[[1]]<-NULL

    for(i in 1:(n-1)){
        rt<-function(p){
          par[[i+1]][par.pos[i]]<-p
          fi<-function(x){f(x,dists[i],par[[i]])}

          fi1<-function(x){f(x,dists[i+1],par[[i+1]])}
          return(grad(fi,theta[i+1])*fi1(theta[i+1])/fi(theta[i+1])-grad(fi1,theta[i+1]))
        }

        par[[i+1]][par.pos[i]]<-uniroot(rt,lower=theta[i]+buffer[1],upper=theta[i+1]-buffer[2])$root

        c1<-F(theta[i+1],dists[i],par[[i]])-F(theta[i],dists[i],par[[i]])
        c2<-F(theta[i+2],dists[i+1],par[[i+1]])-F(theta[i+1],dists[i+1],par[[i+1]])

        phi[i+1]<-phi[i]*f(theta[i+1],dists[i],par[[i]])/f(theta[i+1],dists[i+1],par[[i+1]])*c2/c1

#    cat(theta[i]," ",theta[i+1]," ",theta[i+2],"    ","\n")
#    cat(dists[i]," ",dists[i+1],"\n")
#    cat(par[[i]]," ",par[[i+1]],"\n")
#    cat(phi[i]," ",phi[i+1],"\n")
#    cat(c1," ",c2,"\n")
    }  
  }

  if(missing(borders)){
    theta<-par[[1]]
    phi<-append(1,par[[2]])

    if(is.unsorted(theta)){
      stop()
      print("Theta not increasing")
    }

    if(any(phi<0) | any(phi==0)){
      stop()
      print("Non-positive Phi")
    }

    par[[1]]<-NULL
    par[[1]]<-NULL
  }

  cdf<-function(y){

    eval<-function(x){

      eta<-1/(sum(phi))

      sec<-which.max(theta>x)-1

      if(sec==1){
        res=phi[sec]*(F(x,dists[sec],par[[sec]])-F(theta[sec],dists[sec],par[[sec]]))/(F(theta[sec+1],dists[sec],par[[sec]])-F(theta[sec],dists[sec],par[[sec]]))
      }
      if(sec>1){
        res=sum(phi[1:sec-1])+phi[sec]*(F(x,dists[sec],par[[sec]])-F(theta[sec],dists[sec],par[[sec]]))/(F(theta[sec+1],dists[sec],par[[sec]])-F(theta[sec],dists[sec],par[[sec]]))
      }

      res=res*eta
      return(res)
    }
    return(as.numeric(unlist(lapply(y,eval))))
  }

  return(cdf(xx))
}



# ================================================================
#   
#   EXAMPLE 2
# 
# ================================================================
  
#   
# par<-list()
# distvec<-c("lnorm","gamma")
# par[[1]]<-c(0,1,Inf)
# par[[2]]<-c(1)
# par[[3]]<-c(0,1)
# par[[4]]<-c(1,1)
# 
# 
# 
# x<-seq(0,3,0.01)
# # non-continuous case
# y1<-pcomp(x,distvec,par)
# # continuous case
# y2<-pcomp(x,distvec,par,borders=list(c(0.00001,10)),buffer=c(10e-5,0))
# 
# 
# 
# 
# par(mfrow=c(1,2),oma=rep(0,4))
# xrange<-range(x)
# yrange<-range(y1,y2)
# plot(x,y1,type="l",xlab="x",ylab="Distribution function",xlim=xrange,ylim=yrange)
# abline(v=1,lty=2)
# plot(x,y2,type="l",xlab="x",ylab="Distribution function",xlim=xrange,ylim=yrange)
# abline(v=1,lty=2)


#------------------  Quantile Function  -----------------

qcomp<-function(xx,dists,par,borders,par.pos,buffer=c(0,0)){



f<-function(y,string,par){
  if (string=="lnorm"||string=="f"||string=="burr"||string=="gompertz"||string=="levy"||string=="pareto"||string=="gpareto"||string=="lomax"||string=="gamma"||string=="exp"||string=="weibull")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2]))}
  if (string=="chisq"||string=="rayleigh")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1]))}
  if (string=="gev"||string=="loglap")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2], par[3]))}
  if (string=="beta"||string=="pearson")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2], par[3], par[4]))}
  yy=eee(y)
  return(yy)
}


F<-function(y,string,par){
  if (string=="lnorm"||string=="f"||string=="burr"||string=="gompertz"||string=="levy"||string=="pareto"||string=="gpareto"||string=="lomax"||string=="gamma"||string=="exp"||string=="weibull")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1], par[2]))}
  if (string=="chisq"||string=="rayleigh")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1]))}
  if (string=="gev"||string=="loglap")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1], par[2], par[3]))}
  if (string=="beta"||string=="pearson")
  eee<-function(x){do.call(paste("p",string,sep=""),list(x, par[1], par[2], par[3], par[4]))}
  yy=eee(y)
  return(yy)
}



Finv<-function(y,string,par){
  if (string=="lnorm"||string=="f"||string=="burr"||string=="gompertz"||string=="levy"||string=="pareto"||string=="gpareto"||string=="lomax"||string=="gamma"||string=="exp"||string=="weibull")
  eee<-function(x){do.call(paste("q",string,sep=""),list(x, par[1], par[2]))}
  if (string=="chisq"||string=="rayleigh")
  eee<-function(x){do.call(paste("q",string,sep=""),list(x, par[1]))}
  if (string=="gev"||string=="loglap")
  eee<-function(x){do.call(paste("q",string,sep=""),list(x, par[1], par[2], par[3]))}
  if (string=="beta"||string=="pearson")
  eee<-function(x){do.call(paste("q",string,sep=""),list(x, par[1], par[2], par[3], par[4]))}
  yy=eee(y)
  return(yy)
}



if(missing(par.pos)){par.pos<-rep(1,length(dists)-1)}

  #  Smooth Distribution

  if(!missing(borders)){
    n<-length(dists)
    theta<-par[[1]]
    phi<-append(1,par[[2]])

    if(is.unsorted(theta)){
      stop()
      print("Theta not increasing")
    }

    if(any(phi<0) | any(phi==0)){
      stop()
      print("Non-positive Phi")
    }

    par[[1]]<-NULL
    par[[1]]<-NULL

    for(i in 1:(n-1)){
        rt<-function(p){
          par[[i+1]][par.pos[i]]<-p
          fi<-function(x){f(x,dists[i],par[[i]])}

          fi1<-function(x){f(x,dists[i+1],par[[i+1]])}
#          cat(fi(theta[i+1]),"\n")
#          cat(fi1(theta[i+1]),"\n")
          return(grad(fi,theta[i+1])*fi1(theta[i+1])/fi(theta[i+1])-grad(fi1,theta[i+1]))
        }

        borders[[i]][1]=borders[[i]][1]+buffer[1]
        borders[[i]][2]=borders[[i]][2]-buffer[2]
        
        par[[i+1]][par.pos[i]]<-uniroot(rt,interval=borders[[i]])$root

        c1<-F(theta[i+1],dists[i],par[[i]])-F(theta[i],dists[i],par[[i]])
        c2<-F(theta[i+2],dists[i+1],par[[i+1]])-F(theta[i+1],dists[i+1],par[[i+1]])

        phi[i+1]<-phi[i]*f(theta[i+1],dists[i],par[[i]])/f(theta[i+1],dists[i+1],par[[i+1]])*c2/c1

    }

  n<-length(par)
  q<-pcomp(head(theta,-1),dists,c(list(theta),list(phi[-1]),par),borders=borders,buffer=buffer)

  }

  #  Non-Smooth

  if(missing(borders)){
    theta<-par[[1]]
    phi<-append(1,par[[2]])

    if(is.unsorted(theta)){
      stop()
      print("Theta not increasing")
    }

    if(any(phi<0) | any(phi==0)){
      stop()
      print("Non-positive Phi")
    }

    par[[1]]<-NULL
    par[[1]]<-NULL

    n<-length(par)
    q<-pcomp(head(theta,-1),dists,c(list(theta),list(phi[-1]),par))
    
  }

  #  Return function

  qdf<-function(y){

    eval<-function(x){
      if(x==0){return(min(theta))}
      if(x==1){return(max(theta))}

      if(x<1){

        eta<-1/(sum(phi))
        sec<-max(which(q<x))

        if(sec==1){res=Finv(x*F(theta[sec+1],dists[sec],par[[sec]])/eta,dists[sec],par[[sec]])}
        if(sec>1){res=Finv((x/eta-sum(phi[1:sec-1]))*(F(theta[sec+1],dists[sec],par[[sec]])-F(theta[sec],dists[sec],par[[sec]]))/phi[sec]+F(theta[sec],dists[sec],par[[sec]]),dists[sec],par[[sec]])}

        return(res)
      }
    }
    return(as.numeric(unlist(lapply(y,eval))))
  }

  return(qdf(xx))
}


# ================================================================
#   
#   EXAMPLE 3
# 
# ================================================================
#   
# par<-list()
# distvec<-c("lnorm","gamma")
# par[[1]]<-c(0,1,Inf)
# par[[2]]<-c(1)
# par[[3]]<-c(0,1)
# par[[4]]<-c(1,1)
# 
# 
# 
# x<-seq(0.01,0.99,0.01)
# # non-continuous case
# y1<-qcomp(x,distvec,par)
# # continuous case
# y2<-qcomp(x,distvec,par,borders=list(c(0.00001,10)),buffer=c(10e-5,0))
# 
# 
# 
# 
# par(mfrow=c(1,2),oma=rep(0,4))
# xrange<-range(x)
# yrange<-range(y1,y2)
# plot(x,y1,type="l",xlab="x",ylab="Quantile function",xlim=xrange,ylim=yrange)
# abline(h=1,lty=2)
# plot(x,y2,type="l",xlab="x",ylab="Quantile function",xlim=xrange,ylim=yrange)
# abline(h=1,lty=2)




#------------------  Random Sample Generator  -----------------

rcomp<-function(nn,dists,par,borders,par.pos,buffer=c(0,0)){

if(missing(par.pos)){par.pos<-rep(1,length(dists)-1)}
if(!missing(borders)){qc<-function(x){qcomp(x,dists,par,borders,par.pos,buffer)}}
if(missing(borders)){qc<-function(x){qcomp(x,dists,par)}}

  rdf<-function(y){
    eval<-function(x){
      qs<-runif(x,0,1)
      sample<-qc(qs)
      return(sample)
    }
    return(as.numeric(unlist(lapply(y,eval))))
  }

  return(rdf(nn))
}




# ================================================================
#   
#   EXAMPLE 4
# 
# ================================================================

# par<-list()
# distvec<-c("lnorm","gamma")
# par[[1]]<-c(0,1,Inf)
# par[[2]]<-c(1)
# par[[3]]<-c(0,1)
# par[[4]]<-c(1,1)
# 
# 
# 
# n<-1000
# # non-continuous case
# y1<-rcomp(n,distvec,par)
# # continuous case
# y2<-rcomp(n,distvec,par,borders=list(c(0.00001,10)),buffer=c(10e-5,0))
# 
# 
# 
# 
# par(mfrow=c(1,2),oma=rep(0,4))
# hist(y1,nclass=10,xlab="x",ylab="Frequency",main="")
# hist(y2,nclass=10,xlab="x",ylab="Frequency",main="")


#------------------  Data Fitting Function  -----------------


par.fit<-function(data,dists,par,borders,par.pos,optit=5,buffer=c(0,0),cont=FALSE){



f<-function(y,string,par){
  if (string=="lnorm"||string=="f"||string=="burr"||string=="gompertz"||string=="levy"||string=="pareto"||string=="gpareto"||string=="lomax"||string=="gamma"||string=="exp"||string=="weibull")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2]))}
  if (string=="chisq"||string=="rayleigh")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1]))}
  if (string=="gev"||string=="loglap")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2], par[3]))}
  if (string=="beta"||string=="pearson")
  eee<-function(x){do.call(paste("d",string,sep=""),list(x, par[1], par[2], par[3], par[4]))}
  yy=eee(y)
  return(yy)
}


if(missing(par.pos)){par.pos<-rep(1,length(dists)-1)}

  fit<-function(data){

    N<-length(par)
    sec.length<-numeric(0)
    for(i in 1:N){
      sec.length[i]<-length(par[[i]])
    }

    theta<-par[[1]]
    upper.bound<-tail(theta,n=1)
    lower.bound<-theta[1]

    phi<-par[[2]]
    n<-length(theta)-1
    par[[1]]<-NULL
    par[[1]]<-NULL
    
    for(a in 1:optit){
      
      
      # Weight optimization
      
      PhiLL<-function(p){
        param<-c(list(theta),list(p),par)
        if(exists("borders")){comp.dist<-dcomp(data,dists,param,borders,par.pos,buffer=buffer)}
        if(!exists("borders")){comp.dist<-dcomp(data,dists,param)}
        
        R<-comp.dist
        return(-sum(log(R)))
      }
      
      phi<-nlm(PhiLL,phi)$estimate

      
      #section par optimization
      
      for(i in 1:n){

        d<-subset(data,data>theta[i] & data<theta[i+1])
        LL<-function(param){
          R<-f(d,dists[i],param)
          return(-sum(log(R)))
        }

        par[[i]]=optim(par[[i]],LL,method="Nelder-Mead")$par
      }
      
      
      # theta optimization
      
      ThetaLL<-function(t){
        th<-c(lower.bound,t,upper.bound)
        param<-c(list(th),list(phi),par)
        
        if(exists("borders")){comp.dist<-dcomp(data,dists,param,borders,par.pos,buffer=buffer)}
        if(!exists("borders")){comp.dist<-dcomp(data,dists,param)}
        
        R<-comp.dist
        return(-sum(log(R)))
      }
      
      theta<-theta[-1]
      theta<-theta[-length(theta)]

      if(! cont){thetamid<-nlm(ThetaLL,theta)$estimate}
      if(cont){thetamid<-optim(theta,ThetaLL,control=list(maxit=5))$par}

      theta<-c(lower.bound,thetamid,upper.bound)
      
    }

    LL<-function(param){
      if(exists("borders")){comp.dist<-dcomp(data,dists,param,borders,par.pos,buffer=buffer)}
      if(!exists("borders")){comp.dist<-dcomp(data,dists,param)}
      R<-comp.dist
      return(-sum(na.omit(log(R))))
    }

    #optim(par)
    param<-c(list(theta),list(phi),par)

    LogL<-LL(param)

    k<-N-1
    m<-length(data)

    #AIC
    AIC<-2*(k-LogL)

    #BIC
    BIC<-k*log(m)-2*LogL

    #CAIC
    CAIC<-k*(log(m)+1)-2*LogL

    #AICc
    AICc<-AIC+2*k*(k+1)/(m-k-1)

    #HQC
    HQC<-2*k*log(log(m))-2*LogL

    output<-list(Parameter=param,LogLikelihood=LogL,AIC=AIC,BIC=BIC,AICc=AICc,CAIC=CAIC,HQC=HQC)
    return(output)
  }

  return(fit(data))
}



# ================================================================
#   
#   EXAMPLE 5
# 
# ================================================================
  
# Generate Random Data

# par<-list()
# distvec<-c("lnorm","gamma")
# par[[1]]<-c(0,1,Inf)
# par[[2]]<-c(1)
# par[[3]]<-c(0,1)
# par[[4]]<-c(1,1)
# 
# n<-1000
# # non-continuous case
# r1<-rcomp(n,distvec,par)
# # continuous case
# r2<-rcomp(n,distvec,par,borders=list(c(0.00001,10)),buffer=c(10e-5,0))
# 
# 
# 
# # Initial Guess
# 
# par<-list()
# distvec<-c("lnorm","gamma")
# par[[1]]<-c(0,1,Inf)
# par[[2]]<-c(1)
# par[[3]]<-c(0,0.5)
# par[[4]]<-c(0.5,1)
# 
# # Fitting
# 
# # non-continuous case
# estimate1<-par.fit(r1,distvec,par,optit=1)
# # continuous case
# estimate2<-par.fit(r2,distvec,par,borders=list(c(0.00001,10)),optit=1,buffer=c(10e-5,0),cont=TRUE)
# 
# 
# x<-seq(0,30,0.01)
# # non-continuous case
# y1<-dcomp(x,distvec,estimate1$Parameter)
# # continuous case
# y2<-dcomp(x,distvec,estimate2$Parameter,borders=list(c(0.00001,10)),buffer=c(10e-5,0))
#  
# par(mfrow=c(1,2),oma=rep(0,4))
# hist(r1,probability=T,breaks=40)
# lines(x,y1,col="red")
# hist(r2,probability=T,breaks=40)
# lines(x,y2,col="red")
#  
# estimate1
# estimate2






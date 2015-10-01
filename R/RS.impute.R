 RS.impute=function (X, T.seq, M = 10, theta1 = c(4, 4, 4, 4), theta2 = c(2,2), prop.sds1 = c(0.5, 0.5, 0.5, 0.5), prop.sds2 = c(0.1,0.1), transform = c(1), N.updates = 2000, exclude = NULL,bridge.plot = F, plot.chain = T)
{
  X1=X
  d=diff(T.seq)
  dd=rep(d/(M),each=M)
  dddd=rep(d,each=M+1)
  N=length(X1)
  endpts=c(1,(1:(N-1))*M+1)
  endpts[length(endpts)]=max(endpts)-1
  endpts2 =c(1,(1:(N-1))*(M+1))
  excl=0
  tt=rep((0:M)/M,N-1)
  lasts=seq(M+1,(N-1)*(M+1),by=M+1)
  firsts=seq(1,(M+1)*(N-1),by=M+1)
  ttt=cumsum(c(T.seq[1],dd))

  sgs1=c("sigma[1]","sigma[1]sqrt(X_t)","sigma[1]X_t")
  prior.list = "NONE"
  buffer0 = c("================================================================")
  buffer1 = c("----------------------------------------------------------------")
  buffer2 = c("................................................................")
  buffer3 = c("...   ...   ...   ...   ...   ...   ...   ...   ...   ...   ... ")
  buffer4 = c("_______________________ Drift Function _________________________")
  buffer5 = c("_____________________ Diffusion Function _______________________")
  buffer6 = c("_____________________ Prior Distributions ______________________")
  buffer7 = c("_______________________ Model/Chain Info _______________________")
  trim <- function (x) gsub("([[:space:]])", "", x)
  Info = c(buffer0, "Data Imputation", buffer0, buffer4, "",
           trim(paste0(body("mu")[2])), buffer5,
           "", sgs1[transform[1]], buffer6,
           "", prior.list, "")
  Info = data.frame(matrix(Info, length(Info), 1))
  colnames(Info) = ""
  print(Info, row.names = FALSE, right = F)

  alg1=function()
  {
      rnrm=cumsum(rnorm((M+1)*(N-1),sd=sqrt(dddd/(M))))
      rnrm=(rnrm-rep(rnrm[firsts],each=M+1))
      rnrm=rnrm-tt*rep((rnrm)[lasts],each=M+1)
  }

  nu=function(y,rnrm)
  {
     ymin=rep(y[-N],each=(M+1))
      ypl=rep(y[-1],each=(M+1))
      Xdot=rnrm+((1-tt)*ymin+(tt)*ypl)
      return(Xdot[-firsts[-1]])
  }


  #-----------------------------------------------------------------------------
    if((transform[1]==1))
  {
    f1=function(x1,sval){return(x1/sval)}
    h1=function(y1,sval){return(y1*sval)}
    likelihood2=function(XX1,theta,sval)
    {
      mu1=(mu(XX1,ttt,theta)/sval[1])[-NN]
      dY1t=diff(f1(XX1,sval[1]))
      return(mu1*dY1t-0.5*(mu1^2)*dd)
    }
  }
  #-----------------------------------------------------------------------------

   if((transform[1]==2))
  {
    f1=function(x1,sval){return(2*sqrt(x1)/sval)}
    h1=function(y1,sval){return((y1*sval/2)^2)}
    likelihood2=function(XX1,theta,sval)
    {
      mu1=(((mu(XX1,ttt,theta)/sval[1]-sval[1]/4)/sqrt(XX1)))[-NN]
      dY1t=diff(f1(XX1,sval[1]))
      return(mu1*dY1t-0.5*(mu1^2)*dd)
    }
  }
  #-----------------------------------------------------------------------------
   if(transform[1]==3)
  {
    f1=function(x1,sval){return(log(x1)/sval)}
    h1=function(y1,sval){return(exp(y1*sval))}
    likelihood2=function(XX1,theta,sval)
    {
      mu1=(mu(XX1,ttt,theta)/(sval[1]*XX1)-sval[1]/2)[-NN]
      dY1t=diff(f1(XX1,sval[1]))
      return(mu1*dY1t-0.5*(mu1^2)*dd)
    }
  }



  k   = rep(0,1)
  kk  = matrix(0,1,N.updates)
  ba1 = rep(1,N-1)
  ba2 = ba1
  ll  = rep(0,N.updates)
  sss = 1:(length(X1)*M)
  ss  = seq(1,length(sss),by=M)
  lastss = seq(M,(M)*(N-1),by=M)
  pers1  = matrix(0,length(theta1),N.updates)
  pers2  = matrix(0,length(theta2),N.updates)
  pers1[,1] = theta1
  pers2[,1] = theta2
  n1 = length(theta1)
  n2 = length(theta2)
  Z1.=alg1()

  X1. = h1(nu(f1(X1,pers2[1,1]),Z1.),pers2[1,1])
  #excl=rep(1,(M+1)*(N-1))
  #if(!is.null(exclude))
  #{
  # for(i in exclude)
  # {
  #   excl[((i-1)*(M+1)+1):((i-1)*(M+1)+M)]=0
  # }
  #}
  #excl = excl[-firsts[-1]]

  #if (is.null(exclude))
  #{
  #  exclude = N + 200
  #}

  pb <- txtProgressBar(1, N.updates, 1, style = 3, width = 56)
  tme = Sys.time()
  NN = length(X1.)

  llold.1 = sum(likelihood2(X1., pers1[, 1], pers2[, 1]))
  ll[1] = llold.1
  tme = Sys.time()

  for (i in 2:N.updates)
  {
    theta1.new=pers1[,i-1]+rnorm(n1,sd=prop.sds1)
    sgnew=(pers2[,i-1]+rnorm(n2,sd=prop.sds2))

    measure1 = sum(dnorm(diff(f1(X1,sgnew[1]))/sqrt(d),log=T)-dnorm(diff(f1(X1,pers2[1,i-1]))/sqrt(d),log=T))

    jac1 =(N-1)*(log(sgnew[1]/pers2[1,i-1]))# -log())
    llnew.1=sum(likelihood2(h1(nu(f1(X1,sgnew[1]),Z1.),sgnew[1]),theta1.new,sgnew))
    rat1=min(exp(-jac1+measure1+llnew.1-llold.1),1)
    u=runif(1)
    indtrue=(rat1>u)
    indfalse=(rat1<=u)
    pers1[,i]=theta1.new*indtrue+pers1[,i-1]*indfalse
    pers2[,i]=sgnew*indtrue+pers2[,i-1]*indfalse
    llold.1=llnew.1*indtrue+llold.1*indfalse
    k[1]=k[1]+(rat1>u)
    kk[1,i]=k[1]/i


    Z1.new=alg1()
    ll1=cumsum(likelihood2(h1(nu(f1(X1,pers2[1,i]),Z1.new),pers2[1,i]),pers1[,i],pers2[,i]))
    ll1=diff(ll1[endpts])
    ll2=cumsum(likelihood2(h1(nu(f1(X1,pers2[1,i]),Z1.),pers2[1,i]),pers1[,i],pers2[,i]))
    ll2=diff(ll2[endpts])

    rats=pmin(exp(ll1-ll2),1)
    us=runif(N-1)

    indexes=c(rep((rats>us),each=M+1))
    #print(length(indexes))
    #print(length(Z1.))

    Z1.=Z1.new*indexes+Z1.*(1-indexes)

    ba1=ba1+(rats>us)

    llold.1=sum(likelihood2(h1(nu(f1(X1,pers2[1,i]),Z1.),pers2[1,i]),pers1[,i],pers2[,i]))
    ll[i]=llold.1

     if(bridge.plot)
    {
      plot(X1[1:50]~T.seq[1:50],type="n",pch="+",col="red",ylim=range(h1(nu(f1(X1,pers2[1,i]),Z1.),pers2[1,i])))
      lines(h1(nu(f1(X1,pers2[1,i]),Z1.),pers2[1,i])~ttt,col="blue")
      points(X1.[endpts]~ttt[endpts],pch=19,col='purple')
      #Sys.sleep(1)
     # abline(v=greens,col="purple")
     #Sys.sleep(0.2)
    }
    if(is.na(pers1[1,i]))
    {
      break;
      #stop("Parameters evaluated as NA!")
    }
    if(sum(is.na(ba1))>0)
    {
      break;
      #stop("BridgesevaluatedasNA!")
    }
    setTxtProgressBar(pb,i)
  }
  close(pb)

  tme.eval = function(start_time) {
    start_time = as.POSIXct(start_time)
    dt = difftime(Sys.time(), start_time, units = "secs")
    format(.POSIXct(dt, tz = "GMT"), "%H:%M:%S")
  }
  tme = tme.eval(tme)
  homo = T
  homo.res = T
  if (sum(round(diff(T.seq) - diff(T.seq)[1], 10) == 0) !=
        length(T.seq) - 1) {
    homo.res = F
  }
  homo = T
  Info2 = c(buffer7, paste0("Chain Updates       : ", N.updates),
            paste0("Time Homogeneous    : ", c("Yes", "No")[2 - homo]),
            paste0("Data Resolution     : ", c(paste0("Homogeneous: dt=",round(max(diff(T.seq)), 4)), paste0("Variable: min(dt)=",round(min(diff(T.seq)), 4), ", max(dt)=",round(max(diff(T.seq)), 4)))[2 - homo.res]), paste0("Imputed Resolution  : ",c(paste0("Homogeneous: dt=", round(max(diff(T.seq)/M),4)), paste0("Variable: min(dt)=", round(min(diff(T.seq)/M),4), ", max(dt)=", round(max(diff(T.seq)/M), 4)))[2 -homo.res]), paste0("Elapsed time        : ",tme), paste0("dim(theta)          : ", round(n1 +n2, 3)), buffer1)
  Info2 = data.frame(matrix(Info2, length(Info2), 1))
  colnames(Info2) = ""
  print(Info2, row.names = FALSE, right = F)

    acc=kk
    nper=n1+n2
    d1=1:((n1+n2)+2)
    d2=d1
    O=outer(d1,d2)
    test=O-((n1+n2)+2)
    test[test<0]=100
    test=test[1:4,1:4]
    wh=which(test==min(test))
    wh
    d1=d1[col(test)[wh[1]]]
    d2=d2[row(test)[wh[1]]]
    par(mfrow=c(d1,d2))
    cols=rainbow(nper)
    j=0
    for(i in 1:n1)
    {
      j=j+1
      plot(pers1[i,],col=cols[j],type="l",ylab="",
      xlab="Iteration",main=paste0("theta[",i,"]"))
      abline(v=min(10000,N.updates/2),lty="dotdash")
    }
    for(i in 1:n2)
    {
     j=j+1
     plot(pers2[i,],col=cols[j],type="l",ylab="",
     xlab="Iteration",main=paste0("sigma[",i,"]"))
     abline(v=min(10000,N.updates/2),lty="dotdash")
    }
   plot(acc[1,],type="n",col="blue",ylim=c(0,1),main="AcceptanceRates",xlab="Iteration",ylab="Rate")
   polygon(c(0,length(acc[1,]),length(acc[1,]),0),c(0.4,0.4,0,-0),col="lightblue",border=NA)
   #polygon(c(0,length(acc[1,]),length(acc[1,]),0),c(1,1,0.6,0.6),col="lightgreen",border=NA)
   lines(acc[1,],type="l",col="darkblue")
   #lines(acc[2,],type="l",col="darkgreen")
   abline(h=seq(0,1,1/10),lty="dotted")
   plot(ba1/N.updates,pch=4,col=4,ylim=c(0,1),type="n",
   main=paste("BBAcceptanceRates(M=",M,")"),xlab="Transition",ylab="Rate")
   polygon(c(0,length(ba1/N.updates),length(ba1/N.updates),0),c(1,1,0.8,0.8),col="lightblue",border=NA)
   points(ba1/N.updates,pch="-",col=4)
   points(ba2/N.updates,pch="-",col=3)
   abline(h=seq(0,1,1/10),lty="dotted")
   abline(h=c(0.6),lty="dashed",col="red",lwd=1.2)
   if(any(ba1/N.updates<0.6))
  {
    wh=which((ba1/N.updates)<0.6)
    text((ba1/N.updates)[wh]~wh,labels=wh,pch=0.5,pos=1)
  }
  if(any(ba2/N.updates<0.6))
 {
    wh=which((ba2/N.updates)<0.6)
    text((ba2/N.updates)[wh]~wh,labels=wh,pch=0.5,pos=1)
 }
 #abline(v=exclude,col="grey75",lty="dashed")
 #mtext(exclude,side=1,cex=0.5,at=exclude)
 theta=rbind(pers1,pers2)
 return(list(per.matrix=t(theta),acceptence.rate=t(acc),bridge.rate=cbind(ba1/N.updates,ba2/N.updates),run.time=tme))
}
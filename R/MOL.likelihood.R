MOL.likelihood=function(X,time,delt,xlims,N=31,theta=c(0),plt=T,wrt=F)
{
  topmatter='
  #include <RcppArmadillo.h>
  #include <math.h>
  #define pi           3.14159265358979323846  /* pi */
  using namespace arma;
  using namespace Rcpp;
  using namespace R;
  // [[Rcpp::depends("RcppArmadillo")]]
  // [[Rcpp::export]]
  vec prod(vec a,vec b)
{
  return(a%b);
}
  '
  m1='
  vec mu(vec x,double t,vec theta)
{
  vec res='
  
  s11=';
  return(res);
}
  vec sig(vec x,double t,vec theta)
{
  vec res=0*x+'
  
  
  
  fsolver=
    ';
  return(res);
}
  // [[Rcpp::export]]
  vec f(vec y,vec x,double t,double hh,int N,vec theta)
{
  vec atemp(N);
  atemp.fill(0);
  atemp(span(1,N-2))= -(mu(x(span(2,N-1)),t,theta)%y(span(2,N-1))
  -mu(x(span(0,N-3)),t,theta)%y(span(0,N-3)))/(2*hh)
  +(sig(x(span(2,N-1)),t,theta)%y(span(2,N-1))
  -2*sig(x(span(1,N-2)),t,theta)%y(span(1,N-2))
  +sig(x(span(0,N-3)),t,theta)%y(span(0,N-3)))/(2*hh*hh);
  return atemp;
}
  // [[Rcpp::export]]
  mat  solver(vec y,vec x,double d,int N,int N2,double delt,double hh,vec theta)
{
  mat x0(N,N);
  mat fx0(N,N);
  mat fx1(N,N);
  mat fx2(N,N);
  mat fx3(N,N);
  mat x1(N,N);
  mat x2(N,N);
  mat x3(N,N);
  x0=y;
  for (int i = 1; i < N2+1; i++)
{
  fx0 =f(x0,x,d,hh,N,theta);
  x1  =x0+delt*(0.5*fx0);
  fx1 =f(x1,x,d+0.5*delt,hh,N,theta);
  x2  =x0+delt*(0.5*fx1);
  fx2 =f(x2,x,d+0.5*delt,hh,N,theta);
  x3  =x0+delt*(fx2);
  fx3 =f(x3,x,d+delt,hh,N,theta);
  x0=x0+(0.1666667*fx0+0.3333333*fx1+0.3333333*fx2+0.1666667*fx3)*delt;
  d=d+delt;
}
  return(x0);
}
  '
  
  txt=paste0(topmatter,m1,body('mu')[2],s11,body('sig')[2],fsolver)
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp(code=txt)
  if(wrt){write(txt,'MOL.likelihood.cpp')}
  M=length(X)
  vals=rep(0,M-1)
  for(i in 1:(M-1))
  {
    x=seq(xlims[1],xlims[2],length=N)
    whx=min(which(abs(x-X[i])==min(abs(x-X[i]))))
    x=x+(X[i]-x[whx])
    z=rep(0,N)
    z[whx]=(1/diff(x)[1])
    Nstop = length(seq(round(time[i],3),round(time[i+1],3),delt))
    res=solver(z,x,time[i],N,Nstop,delt,diff(x)[1],c(0,theta))
    test1x=which(x<=X[i+1])
    test2x=which(x>X[i+1])
    
    whminx = max(test1x)
    whmaxx = min(test2x)
    xl=x[whminx]
    xu=x[whmaxx]
    
    vals[i] = log(res[whminx])+(X[i+1]-x[whminx])/(x[whmaxx]-x[whminx])*(log(res[whmaxx])-log(res[whminx]))
    if(plt)
    {
      par(mfrow=c(1,2))
      plot(X~time,type='l',main='Data',xlab='Time',ylab='State')
      abline(v=(time[i+1]-time[i])/2+time[i],col='red',lwd=2)
      plot(res~x,type='l')
      abline(v=X[i+1],col='grey85',lwd=2,lty='dashed')
      points(exp(vals[i])~X[i+1],pch=21,bg='red')
      axis(1,at=x,labels=NA,tcl=-0.25)
    }
  }
  return(list(likelihood=vals,density=res,Xt=x))
}
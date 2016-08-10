globalVariables(c('mu1','mu2','sig11','sig12','sig21','sig22'))
BiMOL.aic=function(X,time,delt,xlims,ylims,N ,theta , diff.type,plt=TRUE, wrt= FALSE, border = NA)
{
  Y =  X[,2]
  X =  X[,1]
  if(!missing(diff.type))
  {
    #print('yeah')
    sig11 = function(x,y,t){}
    sig22 = function(x,y,t){}
    if(diff.type[1] == 1){body(sig11) = parse(text=paste0('(theta[',length(theta)-1,'])'))}
    if(diff.type[1] == 2){body(sig11) = parse(text=paste0('(theta[',length(theta)-1,']*sqrt(X))'))}
    if(diff.type[1] == 3){body(sig11) = parse(text=paste0('(theta[',length(theta)-1,']*X)'))}
    if(diff.type[2] == 1){body(sig22) = parse(text=paste0('(theta[',length(theta),'])'))}
    if(diff.type[2] == 2){body(sig22) = parse(text=paste0('(theta[',length(theta),']*sqrt(Y))'))}
    if(diff.type[2] == 3){body(sig22) = parse(text=paste0('(theta[',length(theta),']*Y)'))}
    #print(body(sig)[2])
  }
   solver   =function(Xs, Xt, theta, N , delt , N2, tt  , P , alpha, lower , upper, tro  ){}
  rm(list =c('solver'))
  topmatter='
  #include <RcppArmadillo.h>
  #include <math.h>
  #define pi           3.14159265358979323846  /* pi */
  using namespace arma;
  using namespace Rcpp;
  using namespace R;
  // [[Rcpp::depends("RcppArmadillo")]]
  // [[Rcpp::export]]
  mat prod(mat a,mat b)
{
  return(a%b);
}
  '
  m1='
  mat mu1(mat X,mat Y,double t,vec theta)
{
  mat res=0*X+'
  
  m2=';
  return(res);
}
  mat mu2(mat X,mat Y,double t,vec theta)
{
  mat res=0*X+'
  
  s11=';
  return(res);
}
  
  mat sig11(mat X,mat Y,double t,vec theta)
{
  mat res=0*X+'
  s22=';
  return(pow(res,2));
}
  
  mat sig22(mat X,mat Y,double t,vec theta)
{
  mat res=0*Y+'
  
  fsolver=
    ';
  return(pow(res,2));
}
  // [[Rcpp::export]]
  mat f(mat z,mat x,mat y,double t,double hh,int N,vec theta)
{
  mat atemp(N,N);
  atemp.fill(0);
  mat D1 = -mu1(x,y,t,theta)%z;
  mat D2 = -mu2(x,y,t,theta)%z;
  mat D3 = 0.5*sig11(x,y,t,theta)%z;
  mat D4 = 0.5*sig22(x,y,t,theta)%z;
  atemp(span(1,N-2),span(1,N-2)) =
  (D1(span(2,N-1),span(1,N-2))-  D1(span(0,N-3),span(1,N-2)))/(2*hh)
  +(D2(span(1,N-2),span(2,N-1))-  D2(span(1,N-2),span(0,N-3)))/(2*hh)
  +(D3(span(2,N-1),span(1,N-2))-2*D3(span(1,N-2),span(1,N-2)) +D3(span(0,N-3),span(1,N-2)))/(hh*hh)
  +(D4(span(1,N-2),span(2,N-1))-2*D4(span(1,N-2),span(1,N-2)) +D4(span(1,N-2),span(0,N-3)))/(hh*hh);
  
  
  return atemp;
}
  // [[Rcpp::export]]
  mat  solver(mat z,mat x,mat y,double d,int N,int N2,double delt,double hh,vec theta)
{
  mat x0(N,N);
  mat fx0(N,N);
  mat fx1(N,N);
  mat fx2(N,N);
  mat fx3(N,N);
  mat x1(N,N);
  mat x2(N,N);
  mat x3(N,N);
  x0=z;
  for (int i = 1; i < N2+1; i++)
{
  fx0 =f(x0,x,y,d,hh,N,theta);
  x1  =x0+delt*(0.5*fx0);
  fx1 =f(x1,x,y,d+0.5*delt,hh,N,theta);
  x2  =x0+delt*(0.5*fx1);
  fx2 =f(x2,x,y,d+0.5*delt,hh,N,theta);
  x3  =x0+delt*(fx2);
  fx3 =f(x3,x,y,d+delt,hh,N,theta);
  x0=x0+(0.1666667*fx0+0.3333333*fx1+0.3333333*fx2+0.1666667*fx3)*delt;
  d=d+delt;
}
  return(x0);
}
  '
  
  txt=paste0(topmatter,m1,body('mu1')[2],m2,body('mu2')[2],s11,body('sig11')[2],s22,body('sig22')[2],fsolver)
  #library(Rcpp)
  #library(RcppArmadillo)
  sourceCpp(code=txt)
  if(wrt){write(txt,'BiMOL.likelihood.cpp')}
  M=length(X)
  vals=rep(0,M-1)
  if(plt)
  {
    par(mfrow=c(1,1))
  }
  for(i in 1:(M-1))
  {
    x=seq(xlims[1],xlims[2],length=N)
    y=seq(ylims[1],ylims[2],length=N)
    whx=min(which(abs(x-X[i])==min(abs(x-X[i]))))
    why=min(which(abs(y-Y[i])==min(abs(y-Y[i]))))
    x=x+(X[i]-x[whx])
    y=y+(Y[i]-y[why])
    xx=outer(x,rep(1,N))
    yy=outer(rep(1,N),y)
    z=matrix(0,N,N)
    z[whx,why]=(1/diff(x)[1])*(1/diff(y)[1])
    Nstop = length(seq(round(time[i],3),round(time[i+1],3),delt))
    res=solver(z,xx,yy,time[i],N,Nstop,delt,diff(x)[1],c(0,theta))
    test1x=which(x<=X[i+1])
    test2x=which(x>X[i+1])
    test1y=which(y<=Y[i+1])
    test2y=which(y>Y[i+1])
    
    whminx = max(test1x)
    whmaxx = min(test2x)
    whminy = max(test1y)
    whmaxy = min(test2y)
    xl=x[whminx]
    xu=x[whmaxx]
    yl=y[whminy]
    yu=y[whmaxy]
    A = (xu-X[i+1])/(xu-xl)*log(res[whminx,whminy])+(X[i+1]-xl)/(xu-xl)*log(res[whmaxx,whminy])
    B = (xu-X[i+1])/(xu-xl)*log(res[whminx,whmaxy])+(X[i+1]-xl)/(xu-xl)*log(res[whmaxx,whmaxy])
    vals[i] =(yu-Y[i+1])/(yu-yl)*A+(Y[i+1]-yl)/(yu-yl)*B

    if(plt)
    {
      #plot(X~time,col='#222299',type='l',main='Time Series',xlab='Time',ylab='State',ylim=range(c(X,Y)))
      #lines(Y~time,col='#BBCCEE')
      #abline(v=(time[i+1]-time[i])/2+time[i],col='darkblue',lwd=1,lty='dashed')
      ress=persp(x,y,res,col='white',main='Density',box=T, xlab='X_t', ylab='Y_t', zlab='Density', shade=0.8, r = sqrt(0.8),phi=30, theta=145, border = border)
      points(trans3d(X[i+1], Y[i+1], exp(vals[i]), pmat = ress), bg = 'red', pch = 21)
      
      #points(trans3d(x, y, exp(diag(res)), pmat = ress), pch = 21)
      legend('topright',bty='n',pch = 21, pt.bg ='red',legend = substitute(list(X[t[i]],Y[t[i]]),list(i=i+1)))
      Sys.sleep(0.1)
    }
  }
  return(list(AIC = -2*sum(vals)+2*length(theta),likelihood=vals,p = length(theta)))
}
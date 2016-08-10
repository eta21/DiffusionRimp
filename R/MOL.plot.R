MOL.plot=function(x)
{

   if(class(x)=='MOL.density')
  {
    if(requireNamespace('rgl', quietly = TRUE))
    {
      open3d(windowRect=c(50,50,640+50,50+640),zoom=0.95)
      persp3d(x=x$Xt,y=x$time,z=x$density,col='white',box=FALSE,xlab='State (X_t)',ylab='Time (t)',zlab='Density f(X_t|X_s)')
      play3d(spin3d(axis=c(0,0,1), rpm=3), duration=10)
    }else
    {
      persp(x=x$Xt,y=x$time,z=x$density,col='white',xlab='State (X_t)',ylab='Time (t)',zlab='Density f(X_t|X_s)',border=NA,shade=0.5,theta=145)
    }
  }
  
   if(class(x)=='MOL.passage')
  {
    if(requireNamespace('rgl', quietly = TRUE))
    {
      open3d(windowRect=c(50,50,640+50,50+640),zoom=0.95)
      persp3d(x=x$Xt,y=x$time,z=x$surface,col='white',box=FALSE,xlab='State (X_t)',ylab='Time (t)',zlab='Survival Probability')
      play3d(spin3d(axis=c(0,0,1), rpm=3), duration=10)
    }else
    {
      persp(x=x$Xt,y=x$time,z=x$surface,col='white',xlab='X_s',ylab='Time (t)',zlab='Survival Probability',border=NA,shade=0.5,theta=145)
    }
    plot(x$dens~x$time,type='l',col="#222299",main='First passage time density',xlab ='Time (t)',ylab='Density')
  }

  
    if(class(x)=='BiMOL.density')
  {
      colpal=function(n){rev(c(sequential_hcl(n-1,power=0.8,l=c(40,100)),'white'))}
    if(dim(x$density)[3]<=25){ss = 1:dim(x$density)[3]}else{ss =round(seq(1,dim(x$density)[3],length=25))}
    for(i in ss)
    {
      # Now illustrate the density:
      filled.contour(x$Xt,x$Yt,x$density[,,i],
                     main=paste0('Transition Density \n (t = ',x$time[i],')'),
                     color.palette=colpal
                     ,xlab='Xt',ylab='Yt')
    }
  }
  
      if(class(x)=='BiMOL.passage')
  {
    colpal=function(n){rev(sequential_hcl(n,power=0.8,l=c(40,100)))}
    if(dim(x$surface)[3]<=25){ss = 1:dim(x$surface)[3]}else{ss =round(seq(1,dim(x$surface)[3],length=25))}
    for(i in ss)
    {
      # Now illustrate the density:
      filled.contour(x$Xt,x$Yt,x$surface[,,i],
                     main=paste0('Survival Probability \n (t = ',x$time[i],')'),
                     color.palette=colpal
                     ,xlab='Xs',ylab='Ys')
    }
     plot(x$dens~x$time,type='l',col="#222299",main='First passage time density',xlab ='Time (t)',ylab='Density')
  }
}
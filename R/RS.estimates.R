RS.estimates = function(x,thin = 100, burns, CI = c(0.05,0.95), corrmat = FALSE, acf.plot =TRUE , palette = 'mono')
{

   if(class(x)=='RS.impute')
  {
     if(missing(burns)){burns =min(round(dim(x$par.matrix)[1]/2),25000)}
     windw = seq(burns,dim(x$par.matrix)[1],thin)
     est = apply(x$par.matrix[windw,], 2, mean)
     CI=t(apply(x$par.matrix[windw,], 2, quantile,probs = CI))
     form = function(x,mm = 2){format(round(x, mm), nsmall = mm)}
     dat=(cbind(form(cbind(est,CI),3)))
     dat = matrix(as.numeric(as.matrix(dat)),dim(dat)[1])
     rownames(dat)=paste0('theta[',1:dim(x$par.matrix)[2],']')
     colnames(dat) = c('Estimate','Lower_CI','Upper_CI')

     dat2=(form(cor(x$par.matrix[windw,])))
     dat2 = matrix(as.numeric(as.matrix(dat2)),dim(dat2)[1])
     rownames(dat2)=paste0('theta[',1:dim(x$par.matrix)[2],']')
     colnames(dat2)=paste0('theta[',1:dim(x$par.matrix)[2],']')
     if(acf.plot)
     {
     nper=dim(x$par.matrix)[2]
      if(nper==1){par(mfrow=c(1,2))}
      if(nper==2){par(mfrow=c(2,2))}
      if(nper==3){par(mfrow=c(2,2))}
      if(nper>3)
      {
        d1=1:((nper)+1)
        d2=d1
        O=outer(d1,d2)
        test=O-((nper)+1)
        test[test<0]=100
        test=test[1:4,1:4]
        test
        wh=which(test==min(test))

        d1=d1[col(test)[wh[1]]]
        d2=d2[row(test)[wh[1]]]
        par(mfrow=c(d1,d2))
      }
          if(palette=='mono')
    {
      cols =rep('#222299',nper)
    }else{
      cols=rainbow_hcl(nper, start = 10, end = 275,c=100,l=70)
    }
     for(i in 1:dim(x$par.matrix)[2])
     {
        acf(x$par.matrix[windw,i],main=paste0('ACF: theta[',i,']\nThin=',thin,', Burns=',burns,', N=',length(windw)),col = cols[i],lwd=2)
     }
     }
     if(corrmat){return(list(estimates = dat, corrmat = dat2))}
     return(dat)
  }

}

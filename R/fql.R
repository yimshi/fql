
 ginv<-function(X, tol = sqrt(.Machine$double.eps)){
   ## Generalized Inverse of a Matrix
   dnx <- dimnames(X)
   if(is.null(dnx)) dnx <- vector("list", 2)
   s <- svd(X)
   nz <- s$d > tol * s$d[1]
   structure(
     if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
     dimnames = dnx[2:1])
 }
#bw bandwidth
 varest<-function(mi,mu,er,bw){
   n=length(mu);ker=exp(-((mu-mi)/bw)^2/2)
   An0=sum(ker)/(n*bw)
   An1=sum(ker*(mu-mi))/(n*bw)
   An2=sum(ker*((mu-mi)^2))/(n*bw)
   WI=(ker*(An2-(mu-mi)*An1)/(An0*An2-An1^2))/(n*bw)
   sigmu=sum(WI*er)
   return(sigmu)
 }

#flexible quasi-likelihood


fql<-function (formula, data,na.action,preci=0.00001,initial.beta = 'Negative Binomial')
{

  n=dim(data)[1]
  p<-2
  i1<-100
  i2<-100
  bw=(4/3)^(1/5)*(n^(-1/5))
  #m=length(all.vars(formula))-1

  mf<-model.frame(formula=formula, data=data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  options(warn=-1)
  x <- model.matrix(mt, mf, contrasts)
  options(warn=1)
  m=dim(x)[2]-1
  #  x2<-matrix(as.numeric(x),ncol=m+1,nrow=n)
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset))  {
    if (length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }

  datae<-data.frame(matrix(0,nrow=n,ncol=2))
  dataee<-datae
  mu<-c(1:n)*0


  #estall=matrix(0, 1, 10);

  datam<-as.matrix(data)



  ###########################################################################################
  ############################################# Pnql ########################################
  ###########################################################################################

  # Initial step: assume constant vi for all subjects
  if(initial.beta == 'Negative Binomial'){
    fitnb<-glm.nb(formula=formula,link="log",data=data)
    bnew<-summary(fitnb)$coefficients[,1]
  }else {  bnew<-c(1:(m+1))*0}
  bold<-bnew+0.1
  Di<-matrix(0,nrow=n,ncol=m+1); Dp<-matrix(0,nrow=(m+1),ncol=(m+1))
  Dp[row(Dp)==col(Dp)] <- 1; Dp[1:(m+1),1:(m+1)]<-0
  lambda<-0.01; iter2<-1


  while (max(abs(bnew-bold))>preci&iter2<i2)
  {
    bold<-bnew
    #offset
    if (!is.null(offset)){
      eta<-x%*%bold+as.matrix(offset,nrow=length(offset))
    }else {eta<-x%*%bold}
    #link function
    #if(link=="log"){
   mu=exp(eta)
    #} else if(link =="logit") {mu= eta/(1+eta) }

    for (i in 1:n){ Di[i,]<- mu[i]*x[i,] }
    DV<-matrix(0,m+1,m+1);
    for (i in 1:n){DV<-DV+Di[i,]%*%t(Di[i,])}

    DY<-c(1:(m+1))*0;
    for (i in 1:n){DY<-DY+Di[i,]*(y[i]-mu[i])}
    if(qr(DV)$rank==(m+1)){
      bnew<-bold+qr.solve(DV,DY)
    }else{
      bnew<-bold+ginv(DV)%*%DY}
    iter2<-iter2+1
  }

  bet<-bnew

  ########
  # Use spline to estimate vi

  res=try(
    {
      bn<-bet; bo<-bn+0.1; iter1<-1; DVR=1

      while (max(abs(bn-bo))>preci&iter1<i1){
        bo<-bn;
        if (!is.null(offset)){
          eta<-x%*%bo+as.matrix(offset,nrow=length(offset))
        }else {eta<-x%*%bo}
        mu<-exp(eta);er=(y-mu)^2

        datae[,1]<-log((y-mu)^2); datae[,2]<-mu

        # Fit the function of V

        fite2<-gam(X1~s(X2,bs="ps"),family=gaussian(link="identity"),data=datae)
        #may have some error on estimation of v
        et<-sort(mu); newd<-data.frame(X2=mu); v<-exp(predict(fite2,newd));
        bnew=bo #bnew<-bet
        bold<-bnew+0.1; iter2<-1
        Di<-matrix(0,nrow=n,ncol=m+1)
        Dp<-matrix(0,nrow=(m+1),ncol=(m+1))
        Dp[row(Dp)==col(Dp)] <- 1
        Dp[1:(m+1),1:(m+1)]<-0

        # vi's are now different for different subjects, re-estimate beta in Step 2
        while (max(abs(bnew-bold))>preci&iter2<i2){

          bold<-bnew
          #offset
          if (!is.null(offset)){
            eta<-x%*%bold+as.matrix(offset,nrow=length(offset))
          }else {eta<-x%*%bold}
          #link function
          #if(link=="log"){
            mu=exp(eta)
          #} else if(link =="logit") {mu= eta/(1+eta) }
          er=(y-mu)^2

          for (i in 1:n){ Di[i,]<- mu[i]*x[i,]}

          DV<-matrix(0,m+1,m+1)
          for (i in 1:n){DV<-DV+Di[i,]%*%t(Di[i,])/v[i]}

          DY<-c(1:(m+1))*0
          for (i in 1:n){DY<-DY+Di[i,]*(y[i]-mu[i])/v[i]}

          if (qr(DV)$rank==(m+1)){
            bnew<-bold+qr.solve(DV,DY)
          }else{
            bnew<-bold+ginv(DV)%*%DY }
          iter2<-iter2+1
        } # end while in line 122

        bn<-bnew
        bn
        iter1<-iter1+1

      } # end while in line 102

      # Sigma is the center piece in Eqn. (5)
      Sigma<-matrix(0,m+1,m+1);
      for (i in 1:n){Sigma<-Sigma+Di[i,]%*%t(Di[i,])*(y[i]-mu[i])^2/(v[i])^2}
      #sigma=(y-mu)/as.vector(v)/(n-m)
      DVR=qr(DV)$rank
      if(DVR==(m+1)){
        DV.I<-qr.solve(DV)
      }else{
        DV.I<-ginv(DV)
      }
      DV1=DV.I%*%Sigma%*%DV.I # Eqn (5)

      step=0
      coefficients=bn
      se=sqrt(diag(DV1)) #s.e.
      p_value=1-pchisq(coefficients^2/se^2,df=1)     # p-value

      meanerror=mean(er)

      estimation=data.frame(coefficients,se,p_value)
      row.names(estimation)=colnames(x)

      estall=list(estimation,covariance=DV1,meanerror=meanerror,step=step)
      estall$fittedvalues=fite2$fitted.values
      estall$muest=t(mu)
      estall$v=v
    }, silent = F

  ) # end of try in line 101

  # Use another method to estimate vi: kernel method

  if(inherits(res, "try-error")==TRUE|DVR<(m+1))
  {
    res1=try(
      {
        bn<-bet; bo<-bn+0.1; iter1<-1; DVR=1

        while (max(abs(bn-bo))>preci&iter1<i1){
          bo<-bn
          if (!is.null(offset)){
            eta<-x%*%bo+as.matrix(offset,nrow=length(offset))
          }else {eta<-x%*%bo}
          mu<-exp(eta);er=(y-mu)^2; v=mu*0

          for (i in 1:n)  {v[i]=varest(mu[i],mu,er,bw+1)}

          bnew=bo #bnew<-bet
          bold<-bnew+0.1; iter2<-1
          Di<-matrix(0,nrow=n,ncol=m+1)
          Dp<-matrix(0,nrow=(m+1),ncol=(m+1))
          Dp[row(Dp)==col(Dp)] <- 1
          Dp[1:(m+1),1:(m+1)]<-0

          while (max(abs(bnew-bold))>preci&iter2<i2){

            bold<-bnew
            #offset
            if (!is.null(offset)){
              eta<-x%*%bo+as.matrix(offset,nrow=length(offset))
            }else {eta<-x%*%bo}
            #link function
            #if(link=="log"){
              mu=exp(eta)
            #} else if(link =="logit") {mu= eta/(1+eta) }

            for (i in 1:n){Di[i,]<- mu[i]*x[i,]}

            DV<-matrix(0,m+1,m+1)
            for (i in 1:n){DV<-DV+Di[i,]%*%t(Di[i,])/v[i]}

            DY<-c(1:(m+1))*0
            for (i in 1:n){DY<-DY+Di[i,]*(y[i]-mu[i])/v[i]}

            if (qr(DV)$rank==(m+1)){
              bnew<-bold+qr.solve(DV,DY)
            }else{
              bnew<-bold+ginv(DV)%*%DY }
            iter2<-iter2+1
          }

          bn<-bnew
          iter1<-iter1+1
        }

        if (!is.null(offset)){
          eta<-x%*%bn+as.matrix(offset,nrow=length(offset))
        }else {eta<-x%*%bn}
        mu<-exp(eta); er=(y-mu)^2; v=mu*0
        for (i in 1:n){v[i]=varest(mu[i],mu,er,bw+1)}
        Di<-matrix(0,nrow=n,ncol=m+1); for (i in 1:n){Di[i,]<- mu[i]*x[i,]}
        DV<-matrix(0,m+1,m+1);
        for (i in 1:n){DV<-DV+Di[i,]%*%t(Di[i,])/v[i]}

        Sigma<-matrix(0,m+1,m+1);
        for (i in 1:n){Sigma<-Sigma+Di[i,]%*%t(Di[i,])*(y[i]-mu[i])^2/(v[i])^2}
        #sigma=(y-mu)/as.vector(v)/(n-m)
        DVR=qr(DV)$rank
        if(DVR==(m+1)){
          DV.I<-qr.solve(DV)
        }
        else{
          DV.I<-ginv(DV)
        }
        DV1=DV.I%*%Sigma%*%DV.I

        step=1
        coefficients=bn
        se=sqrt(diag(DV1)) #s.e.
        p_value=1-pchisq(coefficients^2/se^2,df=1)     # p-value

        meanerror=mean(er)

        estimation=data.frame(coefficients,se,p_value)
        row.names(estimation)=colnames(x)

        estall=list(estimation,covariance=DV1,meanerror=meanerror,step=step)
        estall$muest=t(mu)
        estall$v=v
      }, silent = TRUE
    ) # end of try in line 172

    # Use bw+3 instead of bw+1 in varest(), and rerun

    if(inherits(res1, "try-error")==TRUE|DVR<(m+1)){
      res2=try(
        {
          bn<-bet; bo<-bn+0.1; iter1<-1; DVR=1

          while (max(abs(bn-bo))>preci&iter1<i1)
          {
            bo<-bn
            if (!is.null(offset)){
              eta<-x%*%bo+as.matrix(offset,nrow=length(offset))
            }else {eta<-x%*%bo}
            mu<-exp(eta); er=(y-mu)^2; v=mu*0


            for (i in 1:n)  {v[i]=varest(mu[i],mu,er,bw+3)}

            bnew=bo #bnew<-bet
            bold<-bnew+0.1; iter2<-1
            Di<-matrix(0,nrow=n,ncol=m+1)
            Dp<-matrix(0,nrow=(m+1),ncol=(m+1))
            Dp[row(Dp)==col(Dp)] <- 1
            Dp[1:(m+1),1:(m+1)]<-0

            while (max(abs(bnew-bold))>preci&iter2<i2)
            {

              bold<-bnew
              #offset
              if (!is.null(offset)){
                eta<-x%*%bold+as.matrix(offset,nrow=length(offset))
              }else {eta<-x%*%bold}
              #link function
              #if(link=="log"){
                mu=exp(eta)
              DV<-matrix(0,m+1,m+1)
              #} else if(link =="logit") {mu= eta/(1+eta) }

              for (i in 1:n){Di[i,]<- mu[i]*x[i,]}

              for (i in 1:n){DV<-DV+Di[i,]%*%t(Di[i,])/v[i]}

              DY<-c(1:(m+1))*0
              for (i in 1:n){DY<-DY+Di[i,]*(y[i]-mu[i])/v[i]}

              if (qr(DV)$rank==(m+1)){
                bnew<-bold+qr.solve(DV,DY)
              }else{
                bnew<-bold+ginv(DV)%*%DY }
              iter2<-iter2+1
            } # end of while in line 255

            bn<-bnew
            iter1<-iter1+1
          } # end of while in line 240

          if (!is.null(offset)){
            eta<-x%*%bn+as.matrix(offset,nrow=length(offset))
          }else {eta<-x%*%bn}
          mu<-exp(eta); er=(y-mu)^2; v=mu*0
          for (i in 1:n){v[i]=varest(mu[i],mu,er,bw+2)}
          Di<-matrix(0,nrow=n,ncol=m+1); for (i in 1:n){Di[i,]<- mu[i]*x[i,]}
          DV<-matrix(0,m+1,m+1);
          for (i in 1:n){DV<-DV+Di[i,]%*%t(Di[i,])/v[i]}
          Sigma<-matrix(0,m+1,m+1);
          for (i in 1:n){Sigma<-Sigma+Di[i,]%*%t(Di[i,])*(y[i]-mu[i])^2/(v[i])^2}
          #sigma=(y-mu)/as.vector(v)/(n-m)
          DVR=qr(DV)$rank
          if(DVR==(m+1)){
            DV.I<-qr.solve(DV)
          }
          else{
            DV.I<-ginv(DV)
          }

          DV1=DV.I%*%Sigma%*%DV.I

          step=2
          coefficients=bn
          se=sqrt(diag(DV1)) #s.e.
          p_value=1-pchisq(coefficients^2/se^2,df=1)     # p-value

          meanerror=mean(er)

          estimation=data.frame(coefficients,se,p_value)
          row.names(estimation)=colnames(x)

          estall=list(estimation,covariance=DV1,meanerror=meanerror,step=step)
          estall$muest=t(mu)
          estall$v=v
        }, silent = TRUE
      )

      # If problem still exists, use vi=1 for all subjects

      if(inherits(res2, "try-error")==TRUE|DVR<(m+1))
      {
        res3=try(
          {
            bn<-bet; bo<-bn+0.1; iter1<-1; DVR=1

            while (max(abs(bn-bo))>preci&iter1<i1)
            {
              bo<-bn
              if (!is.null(offset)){
                eta<-x%*%bo+as.matrix(offset,nrow=length(offset))
              }else {eta<-x%*%bo}
              mu<-exp(eta); er=(y-mu)^2; v=mu*0+1;

              bnew=bo #bnew<-bet
              bold<-bnew+0.1; iter2<-1
              Di<-matrix(0,nrow=n,ncol=m+1)
              Dp<-matrix(0,nrow=(m+1),ncol=(m+1))
              Dp[row(Dp)==col(Dp)] <- 1
              Dp[1:(m+1),1:(m+1)]<-0

              while (max(abs(bnew-bold))>preci&iter2<i2){

                bold<-bnew
                #offset
                if (!is.null(offset)){
                  eta<-x%*%bold+as.matrix(offset,nrow=length(offset))
                }else {eta<-x%*%bold}
                #link function
                #if(link=="log"){
                  mu=exp(eta)
                #} else if(link =="logit") {mu= eta/(1+eta) }

                for (i in 1:n){Di[i,]<- mu[i]*x[i,]}

                DV<-matrix(0,m+1,m+1)
                for (i in 1:n){DV<-DV+Di[i,]%*%t(Di[i,])/v[i]}

                DY<-c(1:(m+1))*0
                for (i in 1:n){DY<-DY+Di[i,]*(y[i]-mu[i])/v[i]}

                if (qr(DV)$rank==(m+1)){
                  bnew<-bold+qr.solve(DV,DY)
                }else{
                  bnew<-bold+ginv(DV)%*%DY }
                iter2<-iter2+1
              }

              bn<-bnew
              iter1<-iter1+1
            } # end of while in line 308

            if (!is.null(offset)){
              eta<-x%*%bn+as.matrix(offset,nrow=length(offset))
            }else {eta<-x%*%b}
            mu<-exp(eta); er=(y-mu)^2; v=mu*0
            for (i in 1:n){v[i]=varest(mu[i],mu,er,bw+2)}
            Di<-matrix(0,nrow=n,ncol=m+1); for (i in 1:n){Di[i,]<- mu[i]*x[i,]}
            DV<-matrix(0,m+1,m+1);
            for (i in 1:n){DV<-DV+Di[i,]%*%t(Di[i,])/v[i]}
            Sigma<-matrix(0,m+1,m+1);
            for (i in 1:n){Sigma<-Sigma+Di[i,]%*%t(Di[i,])*(y[i]-mu[i])^2/(v[i])^2}
            #sigma=(y-mu)/as.vector(v)/(n-m)
            DVR=qr(DV)$rank
            if(DVR==(m+1)){
              DV.I<-qr.solve(DV)
            }
            else{
              DV.I<-ginv(DV)
            }
            DV1=DV.I%*%Sigma%*%DV.I


            step=3
            coefficients=bn
            se=sqrt(diag(DV1)) #s.e.
            p_value=1-pchisq(coefficients^2/se^2,df=1)     # p-value

            meanerror=mean(er)

            estimation=data.frame(coefficients,se,p_value)
            row.names(estimation)=colnames(x)

            estall=list(estimation,covariance=DV1,meanerror=meanerror,step=step)
            estall$muest=t(mu)
            estall$v=v
          }, silent = TRUE
        ) # end of try in line 304
      } # end of if in line 303
    }




  }
  return(estall)
}


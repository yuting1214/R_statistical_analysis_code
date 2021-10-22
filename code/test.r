mmsimstudy1<-function(family="gaussian",n=200,m=10,beta=c(2,0),tau=2,nd=100,balance=TRUE,long=FALSE)
{
  bps <- FALSE  # Indicate perfect separation for binomial data
  print(paste("m = ",m))

  remlx1 <- matrix(NA,nrow=nd,ncol=6)
  mlx1 <- matrix(NA,nrow=nd,ncol=8) # 3: rsq_V (fixed only); 4: rsq_KL (fixed only);7: rsq_V (fixed only for X); 8: rsq_KL (fixed only for X)
  remlx2 <- matrix(NA,nrow=nd,ncol=6)
  mlx2 <- matrix(NA,nrow=nd,ncol=8) # 3: rsq_V (fixed only); 4: rsq_KL (fixed only);7: rsq_V (fixed only for X); 8: rsq_KL (fixed only for X)
  for(d in 1:nd)
  {

    print(paste("  data set = ",d))
    tfamily <- family
    if( family=="probit" )
      tfamily <- "binomial"
    if( long ) # only binomial or gaussian
    {
      simdata <- simldata(family=tfamily,beta=beta,tau=tau,n=n,m=m,balance=balance)
      simdata$yx$seq <- rep(1:(n/m),m)
    }
    else
      simdata <- simglmm(family=tfamily,beta=beta,tau=tau,n=n,m=m,balance=balance)
    
    if( family=="gaussian")
    {
      # ------------------------------------------
      # Y~X1

      # REML:
      if( long )
        remlfit <- lme(y~x1,random=list(subject=pdDiag(~1)),
                     corr=corAR1(form=~seq|subject),data=simdata$yx)
      else
        remlfit <- lmer(y~x1+(1|subject),data=simdata$yx)
      
      # rsq
      trsq <- try(rsq(remlfit))
      if(inherits(trsq,"try-error"))
        trsq <- list(model=NA,fixed=NA, random=NA)

      remlx1[d,1] = trsq$model
      remlx1[d,2] = trsq$fixed
      remlx1[d,3] = trsq$random
      
      # r2_xu
      trsq <- r2_xu(remlfit)
      remlx1[d,4] = trsq
      
      # r2_nakagawa
      trsq <- r2_nakagawa(remlfit)
      if( check_singularity(remlfit) )
        remlx1[d,5] = 0
      else
        remlx1[d,5] = trsq[[1]]
      
      remlx1[d,6] = trsq[[2]]
      
      # ------------------------------------------
      # ML:
      if(long )
        mlfit <- lme(y~x1,random=list(subject=pdDiag(~1)),
                     corr=corAR1(form=~seq|subject),data=simdata$yx,method="ML")
      else
        mlfit <- lmer(y~x1+(1|subject),data=simdata$yx,REML=F)
      
      # rsq
      trsq <- try(rsq(mlfit))
      if(inherits(trsq,"try-error"))
        trsq <- list(model=NA,fixed=NA, random=NA)

      mlx1[d,1] = trsq$model
      mlx1[d,2] = trsq$fixed
      mlx1[d,3] = trsq$random
      
      # r2_xu
      trsq <- r2_xu(mlfit)
      mlx1[d,4] = trsq
      
      # r2_nakagawa
      trsq <- r2_nakagawa(mlfit)
      if( check_singularity(remlfit) )
        remlx1[d,5] = 0
      else
        mlx1[d,5] = trsq[[1]]
      
      mlx1[d,6] = trsq[[2]]


      # ------------------------------------------
      # Y~X2
      
      # REML:
      if( long )
        remlfit <- lme(y~x2,random=list(subject=pdDiag(~1)),
                       corr=corAR1(form=~seq|subject),data=simdata$yx)
      else
        remlfit <- lmer(y~x2+(1|subject),data=simdata$yx)
      
      # rsq
      trsq <- try(rsq(remlfit))
      if(inherits(trsq,"try-error"))
        trsq <- list(model=NA,fixed=NA, random=NA)

      remlx2[d,1] = trsq$model
      remlx2[d,2] = trsq$fixed
      remlx2[d,3] = trsq$random
      
      # r2_xu
      trsq <- r2_xu(remlfit)
      remlx2[d,4] = trsq
      
      # r2_nakagawa
      trsq <- r2_nakagawa(remlfit)
      if( check_singularity(remlfit) )
        remlx1[d,5] = 0
      else
        remlx2[d,5] = trsq[[1]]
      
      remlx2[d,6] = trsq[[2]]
      
      # ------------------------------------------
      # ML:
      if( long )
        mlfit <- lme(y~x2,random=list(subject=pdDiag(~1)),
                       corr=corAR1(form=~seq|subject),data=simdata$yx,method="ML")
      else
        mlfit <- lmer(y~x2+(1|subject),data=simdata$yx,REML=F)
      
      # rsq
      trsq <- try(rsq(mlfit))
      if(inherits(trsq,"try-error"))
        trsq <- list(model=NA,fixed=NA, random=NA)

      mlx2[d,1] = trsq$model
      mlx2[d,2] = trsq$fixed
      mlx2[d,3] = trsq$random
      
      # r2_xu
      trsq <- r2_xu(mlfit)
      mlx2[d,4] = trsq
      
      # r2_nakagawa
      trsq <- r2_nakagawa(mlfit)
      if( check_singularity(remlfit) )
        remlx1[d,5] = 0
      else
        mlx2[d,5] = trsq[[1]]
      
      mlx2[d,6] = trsq[[2]]
    } else { #"binomial" or "probit" or others
      sdata <<- simdata$yx
      # ------------------------------------------
      # Y~X1
      
      # ML only
      if( family=="probit" )
        mlfit <- try(glmer(y~x1+(1|subject),data=sdata,family=binomial(link="probit")),silent=T)
      else
        mlfit <- try(glmer(y~x1+(1|subject),data=sdata,family=family),silent=T)
      
      if(inherits(mlfit,"try-error"))
      {
        mlx1[d,1] = NA
        mlx1[d,2] = NA
        #mlx1[d,3] = NA
        
        # r2_nakagawa
        mlx1[d,5] = NA
        mlx1[d,6] = NA
      }
      else
      {
        trsq <- try(rsq(mlfit))
        if(inherits(trsq,"try-error"))
          trsq <- list(model=NA,fixed=NA, random=NA)

        mlx1[d,1] = trsq$model
        mlx1[d,2] = trsq$fixed
        #mlx1[d,3] = trsq$random
        
        # r2_nakagawa
        trsq <- r2_nakagawa(mlfit)
        if( check_singularity(mlfit) )
          mlx1[d,5] = 0
        else
          mlx1[d,5] = trsq[[1]]
        
        mlx1[d,6] = trsq[[2]]
      }
      
      # Fixed effects only for rsq_V & rsq_KL
      if( family=="probit" )
        mlfit <- glm(y~x1+as.factor(subject),data=sdata,family=binomial(link="probit"))
      else
        mlfit <- glm(y~x1+as.factor(subject),data=sdata,family=family)
      
      mlx1[d,3] = rsq(mlfit)  # rsq_V
      mlx1[d,4] = rsq(mlfit,type="kl")  # rsq_KL

      # Without subject effects:
      if( family=="probit" )
        mlfit <- glm(y~x1,data=sdata,family=binomial(link="probit"))
      else
        mlfit <- glm(y~x1,data=sdata,family=family)
      
      mlx1[d,7] = rsq(mlfit)  # rsq_V
      mlx1[d,8] = rsq(mlfit,type="kl")  # rsq_KL
      
      # ------------------------------------------
      # Y~X2

      if( family=="probit" )
        mlfit <- try(glmer(y~x2+(1|subject),data=sdata,family=binomial(link="probit")),silent=T)
      else
        mlfit <- try(glmer(y~x2+(1|subject),data=sdata,family=family),silent=T)

      if(inherits(mlfit,"try-error"))
      {
        mlx2[d,1] = NA
        mlx2[d,2] = NA
        #mlx2[d,3] = NA
        
        # r2_nakagawa
        mlx2[d,5] = NA
        mlx2[d,6] = NA
      }
      else
      {
        # rsq
        trsq <- try(rsq(mlfit))
        if(inherits(trsq,"try-error"))
          trsq <- list(model=NA,fixed=NA, random=NA)

        mlx2[d,1] = trsq$model
        mlx2[d,2] = trsq$fixed
        #mlx2[d,3] = trsq$random
        
        # r2_nakagawa
        trsq <- r2_nakagawa(mlfit)
        if( check_singularity(mlfit) )
          mlx2[d,5] = 0
        else
          mlx2[d,5] = trsq[[1]]
        
        mlx2[d,6] = trsq[[2]]
      }
      
      # Fixed effects only for rsq_v & rsq_KL
      if( family=="probit" )
        mlfit <- glm(y~x2+as.factor(subject),data=sdata,family=binomial(link="probit"))
      else
        mlfit <- glm(y~x2+as.factor(subject),data=sdata,family=family)
      
      mlx2[d,3] = rsq(mlfit)  # rsq_V
      mlx2[d,4] = rsq(mlfit,type="kl")  # rsq_KL
      
      # Without subject effects:
      if( family=="probit" )
        mlfit <- glm(y~x2,data=sdata,family=binomial(link="probit"))
      else
        mlfit <- glm(y~x2,data=sdata,family=family)
      
      mlx2[d,7] = rsq(mlfit)  # rsq_V
      mlx2[d,8] = rsq(mlfit,type="kl")  # rsq_KL
    }
  }
  
  if( family=="gaussian")
  {
    remlx1.median <- apply(remlx1,2,median,na.rm=T)
    remlx2.median <- apply(remlx2,2,median,na.rm=T)
  }
  
  mlx1.median <- apply(mlx1,2,median,na.rm=T)
  mlx2.median <- apply(mlx2,2,median,na.rm=T)
  
  retList <- list(mlx1.median=mlx1.median,mlx2.median=mlx2.median)
  if( family=="gaussian")
    retList <- list(remlx1.median=remlx1.median,mlx1.median=mlx1.median,
                    remlx2.median=remlx2.median,mlx2.median=mlx2.median)
  
  retList
}

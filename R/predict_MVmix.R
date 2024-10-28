


summarize_pred <- function(pred){
  summlist <- lapply(pred,function(pMat){
    lapply(pMat,function(pmat){
      summ <- data.frame(mean=apply(pmat,1,mean),
                         lower=apply(pmat,1,function(x)quantile(x,0.025)),
                         upper=apply(pmat,1,function(x)quantile(x,0.975)) )
    })
  })
  return(summlist)
}

## when summing all fns, need to add in intercept
summarize_pred_all <- function(pred,obj,include_intercept,newZ){
  summlist <- lapply(1:length(pred),function(kk){
    pmat <- Reduce("+",pred[[kk]])+
      as.numeric(include_intercept)*matrix(obj$b0[,kk],ncol=length(obj$b0[,kk]),nrow=nrow(newZ))+
      newZ%*%t(obj$betaZk[,(kk-1)*obj$const$pz+(1:obj$const$pz)])
    return(data.frame(mean=apply(pmat,1,mean),
                      lower=apply(pmat,1,function(x)quantile(x,0.025)),
                      upper=apply(pmat,1,function(x)quantile(x,0.975)) ))
  })
  return(summlist)
}


## function for estimating fitted/predicted values
predict_MVmix <- function(obj,
                          newX=NULL,
                          newZ=NULL, ## if null, dont include covariate effects
                          fixomega=FALSE,
                          fixomegaval=NULL,
                          include_intercept=TRUE,
                          allx=FALSE){ ## combine all x simultaneously

  # beta <- assign_betas(obj)
  # omega <- assign_omegas(obj)

  I_b0 <- as.numeric(include_intercept)
  if(allx==TRUE){
    I_b0 <- 0 ## so we dont double count
  }

  if(is.null(newX)){
    newX <- obj$const$X
  }
  if(!is.list(newX)){
    newX <- rep(list(newX), obj$const$p)
  }

  if(is.null(newZ)){
    newZ <- matrix(0,nrow=nrow(newX[[1]]),ncol=ncol(obj$const$Zcovariates))#0*obj$const$Zcovariates
  }

  RR <- nrow(as.matrix(obj$sigma2))
  if(fixomega==FALSE){ ## sample theta

    pred <- lapply(1:obj$const$K,function(kk){ ## loop over outcomes
      lapply(1:obj$const$p,function(jj){ ## loop over exposures
        sapply(1:RR,function(rr){## loop over samples
          I_b0*(obj$b0[rr,kk])+
            get_Btheta(newX[[jj]]%*%c(obj$omega[[jj]][[kk]][rr,]),obj$const,
                                          list(Ztheta=matrix(obj$Ztheta[,,rr],
                                                             nrow=obj$const$K)),kk,jj)%*%
            (obj$beta[[jj]][[kk]][rr,])
        })
      })
    })
  }else{ ## otherwise use fixed theta

    if(!is.null(fixomegaval)){ ## use given value
      omega <- fixomegaval
    }else{## otherwise use mean values
      omega <- lapply(omega,function(omeg){
        lapply(omeg,function(om){
          apply(om,2,mean)})})
    }
    for(jj in 1:obj$const$p){ ## standardizing if necessary
      for(kk in 1:obj$const$K){
        omega[[jj]][[kk]] <- omega[[jj]][[kk]]/sqrt(sum(omega[[jj]][[kk]]^2))
      }
    }

    pred <- lapply(1:obj$const$K,function(kk){ ## loop over outcomes
      lapply(1:obj$const$p,function(jj){ ## loop over exposures
        Btheta <- get_Btheta(newX[[jj]]%*%c(omega[[jj]][[kk]]),
                             obj$const,list(Ztheta=matrix(obj$Ztheta[,,1],
                                                          nrow=obj$const$K)),kk,jj)
        sapply(1:RR,function(rr){## loop over samples
          I_b0*(obj$b0[rr,kk])+Btheta%*%(obj$beta[[jj]][[kk]][rr,])
        })
      })
    })

  }

  if(allx==TRUE){
    summpred <- summarize_pred_all(pred,obj,include_intercept,newZ)
  }else{
    summpred <- summarize_pred(pred)
  }
  return(list(values=pred,
              summary=summpred))
}



#' Summarize predictions
#' @keywords internal
summarize_pred <- function(pred,contrast){

  summlist <- lapply(pred,function(pMat){
    lapply(pMat,function(pmat){

      if(contrast==TRUE){
        pmat <- pmat-matrix(c(pmat[floor(nrow(pmat)/2),]),ncol=ncol(pmat),nrow=nrow(pmat),byrow=TRUE)
      }

      summ <- data.frame(mean=apply(pmat,1,mean),
                         lower=apply(pmat,1,function(x)stats::quantile(x,0.025)),
                         upper=apply(pmat,1,function(x)stats::quantile(x,0.975)) )
    })
  })
  return(summlist)
}


#' Summarize all predictions with intercept handling
#' @keywords internal
summarize_pred_all <- function(pred,obj,include_intercept,newZ){
  summlist <- lapply(1:length(pred),function(kk){
    pmat <- Reduce("+",pred[[kk]])+
      as.numeric(include_intercept)*matrix(obj$b0[,kk],ncol=length(obj$b0[,kk]),nrow=nrow(newZ))+
      newZ%*%t(obj$betaZk[,(kk-1)*obj$const$pz+(1:obj$const$pz)])
    return(data.frame(mean=apply(pmat,1,mean),
                      lower=apply(pmat,1,function(x)stats::quantile(x,0.025)),
                      upper=apply(pmat,1,function(x)stats::quantile(x,0.975)) ))
  })
  return(summlist)
}


#' Predict/estimate exposure response surface
#'
#' Function for estimating fitted/predicted values
#'
#' @param obj fitted model object
#' @param newX optional list of new exposure values
#' @param newZ optional matrix of new covariates; if null, dont include covariate effects
#' @param fixomega hold weights fixed (default FALSE)
#' @param fixomegaval fixed weight values (if fixomega=TRUE)
#' @param include_intercept include intercept in predictions (T/F)
#' @param allx combine all x simultaneously or do individual exposure prediction (T/F)
#' @param contrast report contrasts (T/F)
#'
#' @return list of predictions/fitted values
#'
#' @export
predict_MVmix <- function(obj,
                          newX=NULL,
                          newZ=NULL, ## if null, dont include covariate effects
                          fixomega=FALSE,
                          fixomegaval=NULL,
                          include_intercept=TRUE,
                          allx=FALSE,## combine all x simultaneously
                          contrast=FALSE){ ## report contrasts


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
      omega <- lapply(obj$omega,function(omeg){
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
    summpred <- summarize_pred(pred,contrast)
  }
  return(list(values=pred,
              summary=summpred))
}



#' Estimate lagged effects
#'
#' Estimate effect of 1 unit change in x_lag for each different lag holding others constant
#'
#' @param obj fitted model
#' @param Xhold value to compare to and for other exposures to be set to. NULL defaults to median
#' @param Xshift increase in X for contrasts
#' @param fixomega hold weights fixed (default FALSE)
#' @param fixomegaval fixed weight values (if fixomega=TRUE)
#'
#' @return lagged effect estimates
#'
#' @export
est_lag <- function(obj,
                    Xhold=NULL, ## value to compare to and for other exposures to be set to. NULL defaults to median
                    Xshift=1, ## increase in X for contrasts
                    fixomega=FALSE,
                    fixomegaval=NULL){

  newX <- lapply(1:obj$const$p,function(jj){
    if(is.null(Xhold)){
      Xholdjj <- stats::median(obj$const$X[[jj]])
    }else{
      Xholdjj <- Xhold
    }
    ## make all values equal to hold
    newXjj <- matrix(Xholdjj,nrow=obj$const$L,ncol=obj$const$L)
    diag(newXjj) <- diag(newXjj)+Xshift
    ## make first row of
    newXjj <- as.matrix(rbind(rep(Xholdjj,obj$const$L),newXjj))
  })

  if(obj$const$DLM==TRUE){
    I_b0 <- 1
  }else{
    stop("Only built for DLMs")
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
      omega <- lapply(obj$omega,function(omeg){
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


  summlist <- lapply(pred,function(pMat){
      lapply(pMat,function(pmat){


        pmat <- pmat-matrix(c(pmat[1,]),ncol=ncol(pmat),nrow=nrow(pmat),byrow=TRUE)


        summ <- data.frame(mean=apply(pmat[-1,],1,mean),
                           lower=apply(pmat[-1,],1,function(x)stats::quantile(x,0.025)),
                           upper=apply(pmat[-1,],1,function(x)stats::quantile(x,0.975)) )
      })
    })



  return(list(values=pred,
              summary=summlist))
}

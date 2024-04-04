


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


## function for estimating fitted/predicted values
predict_MVmix <- function(obj,
                          newX=NULL,
                          fixomega=FALSE,
                          fixomegaval=NULL,
                          include_intercept=TRUE){

  beta <- assign_betas(obj)
  omega <- assign_omegas(obj)

  I_b0 <- as.numeric(include_intercept)

  if(is.null(newX)){
    newX <- obj$const$X
  }

  RR <- nrow(obj$Zbeta)
  if(fixomega==FALSE){ ## sample theta

    pred <- lapply(1:obj$const$K,function(kk){ ## loop over outcomes
      lapply(1:obj$const$p,function(jj){ ## loop over exposures
        sapply(1:RR,function(rr){## loop over samples
          I_b0*(obj$b0[rr,kk])+get_Btheta(newX[[jj]]%*%c(omega[[jj]][rr,(kk-1)*obj$const$L+(1:obj$const$L)]),obj$const,list(Ztheta=matrix(obj$Ztheta[rr,],nrow=obj$const$K,ncol=obj$const$p)),kk,jj)%*%(beta[[jj]][rr,(kk-1)*obj$const$d+(1:obj$const$d)])
        })
      })
    })
  }else{ ## otherwise use fixed theta

    if(!is.null(fixomegaval)){ ## use given value
      omega <- fixomegaval
    }else{## otherwise use mean values
      omega <- lapply(omega,function(th){apply(th,2,mean)})
    }
    for(jj in 1:obj$const$p){ ## standardizing if necessary
      for(kk in 1:obj$const$K){
        omega[[jj]][(kk-1)*obj$const$L+(1:obj$const$L)] <- omega[[jj]][(kk-1)*obj$const$L+(1:obj$const$L)]/sqrt(sum(omega[[jj]][(kk-1)*obj$const$L+(1:obj$const$L)]^2))
      }
    }

    pred <- lapply(1:obj$const$K,function(kk){ ## loop over outcomes
      lapply(1:obj$const$p,function(jj){ ## loop over exposures
        Btheta <- get_Btheta(newX[[jj]]%*%c(omega[[jj]][(kk-1)*obj$const$L+(1:obj$const$L)]),obj$const,list(Ztheta=matrix(obj$Ztheta[1,],nrow=obj$const$K,ncol=obj$const$p)),kk,jj)
          sapply(1:RR,function(rr){## loop over samples
            I_b0*(obj$b0[rr,kk])+Btheta%*%(beta[[jj]][rr,(kk-1)*obj$const$d+(1:obj$const$d)])
        })
      })
    })

  }


  return(list(values=pred,
              summary=summarize_pred(pred)))
}



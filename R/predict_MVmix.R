


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
                          fixtheta=FALSE,
                          fixthetaval=NULL){

  beta <- assign_betas(obj)
  theta <- assign_thetas(obj)


  if(is.null(newX)){
    newX <- obj$const$X
  }

  RR <- nrow(obj$Zbeta)
  if(fixtheta==FALSE){ ## sample theta

    pred <- lapply(1:obj$const$K,function(kk){ ## loop over outcomes
      lapply(1:obj$const$p,function(jj){ ## loop over exposures
        sapply(1:RR,function(rr){## loop over samples
          obj$b0[rr,kk]+get_Btheta(newX[[jj]]%*%c(theta[[jj]][rr,(kk-1)*obj$const$L+(1:obj$const$L)]),obj$const)%*%(beta[[jj]][rr,(kk-1)*obj$const$d+(1:obj$const$d)])
        })
      })
    })
  }else{ ## otherwise use fixed theta

    if(!is.null(fixthetaval)){ ## use given value
      theta <- fixthetaval
    }else{## otherwise use mean values
      theta <- lapply(theta,function(th){apply(th,2,mean)})
    }
    for(jj in 1:obj$const$p){ ## standardizing if necessary
      for(kk in 1:obj$const$K){
        theta[[jj]][(kk-1)*obj$const$L+(1:obj$const$L)] <- theta[[jj]][(kk-1)*obj$const$L+(1:obj$const$L)]/sqrt(sum(theta[[jj]][(kk-1)*obj$const$L+(1:obj$const$L)]^2))
      }
    }

    pred <- lapply(1:obj$const$K,function(kk){ ## loop over outcomes
      lapply(1:obj$const$p,function(jj){ ## loop over exposures
        Btheta <- get_Btheta(newX[[jj]]%*%c(theta[[jj]][(kk-1)*obj$const$L+(1:obj$const$L)]),obj$const)
          sapply(1:RR,function(rr){## loop over samples
            obj$b0[rr,kk]+Btheta%*%(beta[[jj]][rr,(kk-1)*obj$const$d+(1:obj$const$d)])
        })
      })
    })

  }


  return(list(values=pred,
              summary=summarize_pred(pred)))
}



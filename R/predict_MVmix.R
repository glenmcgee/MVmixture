


summarize_pred <- function(pred){
  summlist <- lapply(pred,function(pmat){
    summ <- data.frame(mean=apply(pmat,1,mean),
                       lower=apply(pmat,1,function(x)quantile(x,0.025)),
                       upper=apply(pmat,1,function(x)quantile(x,0.975)) )
  })
  return(summlist)
}


## function for estimating fitted/predicted values
predict_MVmix <- function(obj,
                          newX=NULL,
                          fixtheta=FALSE,
                          fixthetaval=NULL){

  if(is.null(newX)){
    newX <- obj$const$X
  }

  RR <- nrow(obj$beta)
  if(fixtheta==FALSE){ ## sample theta

    pred <- lapply(1:obj$const$K,function(kk){ ## loop over clusters
      sapply(1:RR,function(rr){## loop over samples
        get_Btheta(newX%*%c(obj$theta[rr,(kk-1)*obj$const$L+(1:obj$const$L)]),obj$const)%*%(obj$beta[rr,(kk-1)*obj$const$d+(1:obj$const$d)])
      })
    })
  }else{ ## otherwise use fixed theta

    if(!is.null(fixthetaval)){ ## use given value
      theta <- fixthetaval
    }else{## otherwise use mean values
      theta <- apply(obj$theta,2,mean)
    }
    for(kk in 1:obj$const$K){
      theta[(kk-1)*obj$const$L+(1:obj$const$L)] <- theta[(kk-1)*obj$const$L+(1:obj$const$L)]/sqrt(sum(theta[(kk-1)*obj$const$L+(1:obj$const$L)]^2))
    }

    pred <- lapply(1:obj$const$K,function(kk){ ## loop over clusters
      Btheta <- get_Btheta(newX%*%c(theta[(kk-1)*obj$const$L+(1:obj$const$L)]),obj$const)
      sapply(1:RR,function(rr){## loop over samples
        Btheta%*%(obj$beta[rr,(kk-1)*obj$const$d+(1:obj$const$d)])
      })
    })

  }


  return(list(values=pred,
              summary=summarize_pred(pred)))
}



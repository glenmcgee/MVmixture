## 4 options:
#### Vgridsearch: TRUE vs FALSE
#### DLM: TRUE vs FALSE

## define MCMC sampler options
build_sampler <- function(const){

  ## main update function
  update_params <- function(params){

    params <- update_clustMemb(params,const) ## update Zbeta and Ztheta
    params <- update_V(params,const)
    params <- update_betastar(params,const)
    params <- update_thetastar(params,const)
    params <- update_alpha(params,const)
    params <- update_logrho(params,const)
    params <- update_lambda_beta(params,const)
    params <- update_loglambda_theta(params,const) ## update DLM penalty
    params <- update_u(params,const)
    params <- update_xi(params,const)
    params <- update_sigma2_k(params,const)
    params <- update_intercept(params,const)
    params <- update_betaZk(params,const) ## linear confounder coefficients
    # print(rbind(params$Zbeta,params$Ztheta))
    # print(params$lambda_beta)
    return(params)
  }

  ## edit function with other options

  ## single outcome, no random intercepts
  if(const$K==1){
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_xi", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_u", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
    for(ll in 2:(length(body(update_params))-1)){ ## use conjugate
      if(any(grepl( "update_sigma2_k", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- substitute(params <- update_sigma2(params, const))
      }
    }
  }

  if(is.null(const$Z)){
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_betaZk", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
  }

  ## don't use DLM penalty
  if(const$DLM==FALSE | const$DLMpenalty==FALSE){
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_loglambda_theta", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
  }

  ## replace gridsearch with MH
  if(const$Vgridsearch==FALSE){
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_V", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- substitute(params <- update_V_MH(params, const))
      }
    }
  }

  ## if forcing linearity, dont update lambda_beta
  if(const$LM==TRUE){
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_lambda_beta", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
  }

  ## no theta step (or loglambda_theta step) if L==1
  if(const$L==1){
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_thetastar", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }

    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_loglambda_theta", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
  }


  ## replace rfb/approximation with MH
   if(const$approx==FALSE){
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_thetastar", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- substitute(params <- update_thetastar_MH_beta(params, const))
      }
    }
   }

  ## no joint clustering parameter unless clustering both beta and theta
  if(const$clustering!="both"){
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_logrho", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
  }

  ## use only fixed clusters
  if(const$clustering=="neither"){
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_clustMemb", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_V", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
    for(ll in 2:(length(body(update_params))-1)){
      if(any(grepl( "update_alpha", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
  }



  return(update_params)
}

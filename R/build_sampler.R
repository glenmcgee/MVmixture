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
    # params <- update_loglambda_theta(params,const) ## update DLM penalty
    params <- update_u(params,const)
    params <- update_sigma2_u(params,const)
    params <- update_sigma2(params,const)
    params <- update_intercept(params,const)

    return(params)
  }

  lenfun <- length(body(update_params))

  ## edit function with other options
  if(const$DLM==FALSE){ ## don't use DLM penalty
    for(ll in 2:(lenfun-1)){
      if(any(grepl( "update_loglambda_theta", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- NULL
      }
    }
  }

  ## replace gridsearch with MH
  if(const$Vgridsearch==FALSE){
    for(ll in 2:(lenfun-1)){
      if(any(grepl( "update_V", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- substitute(params <- update_V_MH(params, const))
      }
    }
  }

  ## replace rfb/approximation with MH (either with rfb or mvn)
  if(const$thetaMethod=="MH_vmf"){ ##
    for(ll in 2:(lenfun-1)){
      if(any(grepl( "update_thetastar", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- substitute(params <- update_thetastar_MH_vmf(params, const))
      }
    }
  }else if(const$thetaMethod=="MH_mvn"){
    for(ll in 2:(lenfun-1)){
      if(any(grepl( "update_thetastar", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- substitute(params <- update_thetastar_MH_mvn(params, const))
      }
    }
  }else if(const$thetaMethod=="MH_beta"){
    for(ll in 2:(lenfun-1)){
      if(any(grepl( "update_thetastar", as.character(body(update_params)[[ll]]), fixed = TRUE))){
        body(update_params)[[ll]] <- substitute(params <- update_thetastar_MH_beta(params, const))
      }
    }
  }



  return(update_params)
}

## define MCMC sampler options
build_sampler <- function(const){
  ## 4 possibilities:
  #### Vgridsearch: TRUE vs FALSE
  #### DLM: TRUE vs FALSE

  if(const$Vgridsearch==TRUE){ ## grid search
    if(const$DLM==TRUE){ ##
      update_params <- function(params){

        params <- update_clustMemb(params,const) ## update Zbeta and Ztheta
        params <- update_V(params,const)          ## update V_c^theta and V_c^beta
        params <- update_betastar(params,const)
        params <- update_thetastar(params,const)
        params <- update_alpha(params,const)      ## update alpha_beta and alpha_theta
        params <- update_logrho(params,const)
        params <- update_lambda_beta(params,const)
        params <- update_loglambda_theta(params,const) ## update DLM penaty
        params <- update_u(params,const)
        params <- update_sigma2_u(params,const)
        params <- update_sigma2(params,const)

        return(params)
      }
    }else{ ## No DLM penalty on thetas
      update_params <- function(params){

        params <- update_clustMemb(params,const) ## update Zbeta and Ztheta
        params <- update_V(params,const)          ## update V_c^theta and V_c^beta
        params <- update_betastar(params,const)
        params <- update_thetastar(params,const) ## update_thetastar(params,const)
        params <- update_alpha(params,const)      ## update alpha_beta and alpha_theta
        params <- update_logrho(params,const)
        params <- update_lambda_beta(params,const)
        # params <- update_loglambda_theta(params,const) ## No DLM penalty
        # params <- update_u(params,const)
        # params <- update_sigma2_u(params,const)
        params <- update_sigma2(params,const)

        return(params)
      }
    }
  }else{ ## MH step for V
    if(const$DLM==TRUE){ ##
      update_params <- function(params){

        params <- update_clustMemb(params,const)
        params <- update_V_MH(params,const)
        params <- update_betastar(params,const)
        params <- update_thetastar(params,const)
        params <- update_alpha(params,const)
        params <- update_logrho(params,const)
        params <- update_lambda_beta(params,const)
        params <- update_loglambda_theta(params,const) ## update DLM penaty
        params <- update_u(params,const)
        params <- update_sigma2_u(params,const)
        params <- update_sigma2(params,const)

        return(params)
      }
    }else{
      update_params <- function(params){

        params <- update_clustMemb(params,const)
        params <- update_V_MH(params,const)
        params <- update_betastar(params,const)
        params <- update_thetastar(params,const)
        params <- update_alpha(params,const)
        params <- update_logrho(params,const)
        params <- update_lambda_beta(params,const)
        # params <- update_loglambda_theta(params,const) ## No DLM penaty
        params <- update_u(params,const)
        params <- update_sigma2_u(params,const)
        params <- update_sigma2(params,const)

        return(params)
      }
    }
  }
}

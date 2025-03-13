#########################################
## Analyze HELIX data with MVmix
#########################################

## set up
options(echo=TRUE)
options(stringsAsFactors = FALSE)
exp_set="postnatal" ## alternative is "postnatal", "prenatal" or "all"


load("Data/exposome_v2.RData")
load("Data/covariates.rda")
load("Data/index_exposure_list.rda")
load("Data/time_exposure_list.rda")

################
## exposures ###
################
if(exp_set=="postnatal"){
  exposure_list <- index_exposure_list[grepl("Post",names(index_exposure_list))]
  exposure_list <- exposure_list[-c(2,14)]
}else if(exp_set=="prenatal"){
  exposure_list <- index_exposure_list[grepl("Pregn",names(index_exposure_list))]
}else{
  exposure_list <- index_exposure_list
}
exposure_list <- lapply(exposure_list,scale)

END <- c(cumsum(sapply(exposure_list,ncol)))
START <- END-c(4,diff(cumsum(sapply(exposure_list,ncol))))+1
groupIDs <- lapply(1:length(END),function(ll){START[ll]:END[ll]})


################
## covariates ##
################
z <- data.frame(cbind(z_base_preg,z_base_post))
# remove covariates that are already adjusted for in the outcome or related to the outcome
z <- z[,-which(colnames(z)%in%c("hs_child_age_None","e3_sex_Nonemale","hs_c_weight_None","hs_c_height_None"))]
z$e3_yearbir_None <- factor(z$e3_yearbir_None,levels=sort(unique(z$e3_yearbir_None)),labels=2003:2009)
z <- model.matrix(~.,data=z)[,-1]
N <- nrow(z)

## reduced further for preg+post analyses
z_red <- z[,-which(colnames(z)%in%c("hs_wgtgain_None","e3_gac_None","h_cohort2","h_cohort3","h_cohort4","h_cohort5","h_cohort6"))]

###############
## outcomes ###
###############
## 3 continuous outcomes measured at age 6-11

Y <- scale(cbind(phenotype$hs_zbmi_who,     ## BMI
                 phenotype$hs_correct_raven,## IQ
                 phenotype$hs_Gen_Tot))     ## neuro behaviour

#################
## train/test ###
#################
set.seed(0)
ntrain <- 800
trainids <- sort(sample(N,ntrain))
testids <- (1:N)[-trainids]

#######################################
## Fit MIM
#######################################

nit <- 30000
nburn <- 0.5*nit
nthin = 5
set.seed(0)


MIM <- MVmix(as.matrix(as.matrix(Y)[trainids,]),Reduce("cbind",exposure_list)[trainids,],Z=z_red[trainids,],
                niter=nit,nburn=nburn,nthin=nthin,
                Vgridsearch = TRUE,gridsize=10,
                MIM=TRUE,MIMorder=4,
                cluster="both",maxClusters=12,sharedlambda = FALSE,
                DLM=FALSE,approx=TRUE,
                prior_alpha_beta = c(1,5), ## encouraging clustering
                prior_alpha_theta = c(1,5))
save(MIM, file = "helix_MIM.RData")

pred_MIM <- predict_MVmix(MIM,
                          newX = Reduce("cbind",exposure_list)[testids,],
                          newZ = z_red[testids,],
                          include_intercept=TRUE,
                          allx=TRUE)

lapply(1:3,function(jj) {
  plot(pred_MIM$summary[[jj]]$mean~Y[testids,jj])
  abline(0,1,col="red")
  })
sapply(1:3,function(jj) cor(pred_MIM$summary[[jj]]$mean,Y[testids,jj]))
sapply(1:3,function(jj) lm(pred_MIM$summary[[jj]]$mean~Y[testids,jj])$coef)

boxplot(MIM$sigma2)

pc <- pairwise_clusters(MIM)
make_heatplot(pc$beta_y)
make_heatplot(pc$theta_y)

VIMs_MIM <- ExposureImportance(obj=MIM,exposures=groupIDs,nMC=50)



#######################################
## Fit MIM_nocov
#######################################

nit <- 30000
nburn <- 0.5*nit
nthin = 5
set.seed(0)


MIM_nocov <- MVmix(as.matrix(as.matrix(Y)[trainids,]),Reduce("cbind",exposure_list)[trainids,],Z=NULL,
             niter=nit,nburn=nburn,nthin=nthin,
             Vgridsearch = TRUE,gridsize=10,
             MIM=TRUE,MIMorder=4,
             cluster="both",maxClusters=12,sharedlambda = FALSE,
             DLM=FALSE,approx=TRUE,
             prior_alpha_beta = c(1,5), ## encouraging clustering
             prior_alpha_theta = c(1,5))
save(MIM_nocov, file = "helix_MIM_nocov.RData")

pred_MIM_nocov <- predict_MVmix(MIM_nocov,
                          newX = Reduce("cbind",exposure_list)[testids,],
                          # newZ = z_red[testids,],
                          include_intercept=TRUE,
                          allx=TRUE)

lapply(1:3,function(jj) {
  plot(pred_MIM_nocov$summary[[jj]]$mean~Y[testids,jj])
  abline(0,1,col="red")
})
sapply(1:3,function(jj) cor(pred_MIM_nocov$summary[[jj]]$mean,Y[testids,jj]))
sapply(1:3,function(jj) lm(pred_MIM_nocov$summary[[jj]]$mean~Y[testids,jj])$coef)

boxplot(MIM_nocov$sigma2)

pc_nocov <- pairwise_clusters(MIM_nocov)
make_heatplot(pc_nocov$beta_y)
make_heatplot(pc_nocov$theta_y)

VIMs_MIM_nocov <- ExposureImportance(obj=MIM_nocov,exposures=groupIDs,nMC=50)


#######################################
## Fit singleIndex model
## Different singleIndex per outcome
## Not clustered
#######################################

nit <- 10000
nburn <- 0.5*nit
nthin = 5
set.seed(0)


singleIndex <- MVmix(as.matrix(Y),list(sort(MIM2$omega[[1]][[2]]^2)
sort(MIM2$omega[[1]][[2]]^2)
sort(MIM2$omega[[1]][[2]]^2)
),Z=z,
             niter=nit,nburn=nburn,nthin=nthin,
             Vgridsearch = TRUE,gridsize=10,
             MIM=FALSE,
             cluster="neither",maxClusters=3,sharedlambda = FALSE,
             DLM=FALSE,approx=TRUE)
save(singleIndex, file = "helix_singleIndex.RData")

pred_singleIndex <- predict_MVmix(singleIndex,
                                  newX = list(Reduce("cbind",exposure_list)),
                                  newZ = z,
                                  include_intercept=TRUE, allx=TRUE)

lapply(1:3,function(jj) {
  plot(pred_singleIndex$summary[[jj]]$mean~Y[,jj])
  abline(0,1,col="red")
})
sapply(1:3,function(jj) cor(pred_singleIndex$summary[[jj]]$mean,Y[,jj]))
sapply(1:3,function(jj) lm(pred_singleIndex$summary[[jj]]$mean~Y[,jj])$coef)


boxplot(singleIndex$sigma2)

VIMs_singleIndex <- ExposureImportance(obj=singleIndex,exposures=groupIDs,nMC=100)




#######################################
## Fit singleIndex model
## With clustering
#######################################

set.seed(0)


singleIndexClust <- MVmix(as.matrix(as.matrix(Y)[trainids,]),list(Reduce("cbind",exposure_list)[trainids,]),Z=z[trainids,],
                     niter=nit,nburn=nburn,nthin=nthin,
                     Vgridsearch = TRUE,gridsize=10,
                     MIM=FALSE,
                     cluster="both",maxClusters=3,sharedlambda = FALSE,
                     DLM=FALSE,approx=TRUE,
                     prior_alpha_beta = c(1,5), ## encouraging clustering
                     prior_alpha_theta = c(1,5))
save(singleIndexClust, file = "helix_singleIndexClust.RData")

pred_singleIndexClust <- predict_MVmix(singleIndexClust,
                                  newX = list(Reduce("cbind",exposure_list)[testids,]),
                                  newZ = z[testids,],
                                  include_intercept=TRUE, allx=TRUE)

lapply(1:3,function(jj) {
  plot(pred_singleIndexClust$summary[[jj]]$mean~Y[testids,jj])
  abline(0,1,col="red")
})
sapply(1:3,function(jj) cor(pred_singleIndexClust$summary[[jj]]$mean,Y[testids,jj]))
sapply(1:3,function(jj) lm(pred_singleIndexClust$summary[[jj]]$mean~Y[testids,jj])$coef)


boxplot(singleIndexClust$sigma2)

pc_singleIndexClust <- pairwise_clusters(singleIndexClust)
make_heatplot(pc_singleIndexClust$beta_y)
make_heatplot(pc_singleIndexClust$theta_y)

VIMs_singleIndexClust <- ExposureImportance(obj=singleIndexClust,exposures=groupIDs[[7]],nMC=50)


#######################################
## Fit twoIndex (MIM2) model
## With clustering
#######################################

nit <- 4000 ## worked with one outcome pretty well! try with multiple
nburn <- 0.5*nit
nthin = 5
set.seed(0)


MIM2 <- MVmix(as.matrix(Y)[trainids,1:2],Reduce("cbind",exposure_list)[trainids,],Z=z[trainids,],
                          niter=nit,nburn=nburn,nthin=nthin,
                          Vgridsearch = TRUE,gridsize=10,
                          MIM=TRUE,MIMorder = 2,
                          cluster="both",maxClusters=4,sharedlambda = FALSE,
                          DLM=FALSE,approx=TRUE,
                          prior_alpha_beta = c(1,5), ## encouraging clustering
                          prior_alpha_theta = c(1,5))
save(MIM2, file = "helix_MIM2.RData")

pred_MIM2 <- predict_MVmix(MIM2,
                          newX = Reduce("cbind",exposure_list)[testids,],
                          newZ = z[testids,],
                          include_intercept=TRUE, allx=TRUE)

lapply(1:3,function(jj) {
  plot(pred_MIM2$summary[[jj]]$mean~Y[testids,jj])
  abline(0,1,col="red")
})
sapply(1:3,function(jj) cor(pred_MIM2$summary[[jj]]$mean,Y[testids,jj]))
sapply(1:3,function(jj) lm(pred_MIM2$summary[[jj]]$mean~Y[testids,jj])$coef)


boxplot(MIM2$sigma2)

pc_MIM2 <- pairwise_clusters(MIM2)
make_heatplot(pc_MIM2$beta_y)
make_heatplot(pc_MIM2$theta_y)

# top10index <- head(order(apply(singleIndex$omega[[1]][[2]]^2,2,median),decreasing=TRUE),10)
# VIMs_MIM2 <- ExposureImportance(obj=MIM2,exposures=as.list(top10index),nMC=100)

VIMs_MIM2 <- ExposureImportance(obj=MIM2,exposures=groupIDs,nMC=100,nSamp = 100)



#######################################
## Fit PRE and POST Index Model
#######################################
## one index for pregnancy, one for postnatal
Xpreg <- Reduce("cbind",lapply(time_exposure_list,function(x)x[,1]))
Xpost <- Reduce("cbind",lapply(time_exposure_list,function(x)x[,2]))

pregpost <- MVmix(as.matrix(Y)[trainids,],list(Xpreg[trainids,],Xpost[trainids,]),Z=z_red[trainids,],
                  niter=nit,nburn=nburn,nthin=nthin,
                  Vgridsearch = TRUE,gridsize=10,
                  DLM=TRUE,DLMpenalty=FALSE,lagOrder=NULL,diff=1, ## no penalty or basis expansion (since L=2)
                  cluster="both",maxClusters=6,approx=TRUE,
                  prior_alpha_beta = c(1,5), ## encouraging clustering
                  prior_alpha_theta = c(1,5))
save(pregpost, file = "helix_pregpost.RData")

pred_pregpost <- predict_MVmix(pregpost,
                           newX = list(Xpreg[testids,],Xpost[testids,]),
                           newZ = z_red[testids,],
                           include_intercept=TRUE, allx=TRUE)

lapply(1:3,function(jj) {
  plot(pred_pregpost$summary[[jj]]$mean~Y[testids,jj])
  abline(0,1,col="red")
})
sapply(1:3,function(jj) cor(pred_pregpost$summary[[jj]]$mean,Y[testids,jj]))
sapply(1:3,function(jj) lm(pred_pregpost$summary[[jj]]$mean~Y[testids,jj])$coef)

boxplot(pregpost$sigma2)

pc_pregpost <- pairwise_clusters(pregpost)
make_heatplot(pc_pregpost$beta_y)
make_heatplot(pc_pregpost$theta_y)

## important components
whichx <- lapply(1:3,function(kk){lapply(1:2,function(jj){order(apply(pregpost$theta[[jj]][[kk]]^2,2,mean))[1:5]})})
lapply(1:3,function(kk){lapply(1:2,function(jj){
  mdn <- apply(pregpost$theta[[jj]][[kk]][,whichx[[jj]][[kk]]],2,median)
  names(mdn) <- names(time_exposure_list)[whichx[[jj]][[kk]]]
  mdn
})})



#######################################
## Fit DLNM
#######################################
XDLAG_train <- lapply(time_exposure_list,function(x)x[trainids,])
DLAG <- MVmix(as.matrix(Y)[trainids,],XDLAG_train,Z=z[trainids,],
             niter=nit,nburn=nburn,nthin=nthin,
             Vgridsearch = TRUE,gridsize=10,
             DLM=TRUE,DLMpenalty=FALSE,lagOrder=NULL,diff=1, ## no penalty or basis expansion (since L=2)
             cluster="both",maxClusters=20,approx=TRUE,
             prior_alpha_beta = c(1,5), ## encouraging clustering
             prior_alpha_theta = c(1,5))
save(DLAG, file = "helix_DLAG.RData")



#
# cbmim_mod <- cbmim(y=scale(phenotype$hs_zbmi_who),x=exposure_list,z=z, ## y is outcome; ## x is list of exposure index matrices; ## z is matrix of covariates (need at least 1)
#                    niter=nit, ## number of iterations
#                    nburn=burnpct*nit,#0.5*RR, ## burn-in fraction
#                    nthin=thin, ## thinning number
#                    nchains=4, ## number of chains
#                    ncores=4,
#                    ## prior hyperparameters
#                    prior_theta_kappa=5, # previous: 1 # antoniadis: 700?
#                    prior_pi_additive=c(1,1), ## prior for pi in spike & slab on additive component
#                    prior_tau2_additive=c(10,10),#c(3,50),#c(0.001,0.001), ## shape and rate for inverse gamma on spline penalty terms
#                    prior_pi_nonadditive=c(1,1), ## prior for pi in spike & slab on non-additive component
#                    prior_rho_nonadditive=c(5,5), ## shape and rate for gamma
#                    prior_tau2_nonadditive=c(20,20),#c(3,20),## this was for the inverse gamma that matches moments with bkmr, which useswith mean 10 sd 10. #c(0.001,0.001), ## shape and rate for inverse gamma on nonadditive GP component
#                    prior_sigma2=c(1,1),#c(0.001,0.001), ## shape and rate for inverse gamma on \sigma^2
#                    ## MH tuning
#                    stepsize_theta_kappa=100, # larger means smaller steps # previous: 100. antoniadis: 1000
#                    stepsize_tau2_nonadditive=1,# previously 2 #EDIT:5, ##jumpsize/sd for gamma proposal on tau2_nonadditive
#                    stepsize_rho_nonadditive=0.2,# previously 0.5 ##jumpsize/sd for gamma proposal on rho
#                    oversample=TRUE,
#                    ## approximating large inverse via Sherman-Morrison-Woodbury
#                    invmethod=c("GPP"), ## "exact" for full matrix inversion, "rpSMW" for random projection +SMW, "lzSMWmax/min/both" is lanczos approximation with largest/smallest/both m eigencomponents +SMW, "lzDIRECTmax/min/both" is lanczos approximation to full matrix (I+PKP)
#                    # rank=200, ## rank of approximation
#                    knots=startknots,## pre-specified knots from GAM
#                    kernelfun="gaussian", ## choice of kernel function
#                    approxProj = FALSE, ## dont update P
#                    draw_h=FALSE,
#                    centering=TRUE, ## should covariates be centered to have mean=0
#                    scaling=TRUE,
#                    hierarchical=TRUE)
#
# ## save just model output in case of crash
# save(cbmim_mod, file=paste0("CBMIM/Output/cbmim_bmi_base_post",seed,".Rdata"))
#
# cbmim_pred <- pred_surface(cbmim_mod,Xnew=exposure_list,includeInt=TRUE)
# cbmim_pred_summary <- summarize_pred(cbmim_pred,assoc=FALSE,centered=FALSE)
# cbmim_pred_PSR <- summarize_pred_PSR(cbmim_pred,nchains=cbmim_mod$nchains)
#
# cbmim_chain <- list(cbmim_mod,cbmim_pred,cbmim_pred_summary,cbmim_pred_PSR)
# save(cbmim_chain,
#      file=paste0("CBMIM/Output/cbmim_bmi_base_post",".Rdata"))
#
# ## labels
# index_names <- names(exposure_list)
#
# # ## componentwise curves
# # cbmim_pred <- pred_surface(cbmim_mod,includeInt=FALSE) ## excluding intercept
# # cbmim_pred <- summarize_pred(cbmim_pred)
#
# ## indexwise curves
# cbmim_pred_ind_raw <- try(pred_surface_indexwise(cbmim_mod,includeInt=FALSE)) ## catching errors
# if(exists("cbmim_pred_ind_raw")){
#   cbmim_pred_ind <- summarize_pred(cbmim_pred_ind_raw,assoc=F)
#   cbmim_assoc_ind <- summarize_pred(cbmim_pred_ind_raw,assoc=T)
# }else{## error handling
#   cbmim_pred_ind <- NULL
#   cbmim_assoc_ind <- NULL
# }

# ## indexwise interactions
# cbmim_pred_inter_raw <- try(pred_surface_indexwise_interac(cbmim_mod,includeInt=FALSE,gridlen = 11,restrict = FALSE)) ## catching errors
# if(exists("cbmim_pred_ind_raw")){
#   cbmim_pred_inter <- summarize_pred(cbmim_pred_inter_raw,assoc=F,centered=T)
#   cbmim_assoc_inter <- summarize_pred(cbmim_pred_inter_raw,assoc=T)
# }else{## error handling
#   cbmim_pred_inter <- NULL
#   cbmim_assoc_inter <- NULL
# }


# ## indexwise PIPs (main effects+interactions)
# cbmim_PIPs <- get_PIPs(cbmim_mod)
#
# # ## save all results
# cbmim_chain <- list(mod=cbmim_mod,
#                     index_names=index_names,
#                     # pred=cbmim_pred,
#                     pred=cbmim_pred_summary,
#                     pred_ind=cbmim_pred_ind,
#                     # assoc_ind=cbmim_assoc_ind,
#                     # pred_inter=cbmim_pred_inter,
#                     # assoc_inter=cbmim_assoc_inter,
#                     PIPs=cbmim_PIPs)
# # cbmim_chain <- list(cbmim_mod)
# save(cbmim_chain,
#      file=paste0("CBMIM/Output/cbmim_bmi_base_post",".Rdata"))
#
#
#
# ########################
# ## Standard BMIM Fit  ##
# ########################
#
# nit <- 10000
# thin = 5
# burnpct <- 0.5
# set.seed(1234)
# library(bsmim2)
#
# bmim_mod <- bsmim2(y=scale(phenotype$hs_zbmi_who),
#                    x=exposure_list,
#                    z=z,
#                    niter=nit,nburn=burnpct*nit,nthin=thin,
#                    centering=T,scaling=T,
#                    prior_sigma=c(0.001,0.001),
#                    prior_lambda_shaperate=c(1,0.1),
#                    gaussian=TRUE,
#                    spike_slab=TRUE,
#                    gauss_prior=TRUE,
#                    prior_theta_slab_sd=0.25,
#                    stepsize_theta=0.1,#0.35
#                    draw_h=FALSE,
#                    num_theta_steps = 10)



test <- reshape2::melt(MIM2$omega[[1]][[2]],id.var=c("iter","x"),variable.name = "omega")
colnames(test) <- c("iter","x","omega")
ggplot(data=test,aes(x=x,y=omega,group=iter))+
  geom_line(alpha=0.005)+
  ylim(-1,1)

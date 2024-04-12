###########################################
### Prelim Sim for DLM with 2 exposures ###
###########################################

source("MVmix.R")

set.seed(1)
n <- 100
K <- 6
L <- 20
RR <- 2000
n_sims <- 500

results <- c()
for(rr in 1:n_sims){

  ## generate data
  X1 <- matrix(0,ncol=L,nrow=n)
  X1[,1] <- rnorm(n); for(ii in 2:L){X1[,ii] <- 0.5*X1[,ii-1]+rnorm(n)}
  X2 <- matrix(0,ncol=L,nrow=n)
  X2[,1] <- rnorm(n); for(ii in 2:L){X2[,ii] <- 0.5*X2[,ii-1]+rnorm(n)}
  w1 <- rep(1,L); w1 <- w1/sqrt(sum(w1^2)) ##
  w2 <- pnorm(-seq(-1,5,length=L)); w2 <- w2/sqrt(sum(w2^2)) ## DLM weights
  x11th <- X1%*%w1 ## k=1 j=1
  x12th <- X2%*%w1 ## k=1 j=2
  x21th <- X1%*%w2 ## k=2 j=1
  x22th <- X2%*%w2 ## k=2 j=2
  f11 <- sin(x11th) ## k=1 j=1
  f12 <- sin(x12th) ## k=1 j=2
  f21 <- sin(x21th) ## k=2 j=1
  f22 <- cos(x22th) ## k=2 j=2
  wmat <- rbind(matrix(w1,ncol=2*L,nrow=K/2,byrow=TRUE),
                matrix(w2,ncol=2*L,nrow=K/2,byrow=TRUE)) ## for testing only
  fmat <- matrix(c(rep(f11,K/2),rep(f21,K/2))+c(rep(f12,K/2),rep(f22,K/2)),ncol=K,nrow=n,byrow=FALSE)
  ranef <- matrix(rnorm(n,0,sqrt(0.2)),ncol=K,nrow=n,byrow=FALSE)
  eps <- matrix(rnorm(n*K,0,sqrt(0.5)),ncol=K)
  Y <- fmat+ranef+eps

  ## fit model
  MVmod <- MVmix(Y,list(X1,X2),Z=NULL,
                    niter=RR,nburn=0.5*RR,nthin=2,
                    Vgridsearch = TRUE,gridsize=10,
                    maxClusters=5,
                    DLM=TRUE,
                    approx=TRUE) ## MVN approximation


  ## estimating exposure-response surface
  pred <- predict_MVmix(MVmod,allx=TRUE)

  est_h <- sapply(1:K,function(kk){return(pred$summary[[kk]]$mean)})
  lci_h <- sapply(1:K,function(kk){return(pred$summary[[kk]]$lower)})
  uci_h <- sapply(1:K,function(kk){return(pred$summary[[kk]]$upper)})

  bias_h <- mean((est_h-fmat))
  mse_h <- mean((est_h-fmat)^2)
  width_h <- mean(uci_h-lci_h)
  cvg_h <- mean((lci_h<=fmat) & (fmat<=uci_h) )


  ## estimating weights
  MVmod$omega <- assign_omegas(MVmod,labelswitch = TRUE)

  est_w <- Reduce("cbind",lapply(MVmod$omega,function(omMat){matrix(apply(omMat,2,mean),ncol=L,nrow=K,byrow=TRUE)}))
  lci_w <- Reduce("cbind",lapply(MVmod$omega,function(omMat){matrix(apply(omMat,2,function(x)quantile(x,0.025)),ncol=L,nrow=K,byrow=TRUE)}))
  uci_w <- Reduce("cbind",lapply(MVmod$omega,function(omMat){matrix(apply(omMat,2,function(x)quantile(x,0.975)),ncol=L,nrow=K,byrow=TRUE)}))

  bias_w <- mean((est_w-wmat))
  mse_w <- mean((est_w-wmat)^2)
  width_w <- mean(uci_w-lci_w)
  cvg_w <- mean((lci_w<=wmat) & (wmat<=uci_w) )


  ## evaluating clustering
  Zmat <- cbind(MVmod$Ztheta,MVmod$Zbeta)
  agreemat <- matrix(nrow=ncol(Zmat), ncol=ncol(Zmat))
  for (ii in 1:ncol(Zmat)) {
    for (jj in 1:ncol(Zmat)) {
      agreemat[jj, ii] <- mean(Zmat[,ii] == Zmat[,jj])
    }
  }
  agreetheta <- agreemat[1:(K*2),1:(K*2)]               ## how often do theta_i and theta_j match
  agreebeta <- agreemat[K*2+(1:(K*2)),K*2+(1:(K*2))]    ## how often do beta_i and beta_j match
  agreecross <- diag(agreemat[(1:(K*2)),K*2+(1:(K*2))]) ## how often do theta_i beta_i match?


  ## save results
  res <- c(bias_h,mse_h,width_h,cvg_h,
           bias_w,mse_w,width_w,cvg_w,
           agreetheta,agreetheta,agreebeta)
  results <- rbind(results,res)

  print(c("simulation number: ",rr))
}








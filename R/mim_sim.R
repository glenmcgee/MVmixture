###########################################
### Prelim Sim for DLM with 2 exposures ###
###########################################

source("MVmix.R")

set.seed(1)
n <- 100
K <- 6
L <- 5
RR <- 2000
n_sims <- 500

results <- c()
for(rr in 1:n_sims){

  ## generate data
  X <- matrix(rnorm(L*n),ncol=L,nrow=n)
  w1 <- c(1,1,0,0,0); w1 <- w1/sqrt(sum(w1^2)) ##
  w2 <- c(0,0,1,1,0); w2 <- w2/sqrt(sum(w2^2)) ##
  x11th <- X%*%w1 ## k=1 j=1
  x12th <- X%*%w2 ## k=1 j=2
  x21th <- X%*%w1 ## k=2 j=1
  x22th <- X%*%w2 ## k=2 j=2
  f11 <- sin(x11th) ## k=1 j=1
  f12 <- (x12th)^2  ## k=1 j=2
  f21 <- sin(x21th) ## k=2 j=1
  f22 <- 0*(x22th)  ## k=2 j=2
  wmat <- cbind(matrix(w1,ncol=L,nrow=K,byrow=TRUE),
                matrix(w2,ncol=L,nrow=K,byrow=TRUE)) ## for testing only
  fmat <- matrix(c(rep(f11,K/2),rep(f21,K/2))+c(rep(f12,K/2),rep(f22,K/2)),ncol=K,nrow=n,byrow=FALSE)
  ranef <- matrix(rnorm(n,0,sqrt(0.2)),ncol=K,nrow=n,byrow=FALSE)
  eps <- matrix(rnorm(n*K,0,sqrt(0.5)),ncol=K)
  Y <- fmat+ranef+eps

  ## fit model
  MVmod <- MVmix(Y,list(X),Z=NULL,
                    niter=RR,nburn=0.5*RR,nthin=2,
                    Vgridsearch = TRUE,gridsize=10,
                    maxClusters=3,
                    DLM=FALSE,
                    approx=TRUE,
                    MIM = TRUE,
                    MIMorder=3) ## MVN approximation


  ## estimating exposure-response surface
  pred <- predict_MVmix(MVmod,allx=TRUE)

  est_h <- sapply(1:K,function(kk){return(pred$summary[[kk]]$mean)})
  lci_h <- sapply(1:K,function(kk){return(pred$summary[[kk]]$lower)})
  uci_h <- sapply(1:K,function(kk){return(pred$summary[[kk]]$upper)})

  bias_h <- mean((est_h-fmat))
  mse_h <- mean((est_h-fmat)^2)
  width_h <- mean(uci_h-lci_h)
  cvg_h <- mean((lci_h<=fmat) & (fmat<=uci_h) )


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
           agreetheta,agreetheta,agreebeta)
  results <- rbind(results,res)

  print(c("simulation number: ",rr))
}








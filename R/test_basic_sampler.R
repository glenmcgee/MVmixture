### Conclusions:
## The Taylor approximation works well
## The Fisher Bingham sampling works
## The problem: rFisherBingham fails all the time because of extreme param values
## Possible solution: approximate it via WLS distribution, then scaled.
### It's not completely legit, but if anything it should add a bit of variability
## Works when sampling beta, theta and sigma2 in a standard additive model

## testing theta sampler in simple setting
set.seed(1)
n <- 400
p <- 4
sigma2 <- 0.01
theta1 <- p:1
theta1 <- theta1/sqrt(sum(theta1^2))
X <- matrix(rnorm(n*p),ncol=p)
x1th <- X%*%theta1
f1 <- 3*sin(x1th)
dfdz <- 3*cos(x1th)
eps <- rnorm(n,0,sqrt(sigma2))
Y <- f1+eps

## getting true betas
SS <- mgcv::smoothCon(s(Xtheta,bs="bs"),data=data.frame(Xtheta=x1th),absorb.cons = FALSE)[[1]]
SSderiv <- SS ## make a smooth object for computing first order derivatives
SSderiv$deriv <- 1 ## first order derivatives
Btheta <- SS$X
g <- lm(Y~.-1,data=data.frame(Y,Btheta))
true_beta <- g$coef
plot(predict(g)~x1th)
points(f1~x1th,col="red")

## redraw outcomes with more uncertainty
set.seed(1)
sigma2=1
eps <- rnorm(n,0,sqrt(sigma2))
Y <- f1+eps


## WLS approximation
wlssteps <- 5
start <- runif(4,-1,1)
theta <- c(rFisherBingham(1,start))
for(jj in 1:wlssteps){
  Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
  DerivBtheta <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%theta))

  f <- Btheta%*%true_beta
  fprime <- c(DerivBtheta%*%true_beta)
  W <- diag(fprime^2)
  ytilde <- c(X%*%theta+ (Y-f)/fprime)

  ## WLS
  newtheta <- c(solve(t(X)%*%W%*%X,t(X)%*%W%*%ytilde))

  theta <- newtheta
}
wlsmean <- newtheta
wlsvar <- solve(t(X)%*%W%*%X)
wlsmean/sqrt(sum(wlsmean^2))

## RFB centered on previous value (1-step)
#### appears to work
mcmcsteps <- 2000
mcmcsamples <- c()
start <- runif(4,-1,1)
theta <- c(rFisherBingham(1,start))
for(rr in 1:mcmcsteps){
  Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
  DerivBtheta <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%theta))

  f <- Btheta%*%true_beta
  fprime <- c(DerivBtheta%*%true_beta)

  W <- diag(fprime^2)
  ytilde <- c(X%*%theta+ (Y-f)/fprime)

  mufb <- (1/sigma2)*t(X)%*%W%*%ytilde
  Afb <- -(0.5/sigma2)*t(X)%*%W%*%X

  theta <- c(rFisherBingham(1,mu=mufb,Aplus=Afb,mtop=100000))
  mcmcsamples <- rbind(mcmcsamples,theta)
}
apply(mcmcsamples,2,mean)
var(mcmcsamples)
plot(mcmcsamples[,1],type="l")


## RFB centered on WLS (n steps)
mcmcsteps <- 2000
mcmcsamples2 <- c()
start <- runif(4,-1,1)
theta <- c(rFisherBingham(1,start))
for(rr in 1:mcmcsteps){
  for(jj in 1:wlssteps){
    Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
    DerivBtheta <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%theta))

    f <- Btheta%*%true_beta
    fprime <- c(DerivBtheta%*%true_beta)
    W <- diag(fprime^2)
    ytilde <- c(X%*%theta+ (Y-f)/fprime)

    ## WLS
    newtheta <- c(solve(t(X)%*%W%*%X,t(X)%*%W%*%ytilde))

    theta <- newtheta
  }
  Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
  DerivBtheta <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%theta))

  f <- Btheta%*%true_beta
  fprime <- c(DerivBtheta%*%true_beta)

  W <- diag(fprime^2)
  ytilde <- c(X%*%theta+ (Y-f)/fprime)

  mufb <- (1/sigma2)*t(X)%*%W%*%ytilde
  Afb <- -(0.5/sigma2)*t(X)%*%W%*%X

  mcmcsamples2 <- rbind(mcmcsamples2,
                       c(rFisherBingham(1,mu=mufb,Aplus=Afb,mtop=100000)))
}
apply(mcmcsamples2,2,mean)
var(mcmcsamples2)
plot(mcmcsamples2[,1],type="l")

## comparing the three:
wlsmean/sqrt(sum(wlsmean^2))
apply(mcmcsamples,2,mean)
apply(mcmcsamples2,2,mean)



plot(f1~x1th)

## my estimates
Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%apply(mcmcsamples,2,mean)))
Bbeta <- Btheta%*%true_beta
points(Bbeta~x1th,col="blue")

## wls
Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%c(wlsmean/sqrt(sum(wlsmean^2)))))
Bbeta <- Btheta%*%true_beta
points(Bbeta~x1th,col="purple")

## about the same





##########################################

### sampling both beta and theta

## redraw outcomes with less uncertainty (since we dont need to avoid rfb issues)
set.seed(1)
sigma2=0.01
eps <- rnorm(n,0,sqrt(sigma2))
Y <- f1+eps


## RFB gives bugs
mcmcsteps <- 2000
thetasamples <- betasamples <- c()
start <- runif(4,-1,1)
theta <- c(rFisherBingham(1,start))
beta <- true_beta#c(rnorm(length(true_beta)))
for(rr in 1:mcmcsteps){

  ### sample theta

  Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
  DerivBtheta <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%theta))

  f <- Btheta%*%beta
  fprime <- c(DerivBtheta%*%beta)

  W <- diag(fprime^2)
  ytilde <- c(X%*%theta+ (Y-f)/fprime)

  mufb <- (1/sigma2)*t(X)%*%W%*%ytilde
  Afb <- -(0.5/sigma2)*t(X)%*%W%*%X

  theta <- c(rFisherBingham(1,mu=mufb,Aplus=Afb,mtop=100000))
  thetasamples <- rbind(thetasamples,theta)


  ### sample beta

  Vmat <- solve(0.01*diag(length(beta))+ ## from prior
                  (1/sigma2)*t(Btheta)%*%Btheta)

  beta <- c(mvtnorm::rmvnorm(n=1,
                                mean=Vmat%*%t((1/sigma2)*(t(Y)%*%Btheta)  ),
                                sigma=Vmat))
  betasamples <- rbind(betasamples,beta)

}
apply(thetasamples,2,mean)
apply(betasamples,2,mean)
plot(thetasamples[,1],type="l")
plot(betasamples[,1],type="l")





#### trying with normal approximation to rfb
mcmcsteps <- 4000
thetasamples <- betasamples <- c()
start <- runif(4,-1,1)
theta <- theta1#c(rFisherBingham(1,start))
beta <- true_beta#c(rnorm(length(true_beta)))
for(rr in 1:mcmcsteps){

  ### sample theta

  Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
  DerivBtheta <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%theta))

  f <- Btheta%*%beta
  fprime <- c(DerivBtheta%*%beta)

  W <- diag(fprime^2)
  ytilde <- c(X%*%theta+ (Y-f)/fprime)


  theta <- c(rmvnorm(1,solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%ytilde,sigma2*solve(t(X)%*%W%*%X)))
  theta <- theta/sqrt(sum(theta^2))
  thetasamples <- rbind(thetasamples,theta)


  ### sample beta

  Vmat <- solve(0.01*diag(length(beta))+ ## from prior
                  (1/sigma2)*t(Btheta)%*%Btheta)

  beta <- c(mvtnorm::rmvnorm(n=1,
                                mean=Vmat%*%t((1/sigma2)*(t(Y)%*%Btheta)  ),
                                sigma=Vmat))
  betasamples <- rbind(betasamples,beta)

}
apply(thetasamples,2,mean)
apply(betasamples,2,mean)
plot(thetasamples[,1],type="l")
plot(betasamples[,1],type="l")

plot(f1~x1th)

## my estimates
Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%apply(thetasamples,2,mean)))
Bbeta <- Btheta%*%apply(betasamples,2,mean)
points(Bbeta~x1th,col="blue")

## my betas with true thetas
Bthetawls <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta1))
Bbetawls <- Bthetawls%*%apply(betasamples,2,mean)
points(Bbetawls~x1th,col="purple")

## my thetas with true betas
Bthetafix <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%apply(thetasamples,2,mean)))
gfixtheta <- lm(Y~.-1,data=data.frame(Y,Bthetafix))
points(predict(gfixtheta)~x1th,col="red")

## basically all the same






##############################################


## sampling beta, theta and sigma2

#### trying with normal approximation to rfb
mcmcsteps <- 4000
thetasamples <- betasamples <- c()
sigma2samples <- c()
start <- runif(4,-1,1)
theta <- theta1#c(rFisherBingham(1,start))
beta <- true_beta#c(rnorm(length(true_beta)))
sig2 <- 1
for(rr in 1:mcmcsteps){

  ### sample theta

  Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
  DerivBtheta <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%theta))

  f <- Btheta%*%beta
  fprime <- c(DerivBtheta%*%beta)

  W <- diag(fprime^2)
  ytilde <- c(X%*%theta+ (Y-f)/fprime)


  theta <- c(rmvnorm(1,solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%ytilde,sig2*solve(t(X)%*%W%*%X)))
  theta <- theta/sqrt(sum(theta^2))
  thetasamples <- rbind(thetasamples,theta)


  ### sample beta

  Vmat <- solve(0.01*diag(length(beta))+ ## from prior
                  (1/sig2)*t(Btheta)%*%Btheta)

  beta <- c(mvtnorm::rmvnorm(n=1,
                                mean=Vmat%*%t((1/sig2)*(t(Y)%*%Btheta)  ),
                                sigma=Vmat))
  betasamples <- rbind(betasamples,beta)



  Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
  B_beta <- Btheta%*%beta


  ## sample sigma2

  sig2 <- 1/rgamma(1,shape=1+0.5*n,rate=1+0.5*sum((Y-B_beta)^2) )

  sigma2samples <- c(sigma2samples,sig2)


}
apply(thetasamples,2,mean)
apply(betasamples,2,mean)
mean(sigma2samples)
plot(thetasamples[,1],type="l")
plot(betasamples[,1],type="l")
plot(sigma2samples,type="l")


plot(f1~x1th)

## my estimates
Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%apply(thetasamples,2,mean)))
Bbeta <- Btheta%*%apply(betasamples,2,mean)
points(Bbeta~x1th,col="blue")

## my betas with wls thetas
Bthetawls <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%c(wlsmean/sqrt(sum(wlsmean^2)))))
Bbetawls <- Bthetawls%*%apply(betasamples,2,mean)
points(Bbetawls~x1th,col="purple")

## my thetas with proper betas
Bthetafix <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%apply(thetasamples,2,mean)))
gfixtheta <- lm(Y~.-1,data=data.frame(Y,Bthetafix))
points(predict(gfixtheta)~x1th,col="red")

## the same





#############################################################

##############################################


## sampling beta, theta and sigma2, lambda_beta

#### trying with normal approximation to rfb
set.seed(0)
mcmcsteps <- 4000
thetasamples <- betasamples <- c()
sigma2samples <- c()
lambdasamples <- c()
start <- runif(4,-1,1)
theta <- c(rFisherBingham(1,start))#theta1
beta <- c(rnorm(length(true_beta)))#true_beta
sig2 <- 1
lambda <- 1
for(rr in 1:mcmcsteps){

  ### sample theta

  Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
  DerivBtheta <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%theta))

  f <- Btheta%*%beta
  fprime <- c(DerivBtheta%*%beta)
  W <- diag(fprime^2)
  ytilde <- c(X%*%theta+ (Y-f)/fprime)

  theta <- c(mvtnorm::rmvnorm(1,solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%ytilde,sig2*solve(t(X)%*%W%*%X)))
  theta <- theta/sqrt(sum(theta^2))
  thetasamples <- rbind(thetasamples,theta)


  ### sample beta

  Vmat <- solve(lambda*SS$S[[1]]+ (1/sig2)*t(Btheta)%*%Btheta)
  beta <- c(mvtnorm::rmvnorm(n=1,
                             mean=Vmat%*%t((1/sig2)*(t(Y)%*%Btheta)  ),
                             sigma=Vmat))
  betasamples <- rbind(betasamples,beta)
  Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta))
  B_beta <- Btheta%*%beta

  ## sample lambda_beta

  lambda <- rgamma(1,shape=1+0.5*length(beta),rate=1+0.5*c(t(beta)%*%SS$S[[1]]%*%beta))
  lambdasamples <- c(lambdasamples,lambda)

  ## sample sigma2

  sig2 <- 1/rgamma(1,shape=1+0.5*n,rate=1+0.5*sum((Y-B_beta)^2) )
  sigma2samples <- c(sigma2samples,sig2)


}
apply(thetasamples,2,mean)
apply(betasamples,2,mean)
mean(sigma2samples)
mean(lambdasamples)
plot(thetasamples[,1],type="l")
plot(betasamples[,1],type="l")
plot(sigma2samples,type="l")
plot(lambdasamples,type="l")


plot(f1~x1th)

## my estimates
Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%apply(thetasamples,2,mean)))
Bbeta <- Btheta%*%apply(betasamples,2,mean)
points(Bbeta~x1th,col="blue")

## my betas with true thetas
Bthetawls <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%theta1))
Bbetawls <- Bthetawls%*%apply(betasamples,2,mean)
points(Bbetawls~x1th,col="purple")

## my thetas with proper betas
Bthetafix <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%apply(thetasamples,2,mean)))
gfixtheta <- lm(Y~.-1,data=data.frame(Y,Bthetafix))
points(predict(gfixtheta)~x1th,col="red")

## the same






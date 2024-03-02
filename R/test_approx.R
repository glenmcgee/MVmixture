## testing theta sampler in simple setting
set.seed(1)
n <- 500
p <- 4
X <- matrix(rnorm(n*p),ncol=p)
theta1 <- rep(1,p); theta1 <- theta1/sqrt(sum(theta1^2))
x1th <- X%*%theta1
f1 <- 3*sin(x1th)
dfdz <- 3*cos(x1th)
sigma2 <- 0.01
eps <- rnorm(n,0,sqrt(sigma2))
Y <- f1+eps

SS <- mgcv::smoothCon(s(Xtheta,bs="bs"),data=data.frame(Xtheta=x1th),absorb.cons = FALSE)[[1]]
SSderiv <- SS ## make a smooth object for computing first order derivatives
SSderiv$deriv <- 1 ## first order derivatives
Btheta <- SS$X

# done with n=500
g <- lm(Y~.-1,data=data.frame(Y,Btheta))
true_beta <- g$coef
plot(predict(g)~x1th)
points(f1~x1th,col="red")

Btheta <- mgcv::PredictMat(SS,data=data.frame(Xtheta=x1th))
DerivBtheta <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=x1th))

# check basis expansion
plot(f1~x1th)
points(Btheta%*%true_beta~x1th,col="red")
## deriv
plot(dfdz~x1th)
points(DerivBtheta%*%true_beta~x1th,col="red")




set.seed(1)
n <- 100
p <- 4
X <- matrix(rnorm(n*p),ncol=p)
theta1 <- rep(1,p); theta1 <- theta1/sqrt(sum(theta1^2))
x1th <- X%*%theta1
f1 <- 3*sin(x1th)
dfdz <- 3*cos(x1th)
sigma2 <- 0.01
eps <- rnorm(n,0,sqrt(sigma2))
Y <- f1+eps

thetaold <- c(rFisherBingham(1,1*rep(1,p)))
# theta <- rFisherBingham(1,10*rep(1,p))

Bthetaold <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%thetaold))
DerivBthetaold <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%thetaold))

fold <- Bthetaold%*%true_beta
fprimeold <- c(DerivBthetaold%*%true_beta)

W <- diag(fprimeold^2)
ytilde <- c(X%*%thetaold + (Y-fold)/fprimeold)


mufb <- (1/sigma2)*t(X)%*%W%*%ytilde
Afb <- -(0.5/sigma2)*t(X)%*%W%*%X

## confirmed that mu is a scaled mean direction vector
# draws <- (rFisherBingham(1000,mu=mufb,Aplus=Afb))



## the approximation is good!
nsteps <- 3
oldloss <- c()
newloss <- c()
start <- runif(4,-1,1)
for(ii in 1:1000){
  thetaold <- c(rFisherBingham(1,start))
  oldloss <- c(oldloss,sum((thetaold-theta1)^2))

  for(jj in 1:nsteps){
    Bthetaold <- mgcv::PredictMat(SS,data=data.frame(Xtheta=X%*%thetaold))
    DerivBthetaold <- mgcv::PredictMat(SSderiv,data=data.frame(Xtheta=X%*%thetaold))

    fold <- Bthetaold%*%true_beta
    fprimeold <- c(DerivBthetaold%*%true_beta)


    W <- diag(fprimeold^2)
    ytilde <- c(X%*%thetaold+ (Y-fold)/fprimeold)

    mufb <- (1/sigma2)*t(X)%*%W%*%ytilde
    Afb <- -(0.5/sigma2)*t(X)%*%W%*%X

    # newtheta <- c(mufb/sqrt(sum(mufb^2)))
    ## WLS
    newtheta <- c(solve(t(X)%*%W%*%X,t(X)%*%W%*%ytilde))
    newtheta <- c(newtheta/sqrt(sum(newtheta^2)))
    thetaold <- newtheta
  }


  newloss <- c(newloss,sum((newtheta-theta1)^2))
}
boxplot(oldloss,newloss)


# ## is the approximation good if we use true functions f and f' (obviously yes)
# oldloss <- c()
# newloss <- c()
# start <- runif(4,-1,1)
# for(ii in 1:1000){
#   thetaold <- c(rFisherBingham(1,start))
#   oldloss <- c(oldloss,sum((thetaold-theta1)^2))
#
#   for(jj in 1:nsteps){
#     x1thold <- X%*%thetaold
#     fold <- c(3*sin(x1thold))
#     fprimeold <- c(3*cos(x1thold))
#
#
#     W <- diag(fprimeold^2)
#     ytilde <- c(X%*%thetaold + (Y-fold)/fprimeold)
#
#     mufb <- (1/sigma2)*t(X)%*%W%*%ytilde
#     Afb <- -(0.5/sigma2)*t(X)%*%W%*%X
#
#     # newtheta <- c(mufb/sqrt(sum(mufb^2)))
#     ## WLS
#     newtheta <- c(solve(t(X)%*%W%*%X,t(X)%*%W%*%ytilde))
#     newtheta <- c(newtheta/sqrt(sum(newtheta^2)))
#     thetaold <- newtheta
#   }
#
#
#   newloss <- c(newloss,sum((newtheta-theta1)^2))
# }
# boxplot(oldloss,newloss)

set.seed(1)
source("MVmix.R")
n <- 400
K <- 10
p <- 20
X <- matrix(rnorm(n*p),ncol=p)
theta1 <- rep(1,p)/sqrt(p)
theta2 <- 0.2+dnorm(seq(-3,3,length=p)); theta2 <- theta2/sqrt(sum(theta2^2)) ## DLM weights
x1th <- X%*%theta1 ## TESTING
x2th <- X%*%theta2
f1 <- 3*sin(x1th)
f2 <- 5*cos(x2th)
fmat <- matrix(c(rep(f1,K/2),rep(f2,K/2)),ncol=K,nrow=n,byrow=FALSE)
ranef <- matrix(rnorm(n,0,1),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,1),ncol=K)
Y <- fmat+ranef+eps


set.seed(1)
test2 <- MVmix(Y,X,Z=NULL,niter=200,nburn=0,nthin=1,
              Vgridsearch = TRUE,gridsize=10,
              maxClusters=5,DLM=TRUE,rfbtries = 1000,lagOrder = 2)

pred <- predict_MVmix(test2)
# pred <- predict_MVmix(test2,fixtheta = TRUE,fixthetaval=c(rep(theta1,5),rep(theta2,5)))

plot(f1~x1th)
points(pred$summary[[1]]$mean~x1th,col="blue")
plot(f2~x2th)
points(pred$summary[[6]]$mean~x2th,col="blue")



set.seed(1)
source("MVmix.R")
n <- 100
K <- 8
p <- 4
X <- matrix(rnorm(n*p),ncol=p)
theta1 <- rep(1,p); theta1 <- theta1/sqrt(sum(theta1^2))
theta2 <- c(rep(c(1,0),each=p/2)); theta2 <- theta2/sqrt(sum(theta2^2)) #theta2 <- theta1 #
x1th <- X%*%theta1 ## TESTING
x2th <- X%*%theta2
f1 <- 3*sin(x1th)
f2 <- 5*cos(x2th)
fmat <- matrix(c(rep(f1,K/2),rep(f2,K/2)),ncol=K,nrow=n,byrow=FALSE)
ranef <- matrix(rnorm(n,0,sqrt(0.00000000001)),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,sqrt(0.5)),ncol=K)
Y <- fmat+ranef+eps


set.seed(1)
RR <- 2000
test2 <- MVmix(Y,X,Z=NULL,niter=RR,nburn=0.5*RR,nthin=5,
               Vgridsearch = TRUE,gridsize=10,
               maxClusters=5,DLM=FALSE,sharedlambda=TRUE)

pred <- predict_MVmix(test2)
# pred <- predict_MVmix(test2,fixtheta = TRUE,fixthetaval=c(rep(theta1,4),rep(theta2,4)))

plot(f1~x1th)
points(pred$summary[[1]]$mean~x1th,col="blue")
plot(f2~x2th)
points(pred$summary[[5]]$mean~x2th,col="blue")

for(jj in 1:K){
  boxplot(test2$beta[,(jj-1)*10+1:10])
}

for(jj in 1:K){
  boxplot(test2$theta[,(jj-1)*p+1:p],ylim=c(-1,1))
}


hist(apply(test2$Zbeta,1,function(x)length(unique(x))))
hist(apply(test2$Ztheta,1,function(x)length(unique(x))))



# # ###
# df <- data.frame(yy=c(Y),
#                  Xtheta=rep(x1th,K))
# gg <- gam(yy~s(Xtheta,bs="bs"),data=df)
#
# sss <- smoothCon(s(Xtheta,bs="bs"),data=df,absorb.cons = FALSE)
# df2 <- data.frame(yy=c(Y),sss[[1]]$X)
# truebeta <- lm(yy~.-1,data=df2)$coef
# truebetastar <- rep(truebeta,5)
# truebeta <- rep(truebeta,K)
#

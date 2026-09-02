#########################################
#### R Code for Bayes Lecture 02 ########
#########################################

#########################################
### Visualizing the Bivariate Normal
#########################################

# draw a bivariate normal contour
library(mvtnorm)

# create grid of x and y points
x.points <- seq(-3,3,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100) # save the p(x,y) to z

mu <- c(0,0) # define the mean
sigma <- matrix(c(1,0,0,1),nrow=2) # define the covariance

for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- dmvnorm(c(x.points[i],y.points[j]),mean=mu,sigma=sigma)
  }
}

# draw contours and 3D plot
contour(x.points,y.points,z)
persp(x.points,y.points,z)

#########################################
### Bayesian Linear Regression: Bodyfat Data Example
#########################################

library("MASS")
library("BAS")
data(bodyfat)
summary(bodyfat)

# Plot the relationship between a few variables:
par(mfrow=c(1,2))
plot(bodyfat$Abdomen,bodyfat$Weight,pch=19)
plot(bodyfat$Wrist,bodyfat$Weight,pch=19)

# MLE fit of regression model
model <- lm(Weight~Abdomen+Wrist,data=bodyfat)
# Check residuals
par(mfrow=c(1,1))
plot(model$fitted.values,model$residuals,pch=19)
abline(h=0,col=2)

# Set data values for y and X
y <- bodyfat$Weight
n <- length(y) # Number of data observations n
X <- cbind(rep(1,n),bodyfat$Abdomen,bodyfat$Wrist)
k <- dim(X)[2] # Number of parameters in model k (including intercept)
beta.hat <- model$coef

# Note that MLE fit here is identical to least squares estimate:
beta.hat.LS <- as.vector(solve(t(X)%*%X)%*%t(X)%*%y)
y.hat <- X%*%beta.hat.LS
sigsq.hat <- sum((y-y.hat)^2)/n

# Function that samples from the posterior of a linear model
sample.linearmodel.conj <- function(y,X,mu0,Lambda0,nu0,sigsq0){
  n <- length(y) # Number of observations
  k <- dim(X)[2] # Number of covariates
  beta.hat <- as.vector(solve(t(X)%*%X)%*%t(X)%*%y)
  
  # Step 1: Sample from marginal posterior sigma^2|y,X
  # Note that it's a little easier to sample temp from the gamma, and then take the inverse 1/temp
  nu_n <- nu0 + n
  sigsq_n <- (nu0*sigsq0 + t(y-X%*%beta.hat)%*%(y-X%*%beta.hat) + t(beta.hat - mu0)%*%Lambda0%*%(beta.hat-mu0))/(nu0+n)
  post.alpha <- nu_n/2
  post.beta <- nu_n*sigsq_n/2
  temp <- rgamma(1,shape=post.alpha,rate=post.beta)
  sigsq.samp <- 1/temp
  
  # Step 2: Sample from conditional posterior mu|sigma^2
  # Note that we are plugging in our sigsq.samp from Step 1
  post.mean <- solve(t(X)%*%X + Lambda0)%*%(t(X)%*%X%*%beta.hat + Lambda0%*%mu0)
  post.var <- sigsq.samp*solve(t(X)%*%X + Lambda0)
  beta.samp <- mvrnorm(1,mu=post.mean,Sigma=post.var)
  
  # Return vector containing posterior sample of mu and sigma^2
  out <- c(beta.samp,sigsq.samp)
  out
}

# Try a diffuse prior
mu0 <- rep(0,3)
Lambda0 <- diag(0.0001,3,3)
#Lambda0 <- diag(0.1,3,3) # Stronger prior
nu0 <- 0
sigsq0 <- 0

# Sample from posterior
numsamp <- 1000 # Set number of samples
k <- dim(X)[2] # Number of parameters
beta.samp_all <- matrix(NA,numsamp,k) # Stores samples for beta coefficients
sigsq.samp_all <- rep(NA,numsamp) # Stores samples for sigma^2
for (i in 1:numsamp){
  curr_sample <- sample.linearmodel.conj(y,X,mu0,Lambda0,nu0,sigsq0)
  beta.samp_all[i,] <- curr_sample[1:k]
  sigsq.samp_all[i] <- curr_sample[k+1]
}

# Summarizing posterior distributions: posterior means
postmean.beta <- apply(beta.samp_all,2,mean)
postmean.sigsq <- mean(sigsq.samp_all)
postmean.beta
postmean.sigsq

# 95% posterior intervals ("Credible Intervals")
quantile(beta.samp_all[,1],c(0.025,0.975))
quantile(beta.samp_all[,2],c(0.025,0.975))
quantile(beta.samp_all[,3],c(0.025,0.975))
quantile(sigsq.samp_all,c(0.025,0.975))

# posterior histograms
par(mfrow=c(2,2))
# Intercept (beta0)
hist(beta.samp_all[,1],main="Intercept (beta0)")
abline(v=postmean.beta[1],col=2,lwd=2) # Plot posterior
abline(v=beta.hat[1],col=3,lwd=2,lty=2) # Plot LS estimate
# Abdomen (beta1)
hist(beta.samp_all[,2],main="Abdomen (beta1)")
abline(v=postmean.beta[2],col=2,lwd=2) # Plot posterior
abline(v=beta.hat[2],col=3,lwd=2,lty=2) # Plot LS estimate
# Wrist (beta2)
hist(beta.samp_all[,3],main="Wrist (beta2)")
abline(v=postmean.beta[3],col=2,lwd=2) # Plot posterior
abline(v=beta.hat[3],col=3,lty=2,lwd=2) # Plot LS estimate
# sigma^2
hist(sigsq.samp_all,main="sigma^2")
abline(v=postmean.sigsq[1],col=2,lwd=2) # Plot posterior
abline(v=summary(model)$sigma^2,col=3,lwd=2,lty=2) # Plot LS estimate


# Compare the prior, likelihood, and data for Wrist (beta2)
beta1 <- seq(0,20,0.01)
# Prior
prior.mean <- mu0
prior.var <- postmean.sigsq*solve(Lambda0)
density.prior <- dnorm(beta1,mean=prior.mean[3],sd=sqrt(prior.var[3,3]))
# Data (Likelihood)
data.mean <- beta.hat
data.var <- sigsq.hat*solve(t(X)%*%X)
density.data <- dnorm(beta1,mean=data.mean[3],sd=sqrt(data.var[3,3]))
plot(beta1,density.data,ylim=c(0,2),type="l",lwd=2)
# Posterior
post.mean <- solve(t(X)%*%X + Lambda0)%*%(t(X)%*%X%*%beta.hat + Lambda0%*%mu0)
post.var <- postmean.sigsq*solve(t(X)%*%X + Lambda0)
density.post <- dnorm(beta1,mean=post.mean[3],sd=sqrt(post.var[3,3]))
# Plot together
par(mfrow=c(1,1))
plot(beta1,density.prior,type="l",col="red",lwd=2,main="Densities for beta2",xlab="beta2",ylab="Density",ylim=c(0,0.4))
lines(beta1,density.data,type="l",col="green",lwd=2,lty=2)
lines(beta1,density.post,type="l",col="blue",lwd=2,lty=3)
legend("topright",c("Prior","Data (Likelihood)","Posterior"),col=c("red","green","blue"),lwd=2,lty=c(1,2,3))



# posterior predictive distribution for a new covariate vector Xstar
Xstar <- c(1,130,17) # new person with Abdomen = 130 and Wrist = 17
ystar.samp <- rep(NA,numsamp)
for (i in 1:numsamp){
  ystar.samp[i] <- rnorm(1,mean=Xstar%*%t(t(beta.samp_all[i,])),sd=sqrt(sigsq.samp_all[i]))
}
ystar.postmean <- mean(ystar.samp)
ystar.postmean
# Compare this to the LS estimator
ystar.LS.mean <- beta.hat%*%Xstar
# Plot the distribution
par(mfrow=c(2,1))
xmin <- min(y,ystar.postmean)
xmax <- max(y,ystar.postmean)
hist(y,main="Dataset Weights",xlim=c(xmin,xmax))
hist(ystar.samp,main="Predicted Weight of Person with abdomen = 130 and wrist = 17",xlim=c(xmin,xmax))
abline(v=ystar.postmean,col=2,lwd=2) # Plot posterior prediction
abline(v=ystar.LS.mean,col=3,lwd=2,lty=2) # Plot LS estimator



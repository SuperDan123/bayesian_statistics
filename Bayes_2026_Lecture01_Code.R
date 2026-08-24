#########################################
#### R Code for Bayes Lecture 01 ########
#########################################

#########################################
### Unif (0,1) is equivalent to Beta(1,1)
#########################################

theta <- ppoints(1000)
density.beta <- dbeta(theta,1,1)
plot(theta,density.beta,type="l",col="red",lwd=2,main="Beta(1,1)",ylim=c(0,2),xlab="theta",ylab="Density")

#########################################
### Take a look at the beta distribution
#########################################

# Original course code used manipulate (RStudio-only).
# Use shiny instead so sliders work in Cursor / VS Code / browser.
# install.packages("shiny") if needed; then run this whole block.
library(shiny)

theta <- ppoints(1000)
shinyApp(
  ui = fluidPage(
    titlePanel("Beta(alpha, beta) density"),
    sliderInput("alpha", "alpha", min = 0.01, max = 25, value = 1, step = 0.01),
    sliderInput("beta", "beta", min = 0.01, max = 25, value = 1, step = 0.01),
    plotOutput("beta_plot", height = "400px")
  ),
  server = function(input, output) {
    output$beta_plot <- renderPlot({
      y <- dbeta(theta, input$alpha, input$beta)
      plot(theta, y, type = "l", lwd = 2, col = "red",
           main = sprintf("Beta(%.2f, %.2f)", input$alpha, input$beta),
           xlab = "theta", ylab = "Density", ylim = c(0, max(2, max(y) * 1.05)))
    })
  }
)

#########################################
### Beta-Binomial Example: Coin Flips
#########################################

### Example (1) ###########

# Suppose our prior is alpha = 1, beta = 1 (uniform prior)
alpha <- 1
beta <- 1

# Suppose we flip 10 coins, and observe 8 heads and 2 tails
numheads <- 8
numtails <- 2

# Calculate prior, data (likelihood), and posterior distributions
theta <- ppoints(1000)
density.prior <- dbeta(theta,alpha,beta)
density.data <- dbinom(numheads,numheads+numtails,theta)
density.posterior <- dbeta(theta,alpha+numheads,beta+numtails)

# Plot distributions
plot(theta,density.prior,type="l",col="red",lwd=2,main="Prior p(theta)",ylim=c(0,2),xlab="theta",ylab="Density")
plot(theta,density.data,type="l",col="blue",lwd=2,main="Data Likelihood p(y|theta)",ylim=c(0,max(density.data)),xlab="theta",ylab="Density")
plot(theta,density.posterior,type="l",col="green",lwd=2,main="Posterior (theta|y)",ylim=c(0,max(density.posterior)),xlab="theta",ylab="Density")


### Example (2) ###########

# Suppose our prior is alpha = 5, beta = 5 (mass at 0.5)
alpha <- 5
beta <- 5

# Suppose we flip 10 coins, and observe 8 heads and 2 tails
numheads <- 8
numtails <- 2

# Calculate prior, data (likelihood), and posterior distributions
theta <- ppoints(1000)
density.prior <- dbeta(theta,alpha,beta)
density.data <- dbinom(numheads,numheads+numtails,theta)
density.posterior <- dbeta(theta,alpha+numheads,beta+numtails)

# Plot distributions
plot(theta,density.prior,type="l",col="red",lwd=2,main="Prior p(theta)",ylim=c(0,3),xlab="theta",ylab="Density")
plot(theta,density.data,type="l",col="blue",lwd=2,main="Data Likelihood p(y|theta)",ylim=c(0,max(density.data)),xlab="theta",ylab="Density")
plot(theta,density.posterior,type="l",col="green",lwd=2,main="Posterior (theta|y)",ylim=c(0,max(density.posterior)),xlab="theta",ylab="Density")


#########################################
### Plot the 95% credible intervals for Beta(2,2)
#########################################

theta <- ppoints(1000)
density.beta <- dbeta(theta,2,2)
lower.bound <- qbeta(0.025,2,2)
upper.bound <- qbeta(0.975,2,2)
plot(theta,density.beta,type="l",col="red",lwd=2,main="Beta(1,1)",ylim=c(0,2),xlab="theta",ylab="Density")
abline(v=lower.bound,lty=2,col="black",lwd=2)
abline(v=upper.bound,lty=2,col="black",lwd=2)

#########################################
### Beta-Binomial Example: Different Priors
#########################################

theta <- ppoints(1000)
density.unif <- dbeta(theta,1,1)
density.jeff <- dbeta(theta,1/2,1/2)
density.neut <- dbeta(theta,1/3,1/3)

density.impr <- dbeta(theta,0,0)
density.impr <- dbeta(theta,0.00001,0.00001)

mindensity <- min(density.unif,density.jeff,density.neut,density.impr)
maxdensity <- max(density.unif,density.jeff,density.neut,density.impr)
plot(theta,density.unif,type="l",col="black",lwd=2,main="Prior Distributions of Theta",ylim=c(mindensity,maxdensity),xlab="theta",ylab="density")
lines(theta,density.jeff,type="l",col="blue",lwd=2)
lines(theta,density.neut,type="l",col="red",lwd=2)
lines(theta,density.impr,type="l",col="green",lwd=2)
legend(0.4,25,c("Beta(1,1)","Beta(1/2,1/2)","Beta(1/3,1/3)","Beta(0,0)"),col=c("black","blue","red","green"),lwd=2)

### Comparing Posterior for Different Priors for Small Dataset ###
### Small Dataset: 2 successes, 1 failure

numsucc <- 2
numfail <- 1

theta.mle <- numsucc/(numsucc+numfail)

theta <- ppoints(1000)
posterior.unif <- dbeta(theta,numsucc+1,numfail+1)
posterior.jeff <- dbeta(theta,numsucc+1/2,numfail+1/2)
posterior.neut <- dbeta(theta,numsucc+1/3,numfail+1/3)
posterior.impr <- dbeta(theta,numsucc+0.000001,numfail+0.000001)
minpost <- min(posterior.unif,posterior.jeff,posterior.neut,posterior.impr)
maxpost <- max(posterior.unif,posterior.jeff,posterior.neut,posterior.impr)
plot(theta,posterior.unif,type="l",col="black",lwd=2,main="Posterior Distributions for Theta (2 Successes, 1 Failure)",ylim=c(minpost,maxpost),xlab="theta",ylab="density")
lines(theta,posterior.jeff,type="l",col="blue",lwd=2)
lines(theta,posterior.neut,type="l",col="red",lwd=2)
lines(theta,posterior.impr,type="l",col="green",lwd=2)
legend(0.05,2,c("Beta(1,1)","Beta(1/2,1/2)","Beta(1/3,1/3)","Beta(0,0)"),col=c("black","blue","red","green"),lwd=2)
abline(v=theta.mle)


### Comparing Posterior for Different Priors for Medium Dataset ###
### Medium Dataset: 20 successes, 10 failure

numsucc <- 20
numfail <- 10

theta.mle <- numsucc/(numsucc+numfail)

theta <- ppoints(1000)
posterior.unif <- dbeta(theta,numsucc+1,numfail+1)
posterior.jeff <- dbeta(theta,numsucc+1/2,numfail+1/2)
posterior.neut <- dbeta(theta,numsucc+1/3,numfail+1/3)
posterior.impr <- dbeta(theta,numsucc+0.000001,numfail+0.000001)
minpost <- min(posterior.unif,posterior.jeff,posterior.neut,posterior.impr)
maxpost <- max(posterior.unif,posterior.jeff,posterior.neut,posterior.impr)
plot(theta,posterior.unif,type="l",col="black",lwd=2,main="Posterior Distributions for Theta (20 Successes, 10 Failures)",ylim=c(minpost,maxpost),xlab="theta",ylab="density")
lines(theta,posterior.jeff,type="l",col="blue",lwd=2)
lines(theta,posterior.neut,type="l",col="red",lwd=2)
lines(theta,posterior.impr,type="l",col="green",lwd=2)
legend(0.05,2,c("Beta(1,1)","Beta(1/2,1/2)","Beta(1/3,1/3)","Beta(0,0)"),col=c("black","blue","red","green"),lwd=2)
abline(v=theta.mle)




#########################################
### Mixture of Betas: Credible Interval vs. Highest Density Region
#########################################

##### prepare the values to evaluate #####

# setup value quantiles of the distribution
p=seq(0,1,.002)               # notice that p is a vector
m=length(p)                   # length of m

# compute the density as a mixture of densities (or replace with your posterior density)
dens=0.5*dbeta(p,4,12)+0.5*dbeta(p,16,4)  # compute the density
ndens=(dens/sum(dens))         # normalize the density to it sums to 1



##### plot the density and compute probability of interest #####

# select the region of the area of interest
select.p = (p>0.2 & p<0.4) | (p>0.6 & p<0.8)   # create an indicator of locations to highlight
select.dens = ifelse(select.p,dens,0)    # create a vector of the density that corresponds with selected points

# what is the density selected
sum(ndens[select.p])

# plot the density
plot(p,dens,type='l')
polygon( c(p,rev(p)), c(select.dens,rep(0,m)), col='gray')  # shade area under curve



##### compute the credible region #####

# set the confidence level (say 0.95)
confidence=0.95  # confidence level (say 0.95)
alpha=1-confidence    # complement of region, so 95% would give alpha=0.05

# find the 95% credible region
# by locating the lower and upper quartile that correspond with the region
cdens = cumsum(ndens)   # compute cumulative probability
cregion.lowerindex = which(cdens >= alpha/2)[1]   # index of lower boundary
cregion.lower = p[cregion.lowerindex]  # quantile of the lower boundary
cregion.upperindex = which(cdens >= (1-alpha/2))[1]  # index of upper boundary
cregion.upper = p[cregion.upperindex]  # quantile of the upper boundary

# summarize CR and select the area of interest
paste(round(100*confidence),"% Confidence region: (",cregion.lower,",",cregion.upper,")")
select.pcr = (p>cregion.lower & p<cregion.upper)
prob.pcr=sum(ndens[select.pcr])  # compute the probability of this area (should =confidence)
paste("Probability of area: ",round(100*prob.pcr),"%")
select.denscr = ifelse(select.pcr,dens,0)  # vector of density that corresponds with CR



##### compute the highest credible region #####

# find the highest credible region
sort.dens = order(ndens, decreasing = TRUE)  # save the sorted index order
select.phcr.index = sort.dens[which(cumsum(ndens[sort.dens]) <= (1-alpha))]  # find indices of critcal values
select.phcr = rep(FALSE,m)  # create a vector to identify density to highlight
select.phcr[select.phcr.index] = TRUE

# summarize HCR and select the area of interest
prob.phcr=sum(ndens[select.phcr])
paste("Probability of area: ",round(100*prob.phcr),"%")
segment.phcr = which(diff(select.phcr)!=0)  # find beginning and end of segments
nsegments=(length(segment.phcr)/2)  # number of segments
string.segments=rep("",nsegments)  # create storage for the segments
for (i in 1:nsegments ) {
  ilower=(i-1)*2+1
  ihigh=(i-1)*2+2
  string.segments[i]=paste("Segment: ",i," Region: (",p[segment.phcr[ilower]],", ",p[segment.phcr[ihigh]],")")
}
paste(round(100*confidence),"% Highest Confidence region: ",string.segments)
select.denshcr = ifelse(select.phcr,dens,0)  # vector of density that corresponds with CR



##### plot the credible regions #####

# setup plots
par(mfrow=c(1,3))  # side by side plots
#par(mfrow=c(1,1))  # each plot in a separate panel

# plot the CR and HCR side by side
plot(p,dens,type='l',xlab="Quantile",ylab="Density",main="Credible Interval")  # plot density
polygon( c(p,rev(p)), c(select.denscr,rep(0,m)), col='gray') # color credible region
plot(p,dens,type='l',xlab="Quantile",ylab="Density",main="Highest Density Region")  # plot density
polygon( c(p,rev(p)), c(select.denshcr,rep(0,m)), col='lightblue')  # color highest credible region

# overlay the CR and HCRs
plot(p,dens,type='l',xlab="Quantile",ylab="Density",main="Credible Intervals")  # plot density
polygon( c(p,rev(p)), c(select.denscr,rep(0,m)), col='gray') # color credible region
polygon( c(p,rev(p)), c(select.denshcr,rep(0,m)), col='lightblue', density=50, angle=45)  # color highest credible region
legend("topleft",c("Credible Interval","Highest Density Region"),pch=15,col=c("gray","lightblue"),bty='n')




#########################################
### Full Normal: Baseball "Bayesball" Example
#########################################

# Reading in Data:
data <- read.csv("hitters.post1975.csv")
dim(data)
head(data)

# Reducing data to player-seasons where AB >= 200 ("at bat")
data <- data[data$AB >= 200,]
dim(data)

# Calculating batting average
Batting.Average <- data$H/data$AB # Hits divided by At Bat
n <- length(Batting.Average)
hist(Batting.Average,main="Histogram of Batting Average")

# Function that samples from the Posterior of a Normal Distribution
sample.norm.conj <- function(y,mu0,kappa0,nu0,sigsq0,numsamp){
  n <- length(y) # Number of observations
  
  # Step 1: Sample from marginal posterior sigma^2|y
  # Note that it's a little easier to sample x from the gamma, and then take the inverse 1/x
  post.alpha <- (nu0+n)/2
  post.beta <- (nu0*sigsq0+var(y)*(n-1)+(n*kappa0/(n+kappa0))*(mean(y)-mu0)^2)/2
  x <- rgamma(numsamp,shape=post.alpha,rate=post.beta)
  sigsq.samp <- 1/x
  
  # Step 2: Sample from conditional posterior mu|sigma^2
  # Note that we are plugging in our sigsq.samp from Step 1
  post.mean <- (n*mean(y)/sigsq.samp + kappa0*mu0/sigsq.samp)/(n/sigsq.samp + kappa0/sigsq.samp)
  post.var <- 1/(n/sigsq.samp + kappa0/sigsq.samp)
  mu.samp <- rnorm(numsamp,mean=post.mean,sd=sqrt(post.var))
  
  # Return vector containing posterior sample of mu and sigma^2
  out <- cbind(mu.samp,sigsq.samp)
  out
}

# checking posteriors for different conjugate priors and non-informative priors
# (different values of mu0,kappa0)

theta1 <- sample.norm.conj(Batting.Average,0.24,100,10,10,1000) # mu0 = 0.24, kappa0 = 100, etc.
theta2 <- sample.norm.conj(Batting.Average,0.24,1000,10,10,1000) # mu0 = 0.24, kappa0 = 1000, etc.
theta3 <- sample.norm.conj(Batting.Average,0.24,10000,10,10,1000) # mu0 = 0.24, kappa0 = 10000, etc.
theta4 <- sample.norm.conj(Batting.Average,0.24,0,0,10,1000) # mu0 = 0.24, kappa0 = 0, etc. (non-informative)

par(mfrow=c(2,2))
minmu <- min(theta1[,1],theta2[,1],theta3[,1],theta4[,1],mean(Batting.Average))
maxmu <- max(theta1[,1],theta2[,1],theta3[,1],theta4[,1],mean(Batting.Average))
hist(theta1[,1],main="Mu: mu0=0.2, kappa0=100",xlim=c(0.24,0.28),xlab="mu")
abline(v=mean(Batting.Average),col=2,lwd=2) # Data mean
abline(v=0.24,col=3,lwd=2) # Prior mean
hist(theta2[,1],main="Mu: mu0=0.2, kappa0=1000",xlim=c(0.24,0.28),xlab="mu")
abline(v=mean(Batting.Average),col=2,lwd=2) # Data mean
abline(v=0.24,col=3,lwd=2) # Prior mean
hist(theta3[,1],main="Mu: mu0=0.2, kappa0=10000",xlim=c(0.24,0.28),xlab="mu")
abline(v=mean(Batting.Average),col=2,lwd=2) # Data mean
abline(v=0.24,col=3,lwd=2) # Prior mean
hist(theta4[,1],main="Mu: non-informative",xlim=c(0.24,0.28),xlab="mu")
abline(v=mean(Batting.Average),col=2,lwd=2) # Data mean
abline(v=0.24,col=3,lwd=2) # Prior mean

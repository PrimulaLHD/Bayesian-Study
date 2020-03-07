#Bayesian statistics uisng R and jags

# From MOOC: https://www.coursera.org/learn/mcmc-bayesian-statistics
# as well as in http://www4.stat.ncsu.edu/~reich/ST590/

library(rjags)

dat   <- read.csv("http://www4.stat.ncsu.edu/~reich/ST590/assignments/Obama2012.csv")
Y     <- 100*dat[,2]
Y     <- (Y-mean(Y))/sd(Y)
white <- dat[,7]
white <- (white-mean(white))/sd(white)
unemp <- dat[,18]
unemp <- (unemp-mean(unemp))/sd(unemp)
n     <- 100

model_string <- "model{
  # Likelihood
  for(i in 1:n){
    Y[i]   ~ dnorm(mu[i],inv.var)
    mu[i] <- beta[1] + beta[2]*white[i] + beta[3]*unemp[i]
  }

  # Prior for beta
  for(j in 1:3){
    beta[j] ~ dnorm(0,0.0001)
  }

  # Prior for the inverse variance
  inv.var   ~ dgamma(0.01, 0.01)
  sigma     <- 1/sqrt(inv.var)

}"

model <- jags.model(textConnection(model_string), 
                    data = list(Y=Y,n=n,white=white,unemp=unemp),n.chains = 3)

update(model, 10000, progress.bar="none"); # Burnin for 10000 samples

samp <- coda.samples(model, 
                     variable.names=c("beta","sigma"), 
                     n.iter=20000, progress.bar="none")

summary(samp)

plot(samp)

traceplot(samp)



####################################################################
#
# Example of joint, marginal, and conditional distributions (sample)
# Using the Obama voting data.
#
####################################################################


#Load data

dat      <- read.csv("http://www4.stat.ncsu.edu/~reich/ST590/assignments/Obama2012.csv")
pctObama <- 100*dat[,2]
pctUnEmp <- dat[,18]



############################################
# Convert to discrete variables
############################################

X        <- ifelse(pctObama>50,1,0)
Y        <- ifelse(pctUnEmp>10,1,0)+
    ifelse(pctUnEmp>15,1,0)
plot(Y,pctUnEmp)

# Compute the sample joint distribution

table(X,Y)/100

# Compute the sample marginal distributions

table(X)/100
table(Y)/100

# Compute the conditional probabilities

mean(X[Y==0])
mean(X[Y==1])
mean(X[Y==2])



############################################
# Plot for continuous variables
############################################

X <- pctObama
Y <- pctUnEmp

#Joint

plot(X,Y,main="Joint distribution")

# Marginals

hist(X,main="Marginal distribution of X")
hist(Y,main="Marginal distribution of Y")

# Probability in a set

inA <- (X>50) & (Y>10) & (Y<15)  
plot(X,Y,col=ifelse(inA,2,1),main="Prob in set A")
polygon(c(50,50,100,100,50),c(10,15,15,10,10))

# Approximate conditional pdf

Y10 <- Y>9.5 & Y<10.5
plot(X,Y,col=ifelse(Y10,2,1))
abline(9.5,0)
abline(10.5,0)

X10 <- X[Y10]
hist(X10,main="f(x|Y=10)",xlim=range(X),prob=TRUE)

binorm<-function(x,y,muX=0,muY=0,sigmaX=1,sigmaY=1,rho=0){
    
    c    <- sigmaX*sigmaY*sqrt(1-rho^2)*2*pi
    d    <- 1/(1-rho^2)
    z_x  <- (x-muX)/sigmaX
    z_y  <- (y-muY)/sigmaY
    
    pdf <- (1/c)*exp(-0.5*d*(z_x^2+z_y^2-2*rho*z_x*z_y))
    
    return(pdf)
    }

m    <- 100
pts  <- seq(-3,3,length=m)
grid <- expand.grid(pts,pts)
plot(grid)


muX    <- 0
muY    <- 0
sigmaX <- 1
sigmaY <- 1
rho    <- 0.9

pdf    <- binorm(grid[,1],grid[,2],muX,muY,sigmaX,sigmaY,rho)
pdf    <- matrix(pdf,m,m)

library(fields)
image.plot(pts,pts,pdf,
           xlab="x",ylab="y",
           main="Bivariate normal PDF",
           col=gray(1-seq(0,1,.05)))


######################
post_prob<-function(p,q0,q1){
    p*q1/(p*q1+(1-p)*q0)  
}

p  <- 0.50   # Prior probability
q0 <- 0.01   # False positive probability
q1 <- 0.90   # True positive probability

post_prob(p,q0,q1)

grid  <- seq(0.01,0.99,.01)

plot(grid,post_prob(grid,q0,q1),
     type="l",
     xlab="Prior probability",
     ylab="Posterior probability")


plot(grid,post_prob(p,grid,q1),
     type="l",
     xlab="False positive rate",
     ylab="Posterior probability") 

plot(grid,post_prob(p,q0,grid),
     type="l",
     xlab="True positive rate",
     ylab="Posterior probability") 

##MCMC
n     <- 10000
theta <- NULL
Y     <- NULL

#start sampling
for(i in 1:n){
    theta[i] <- rbinom(1,1,p)
    prob     <- ifelse(theta[i]==1,q1,q0)
    Y[i]     <- rbinom(1,1,prob)
}

table(Y,theta)/n  # Approximate joint distribution

mean(theta[Y==1]) # Approximate conditional probability

################################################
#
# Complete analysis of a beta-binomial model
#
################################################


# Data:

n <- 27  # Number of trials
Y <- 15  # Number of successes


# Set working directory where all output will be stored:




# Define three priors:

a   <- c(0.1, 1, 10)
b   <- c(0.1, 1, 10)
lab <- c("Bathtub, a=b=0.1", "Uniform, a=b=1", "Informative, a=b=10")


# Plot the posterior on a grid for the three priors:

theta <- seq(0,1,length=100)
post1 <- dbeta(theta,Y+a[1],n-Y+b[1])
post2 <- dbeta(theta,Y+a[2],n-Y+b[2])
post3 <- dbeta(theta,Y+a[3],n-Y+b[3])

cex <- 1.25



plot(theta,post3,type="l",lty=3,
     cex.lab=cex,cex.axis=cex,
     xlab=expression(theta),ylab="Posterior")
lines(theta,post1,lty=1)
lines(theta,post2,lty=2)

legend("topright",lab,inset=0.05,lty=1:3,cex=cex)



# Summarize the posterior in a table:

A    <- Y+a
B    <- n-Y+b 
Mean <- A/(A+B)
Var  <- A*B/((A+B)*(A+B)*(A+B+1))
SD   <- sqrt(Var)
Q05  <- qbeta(0.05,A,B)
Q95  <- qbeta(0.95,A,B)
P50  <- pbeta(0.50,A,B)

output <- cbind(Mean,SD,Q05,Q95,P50)
output <- round(output,2)

rownames(output)<-lab

Fagaceae <- name_lookup(query='Fagaceae', rank="family", return="data")

keys <- name_suggest(q='Fagaceae', rank="family")$key[1]
occ <- occ_search(taxonKey=keys, return='data')

str(occ)

x <- map_fetch(taxonKey = keys)
library(raster)
plot(x)


set.seed(32)
m = 100
a = 2.0
b = 1.0/3.0
theta=rgamma(n=m,shape = a,rate = b)
hist(theta,freq=FALSE)
curve(dgamma(x,shape = a,rate = b),col="blue",add = T)

x=c(.612,.388)
pi = matrix(c(.365,.635,1,0),nrow  = 2,ncol = 2,byrow = T)

x%*%pi

#########Day 2
library(rjags)
#1. specify the model
model_string = "model{
    for (i in 1:n) {
       y[i]~dnorm(mu,1/sig2) 
    }
    mu ~ dt(0,1/1,1)
    sig2 = 1
}"

#2. Ste up the model
set.seed(50)
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
n = length(y)

data_jags = list(y=y,n=n)
params = c("mu")

inits = function(){
     inits = list("mu"=0)
 }

mod = jags.model(textConnection(model_string),
                 data=data_jags,inits = inits)


# 3. Run the MMCMC samplers
update(mod,500)
mod_sim = coda.samples(model = mod,variable.names = params,
                       n.iter = 1000)

# 4. Post processing
library("coda")
plot(mod_sim)
summary(mod_sim)

###########################################################
#Gibbs sampler in R

update_mu = function(n, ybar, sig2, mu_0, sig2_0) {
    sig2_1 = 1.0 / (n / sig2 + 1.0 / sig2_0)
    mu_1 = sig2_1 * (n * ybar / sig2 + mu_0 / sig2_0)
    rnorm(n=1, mean=mu_1, sd=sqrt(sig2_1))
}

gibbs = function(y, n_iter, init, prior) {
    ybar = mean(y)
    n = length(y)
    
    ## initialize
    mu_out = numeric(n_iter)
    sig2_out = numeric(n_iter)
    
    mu_now = init$mu
    
    ## Gibbs sampler
    for (i in 1:n_iter) {
        sig2_now = update_sig2(n=n, y=y, mu=mu_now, nu_0=prior$nu_0, beta_0=prior$beta_0)
        mu_now = update_mu(n=n, ybar=ybar, sig2=sig2_now, mu_0=prior$mu_0, sig2_0=prior$sig2_0)
        
        sig2_out[i] = sig2_now
        mu_out[i] = mu_now
    }
    
    cbind(mu=mu_out, sig2=sig2_out)
}

y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
ybar = mean(y)
n = length(y)

## prior
prior = list()
prior$mu_0 = 0.0
prior$sig2_0 = 1.0
prior$n_0 = 2.0 # prior effective sample size for sig2
prior$s2_0 = 1.0 # prior point estimate for sig2
prior$nu_0 = prior$n_0 / 2.0 # prior parameter for inverse-gamma
prior$beta_0 = prior$n_0 * prior$s2_0 / 2.0 # prior parameter for inverse-gamma

hist(y, freq=FALSE, xlim=c(-1.0, 3.0)) # histogram of the data
curve(dnorm(x=x, mean=prior$mu_0, sd=sqrt(prior$sig2_0)), lty=2, add=TRUE) # prior for mu
points(y, rep(0,n), pch=1) # individual data points
points(ybar, 0, pch=19) # sample mean

set.seed(53)

init = list()
init$mu = 0.0

post = gibbs(y=y, n_iter=1e3, init=init, prior=prior)

head(post)

library("coda")
plot(as.mcmc(post))
summary(as.mcmc(post))

library(bipartite)
testdata <- data.frame(higher = c("bee1","bee1","bee1","bee2","bee1","bee3"),
                       lower = c("plant1","plant2","plant1","plant2","plant3","plant4"),
                       webID = c("meadow","meadow","meadow","meadow","bog","bog"), freq=c(5,1,1,1,3,7))
networkss <- frame2webs(testdata,type.out="list",emptylist = F)

degreedistr(Safariland)
betanew.dist(networkss)



#######Day3
library("car")
data("Leinhardt")
?Leinhardt
head(Leinhardt)
str(Leinhardt)
pairs(Leinhardt)
plot(infant ~ income, data=Leinhardt)
hist(Leinhardt$infant)

Leinhardt$loginfant = log(Leinhardt$infant)
Leinhardt$logincome = log(Leinhardt$income)

plot(loginfant ~ logincome, data=Leinhardt)

#Modling
lmod = lm(loginfant ~ logincome, data=Leinhardt)
summary(lmod)
# Modeling in Jags
dat = na.omit(Leinhardt)

library("rjags")

mod1_string = " model {
    for (i in 1:n) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = b[1] + b[2]*log_income[i] 
}

for (i in 1:2) {
b[i] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*10.0/2.0)
sig2 = 1.0 / prec
sig = sqrt(sig2)
} "

set.seed(72)
data1_jags = list(y=dat$loginfant, n=nrow(dat), 
                  log_income=dat$logincome)

params1 = c("b", "sig")

inits1 = function() {
    inits = list("b"=rnorm(2,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod1 = jags.model(textConnection(mod1_string), data=data1_jags, inits=inits1, n.chains=3)
update(mod1, 1000) # burn-in

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=50000)

mod1_csim = do.call(rbind, mod1_sim) # combine multiple chains



###MCMC CHECK/ MCMC convergence

plot(mod1_sim)

gelman.diag(mod1_sim)

autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

summary(mod1_sim)


#Residual checks

lmod0 = lm(infant ~ income, data=Leinhardt)
plot(resid(lmod0)) # to check independence (looks okay)

plot(predict(lmod0), resid(lmod0)) # to check for linearity, constant variance (looks bad)

qqnorm(resid(lmod0)) # to check Normality assumption (we want this to be a straight line)
X = cbind(rep(1.0, data1_jags$n), data1_jags$log_income)
head(X)
(pm_params1 = colMeans(mod1_csim)) # posterior mean

yhat1 = drop(X %*% pm_params1[1:2])
resid1 = data1_jags$y - yhat1
plot(resid1) # against data index
plot(yhat1, resid1) # against predicted values
qqnorm(resid1) # checking normality of residuals
plot(predict(lmod), resid(lmod)) # to compare with reference linear model
rownames(dat)[order(resid1, decreasing=TRUE)[1:5]] # which countries have the largest positive residuals?


#### Linear regression, Additional covariates
library("rjags")

mod2_string = " model {
for (i in 1:length(y)) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = b[1] + b[2]*log_income[i] + b[3]*is_oil[i]
}

for (i in 1:3) {
b[i] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*10.0/2.0)
sig = sqrt( 1.0 / prec )
} "


set.seed(73)
data2_jags = list(y=dat$loginfant, log_income=dat$logincome,
                  is_oil=as.numeric(dat$oil=="yes"))
data2_jags$is_oil

params2 = c("b", "sig")

inits2 = function() {
    inits = list("b"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, inits=inits2, n.chains=3)
update(mod2, 1e3) # burn-in

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)

mod2_csim = as.mcmc(do.call(rbind, mod2_sim)) # combine multiple chains
plot(mod2_sim)
gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)
summary(mod2_sim)

(pm_params2 = colMeans(mod2_csim)) # posterior mean
yhat2 = drop(X2 %*% pm_params2[1:3])
resid2 = data2_jags$y - yhat2
plot(resid2) # against data index

plot(yhat2, resid2) # against predicted values
plot(yhat1, resid1) # residuals from the first model
sd(resid2) # standard deviation of residuals
 

# t likelihood
mod3_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dt( mu[i], tau, df )
mu[i] = b[1] + b[2]*log_income[i] + b[3]*is_oil[i]
}

for (i in 1:3) {
b[i] ~ dnorm(0.0, 1.0/1.0e6)
}

df = nu + 2.0 # we want degrees of freedom > 2 to guarantee existence of mean and variance
nu ~ dexp(1.0)

tau ~ dgamma(5/2.0, 5*10.0/2.0) # tau is close to, but not equal to the precision
sig = sqrt( 1.0 / tau * df / (df - 2.0) ) # standard deviation of errors
} "

#Compare models using Deviance Information Criterion
dic.samples(mod1, n.iter=1e3)
## Mean deviance:  231.2 
## penalty 2.911 
## Penalized deviance: 234.1

dic.samples(mod2, n.iter=1e3)


###ANOVA
data("PlantGrowth")
?PlantGrowth
head(PlantGrowth)

boxplot(weight ~ group, data=PlantGrowth)

lmod = lm(weight ~ group, data=PlantGrowth)
summary(lmod)

anova(lmod)

library("rjags")
## Loading required package: coda
## Linked to JAGS 4.2.0
## Loaded modules: basemod,bugs
mod_string = " model {
for (i in 1:length(y)) {
y[i] ~ dnorm(mu[grp[i]], prec)
}

for (j in 1:3) {
mu[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*1.0/2.0)
sig = sqrt( 1.0 / prec )
} "

set.seed(82)
str(PlantGrowth)
data_jags = list(y=PlantGrowth$weight, 
              grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
    inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains

plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
effectiveSize(mod_sim)

(pm_params = colMeans(mod_csim))

yhat = pm_params[1:3][data_jags$grp]
resid = data_jags$y - yhat
plot(resid)
plot(yhat, resid)

summary(mod_sim)

HPDinterval(mod_csim)

mean(mod_csim[,3] > mod_csim[,1])
## [1] 0.9397333
mean(mod_csim[,3] > 1.1*mod_csim[,1])



##########Day 4
## Logistic regression  
library("boot")
data("urine")
?urine
head(urine)
dat = na.omit(urine)
pairs(dat)

library("corrplot")
Cor = cor(dat)
corrplot(Cor, type="upper", method="ellipse", tl.pos="d")
corrplot(Cor, type="lower", method="number", col="black", 
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")


X = scale(dat[,-1], center=TRUE, scale=TRUE)
head(X[,"gravity"])

colMeans(X)

apply(X, 2, sd)

ddexp = function(x, mu, tau) {
    0.5*tau*exp(-tau*abs(x-mu)) 
}
curve(ddexp(x, mu=0.0, tau=1.0), from=-5.0, to=5.0, ylab="density", main="Double exponential\ndistribution") # double exponential distribution
curve(dnorm(x, mean=0.0, sd=1.0), from=-5.0, to=5.0, lty=2, add=TRUE) # normal distribution
legend("topright", legend=c("double exponential", "normal"), lty=c(1,2), bty="n")

library("rjags")
## Loading required package: coda
## Linked to JAGS 4.2.0
## Loaded modules: basemod,bugs
mod1_string = " model {
for (i in 1:length(y)) {
y[i] ~ dbern(p[i])
logit(p[i]) = int + b[1]*gravity[i] + b[2]*ph[i] + b[3]*osmo[i] + b[4]*cond[i] + b[5]*urea[i] + b[6]*calc[i]
}
int ~ dnorm(0.0, 1.0/25.0)
for (j in 1:6) {
b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
}
} "

set.seed(92)
head(X)

data_jags = list(y=dat$r, gravity=X[,"gravity"], ph=X[,"ph"], osmo=X[,"osmo"], cond=X[,"cond"], urea=X[,"urea"], calc=X[,"calc"])

params = c("int", "b")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

## convergence diagnostics
plot(mod1_sim, ask=TRUE)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)

summary(mod1_sim)

par(mfrow=c(3,2))
densplot(mod1_csim[,1:6], xlim=c(-3.0, 3.0))

mod2_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dbern(p[i])
logit(p[i]) = int + b[1]*gravity[i] + b[2]*cond[i] + b[3]*calc[i]
}
int ~ dnorm(0.0, 1.0/25.0)
for (j in 1:3) {
b[j] ~ dnorm(0.0, 1.0/25.0) # noninformative for logistic regression
}
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)

update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim, ask=TRUE)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

dic2 = dic.samples(mod2, n.iter=1e3)


dic1

dic2

summary(mod2_sim)

HPDinterval(mod2_csim)

par(mfrow=c(3,1))
densplot(mod2_csim[,1:3], xlim=c(-3.0, 3.0))

(pm_coef = colMeans(mod2_csim))

pm_Xb = pm_coef["int"] + X[,c(1,4,6)] %*% pm_coef[1:3]
phat = 1.0 / (1.0 + exp(-pm_Xb))
head(phat)


(tab0.5 = table(phat > 0.5, data_jags$y))

sum(diag(tab0.5)) / sum(tab0.5)

(tab0.3 = table(phat > 0.3, data_jags$y))

sum(diag(tab0.3)) / sum(tab0.3)


#Poisson regression

library("COUNT")
## Loading required package: msme
## Loading required package: MASS
## Loading required package: lattice
## Loading required package: sandwich
data("badhealth")
?badhealth
head(badhealth)

any(is.na(badhealth))
## [1] FALSE
As usual, let’s visualize these data.

hist(badhealth$numvisit, breaks=20)
plot(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==0, xlab="age", ylab="log(visits)")
points(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==1, col="red")

library("rjags")

mod_string = " model {
    for (i in 1:length(numvisit)) {
numvisit[i] ~ dpois(lam[i])
log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
}

int ~ dnorm(0.0, 1.0/1e6)
b_badh ~ dnorm(0.0, 1.0/1e4)
b_age ~ dnorm(0.0, 1.0/1e4)
b_intx ~ dnorm(0.0, 1.0/1e4)
} "

set.seed(102)

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)


X = as.matrix(badhealth[,-1])
X = cbind(X, with(badhealth, badh*age))
head(X)


(pmed_coef = apply(mod_csim, 2, median))


llam_hat = pmed_coef["int"] + X %*% pmed_coef[c("b_badh", "b_age", "b_intx")]
lam_hat = exp(llam_hat)

hist(lam_hat)



resid = badhealth$numvisit - lam_hat
plot(resid) # the data were ordered

plot(lam_hat, badhealth$numvisit)
abline(0.0, 1.0)


plot(lam_hat[which(badhealth$badh==0)], resid[which(badhealth$badh==0)], xlim=c(0, 8), ylab="residuals", xlab=expression(hat(lambda)), ylim=range(resid))
points(lam_hat[which(badhealth$badh==1)], resid[which(badhealth$badh==1)], col="red")


var(resid[which(badhealth$badh==0)])

var(resid[which(badhealth$badh==1)])


summary(mod_sim)


#Hierarchical modeling

dat = read.table(file="cookies.dat", header=TRUE)
head(dat)

table(dat$location)

hist(dat$chips)

boxplot(chips ~ location, data=dat)

#Prior predictive checks

set.seed(112)
n_sim = 500
alpha_pri = rexp(n_sim, rate=1.0/2.0)
beta_pri = rexp(n_sim, rate=5.0)
mu_pri = alpha_pri/beta_pri
sig_pri = sqrt(alpha_pri/beta_pri^2)

summary(mu_pri)

summary(sig_pri)

lam_pri = rgamma(n=n_sim, shape=alpha_pri, rate=beta_pri)
summary(lam_pri)
(lam_pri = rgamma(n=5, shape=alpha_pri[1:5], rate=beta_pri[1:5]))

(y_pri = rpois(n=150, lambda=rep(lam_pri, each=30)))

library("rjags")
## Loading required package: coda
## Linked to JAGS 4.2.0
## Loaded modules: basemod,bugs
mod_string = " model {
for (i in 1:length(chips)) {
chips[i] ~ dpois(lam[location[i]])
}

for (j in 1:max(location)) {
lam[j] ~ dgamma(alpha, beta)
}

alpha = mu^2 / sig^2
beta = mu / sig^2

mu ~ dgamma(2.0, 1.0/5.0)
sig ~ dexp(1.0)

} "

set.seed(113)

data_jags = as.list(dat)

params = c("lam", "mu", "sig")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)


(pm_params = colMeans(mod_csim))

yhat = rep(pm_params[1:5], each=30)
resid = dat$chips - yhat
plot(resid)

plot(jitter(yhat), resid)

var(resid[yhat<7])
## [1] 6.447126
var(resid[yhat>11])

## location level residuals
lam_resid = pm_params[1:5] - pm_params["mu"]
plot(lam_resid)
abline(h=0, lty=2)

summary(mod_sim)
(n_sim = nrow(mod_csim))
## [1] 15000
lam_pred = rgamma(n=n_sim, shape=mod_csim[,"mu"]^2/mod_csim[,"sig"]^2, 
                  rate=mod_csim[,"mu"]/mod_csim[,"sig"]^2)
hist(lam_pred)

mean(lam_pred > 15)

y_pred = rpois(n=n_sim, lambda=lam_pred)
hist(y_pred)

mean(y_pred > 15)

hist(dat$chips)

y_pred1 = rpois(n=n_sim, lambda=mod_csim[,"lam[1]"])
hist(y_pred1)

library("car")
data("Leinhardt")
?Leinhardt
str(Leinhardt)
pairs(Leinhardt)
head(Leinhardt)

dat = na.omit(Leinhardt)
dat$logincome = log(dat$income)
dat$loginfant = log(dat$infant)
str(dat)

library("rjags")

mod_string = " model {
for (i in 1:length(y)) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = a[region[i]] + b[1]*log_income[i] + b[2]*is_oil[i]
}

for (j in 1:max(region)) {
a[j] ~ dnorm(a0, prec_a)
}

a0 ~ dnorm(0.0, 1.0/1.0e6)
prec_a ~ dgamma(1/2.0, 1*10.0/2.0)
tau = sqrt( 1.0 / prec_a )

for (j in 1:2) {
b[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*10.0/2.0)
sig = sqrt( 1.0 / prec )
} "

set.seed(116)
data_jags = list(y=dat$loginfant, log_income=dat$logincome,
                 is_oil=as.numeric(dat$oil=="yes"), region=as.numeric(dat$region))
data_jags$is_oil
table(data_jags$is_oil, data_jags$region)

params = c("a0", "a", "b", "sig", "tau")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3) # burn-in

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)

mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combine multiple chains

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)


dic.samples(mod, n.iter=1e3)

summary(mod_sim)


############GLM with Bayesian approach
library(boot)
data("urine")
head(urine)
dat<- na.omit(urine)
dim(dat)
pairs(dat)
X = scale(dat[,-1],center = T,scale = T)
X
apply(X,2,sd)
library(rjags)
ddexp = function(x, mu, tau) {
  0.5*tau*exp(-tau*abs(x-mu)) 
}
curve(ddexp(x, mu=0.0, tau=1.0), from=-5.0, to=5.0, ylab="density", main="Double exponential\ndistribution") # double exponential distribution
curve(dnorm(x, mean=0.0, sd=1.0), from=-5.0, to=5.0, lty=2, add=TRUE) # normal distribution
legend("topright", legend=c("double exponential", "normal"), lty=c(1,2), bty="n")

mod1_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*ph[i] + b[3]*osmo[i] + b[4]*cond[i] + b[5]*urea[i] + b[6]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:6) {
        b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
    }
} "

set.seed(92)
head(X)

data_jags = list(y=dat$r, gravity=X[,"gravity"], ph=X[,"ph"], osmo=X[,"osmo"], cond=X[,"cond"], urea=X[,"urea"], calc=X[,"calc"])

params = c("int", "b")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

## convergence diagnostics
plot(mod1_sim, ask=TRUE)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)
summary(mod1_sim)

par(mfrow=c(3,2))
densplot(mod1_csim[,1:6], xlim=c(-3.0, 3.0))
colnames(X)

mod2_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*cond[i] + b[3]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:3) {
        b[j] ~ dnorm(0.0, 1.0/25.0) # noninformative for logistic regression
    }
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)

update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim, ask=TRUE)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

dic2 = dic.samples(mod2, n.iter=1e3)
summary(mod2_sim)

################################################################################
n      <- 20
Y      <- 4
a      <- 3
b      <- 1

library(rjags)
#The model specification
model_string <- "model{

  # Likelihood
  Y ~ dbinom(theta,n)

  # Prior
  theta ~ dbeta(a, b)
}"

#Compile the model in JAGS
model <- jags.model(textConnection(model_string), 
                    data = list(Y=Y,n=n,a=a,b=b))
#Draw samples
update(model, 10000, progress.bar="none"); # Burnin for 10000 samples

samp <- coda.samples(model, 
                     variable.names=c("theta"), 
                     n.iter=20000, progress.bar="none")

summary(samp)

plot(samp)


#######################################
#poisson model
N      <- 20
Y      <- 11
a      <- 0.5
b      <- 0.5



model_string_pois <- "model{
  # Likelihood (can't have formulas in distribution functions)
  Y ~ dpois(mu)
  mu <- N*lamda
  # Prior
  lamda ~dgamma(a,b)
}"

#Compile the model in JAGS
model.pois <- jags.model(textConnection(model_string_pois),
                         data = list(Y=Y,N=N,a=a,b=b))


#Draw samples
update(model.pois, 10000, progress.bar="none") # Burnin for 10000 samples



samp.pois <- coda.samples(model.pois, 
                     variable.names=c("lamda"), 
                     n.iter=20000, progress.bar="none")

summary(samp.pois)

plot(samp.pois)



################################################################################
################ocupancy model
setwd("C:/Users/lihai/Documents/GitHub/Bayesian Study/OccupancyData")
library(rjags)
library(vegan)
Paguma_larvata_inner <-read.table("Paguma_larvata_hist_inner.txt",head=TRUE)
str(Paguma_larvata_inner)

CBLcovs80<-read.table("CBLcovs80.txt",head=TRUE)

CBLcov80.std<-decostand(CBLcovs80,method="standardize")

Pagumalarvata80<-unmarkedFrameOccu(y=Paguma_larvata_inner,siteCovs= CBLcov80.std)

###Occupancy model
#1.Specify the model.

occu_model_string <- " model {
  for (i in 1:nsites) {
    z[i] ~ dbern(psi)
    p.eff[i] <- z[i] * p
    for (j in 1:nvisits) {
      y[i,j] ~ dbern(p.eff[i])
    } #j
  } #i
  # Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  # Derived quantities
  occ.fs <- sum(z[])
} "


# 2.Set up the model.


# 3. Run the MCMC sampler.


# 4. Post processing.









































































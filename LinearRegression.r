# http://www.math.uah.edu/stat/data/Challenger2.txt
# 23 previous space shuttle launches before the Challenger disaster
# T is the temperature in Fahrenheit, I is the O-ring damage index

oring=read.table("http://www.math.uah.edu/stat/data/Challenger2.txt",header=T)
attach(oring)
#note: masking T=TRUE 

plot(T,I)

oring.lm=lm(I~T)
summary(oring.lm)

# add fitted line to scatterplot
lines(T,fitted(oring.lm))            
# 95% posterior interval for the slope
-0.24337 - 0.06349*qt(.975,21)
-0.24337 + 0.06349*qt(.975,21)
# note that these are the same as the frequentist confidence intervals

# the Challenger launch was at 31 degrees Fahrenheit
# how much o-ring damage would we predict?
# y-hat
18.36508-0.24337*31
coef(oring.lm)
coef(oring.lm)[1] + coef(oring.lm)[2]*31  

# posterior prediction interval (this implies an non-informative prior, the result is same as frequentist)
predict(oring.lm,data.frame(T=31),interval="predict")  
10.82052-2.102*qt(.975,21)*sqrt(1+1/23+((31-mean(T))^2/22/var(T)))

# posterior probability that damage index is greater than zero
1-pt((0-10.82052)/(2.102*sqrt(1+1/23+((31-mean(T))^2/22/var(T)))),21)


# http://www.math.uah.edu/stat/data/Galton.txt
# Galton's seminal data on predicting the height of children from the 
# heights of the parents, all in inches

heights=read.table("http://www.math.uah.edu/stat/data/Galton.txt",header=T)
attach(heights)
names(heights)

pairs(heights)
summary(lm(Height~Father+Mother+Gender+Kids))
summary(lm(Height~Father+Mother+Gender))
heights.lm=lm(Height~Father+Mother+Gender)

# each extra inch taller a father is is correlated with 0.4 inch extra height in the child
# each extra inch taller a mother is is correlated with 0.3 inch extra height in the child
# a male child is on average 5.2 inches taller than a female child
# 95% posterior interval for the the difference in height by gender
5.226 - 0.144*qt(.975,894)
5.226 + 0.144*qt(.975,894)

# posterior prediction interval (same as frequentist)
predict(heights.lm,data.frame(Father=68,Mother=64,Gender="M"),interval="predict")
predict(heights.lm,data.frame(Father=68,Mother=64,Gender="F"),interval="predict")






#-------------------------------
#Using JAG
#Get data
library("car")
data("Leinhardt")
?Leinhardt
head(Leinhardt)
str(Leinhardt)
pairs(Leinhardt)

#--a simple linear regression model that relates infant mortality to per capita income.
plot(infant ~ income, data=Leinhardt)
hist(Leinhardt$infant)
hist(Leinhardt$income)
Leinhardt$loginfant = log(Leinhardt$infant)
Leinhardt$logincome = log(Leinhardt$income)

plot(loginfant ~ logincome, data=Leinhardt)

#Model implementation
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
                        n.iter=5000)

mod1_csim = do.call(rbind, mod1_sim) # combine multiple chains

#MCMC convergence check
plot(mod1_sim)
gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

#posterior summary of the parameters in our model.
summary(mod1_sim)

#Residual checks
X = cbind(rep(1.0, data1_jags$n), data1_jags$log_income)
head(X)

(pm_params1 = colMeans(mod1_csim)) # posterior mean

yhat1 = drop(X %*% pm_params1[1:2])
resid1 = data1_jags$y - yhat1
plot(resid1) # against data index

plot(yhat1, resid1) # against predicted values

qqnorm(resid1) # checking normality of residuals


#--Multiple covariates
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

X2 = cbind(rep(1.0, data1_jags$n), data2_jags$log_income, data2_jags$is_oil)
head(X2)

(pm_params2 = colMeans(mod2_csim)) # posterior mean

yhat2 = drop(X2 %*% pm_params2[1:3])
resid2 = data2_jags$y - yhat2
plot(resid2) # against data index

plot(yhat2, resid2) # against predicted values
plot(yhat1, resid1) # residuals from the first model
sd(resid2) # standard deviation of residuals


#--Change the likelihood dist to t to accommodate the outliers
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

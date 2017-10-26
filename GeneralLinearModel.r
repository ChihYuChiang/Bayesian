#-------------------------------
#Logistic

#--Data
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


#--Variable selection
#Standardize
X = scale(dat[,-1], center=TRUE, scale=TRUE)
head(X[,"gravity"])

colMeans(X)

apply(X, 2, sd)

#Model
library("rjags")
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

colnames(X) # variable names


#--Logistic model
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

colnames(X)[c(1,4,6)] # variable names


#--Prediciton
(pm_coef = colMeans(mod2_csim))

pm_Xb = pm_coef["int"] + X[,c(1,4,6)] %*% pm_coef[1:3]
phat = 1.0 / (1.0 + exp(-pm_Xb))
head(phat)

plot(phat, jitter(dat$r))

(tab0.5 = table(phat > 0.5, data_jags$y))
sum(diag(tab0.5)) / sum(tab0.5)

(tab0.3 = table(phat > 0.3, data_jags$y))
sum(diag(tab0.3)) / sum(tab0.3)




#-------------------------------
#Poisson

#--Data
library("COUNT")

data("badhealth")
?badhealth
head(badhealth)

any(is.na(badhealth))

hist(badhealth$numvisit, breaks=20)

plot(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==0, xlab="age", ylab="log(visits)")
points(jitter(log(numvisit)) ~ jitter(age), data=badhealth, subset=badh==1, col="red")


#--Model
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


#--Model checking
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


#--Results
summary(mod_sim)


#--Predictive distribution
'
we have two people aged 35, one in good health and the other in poor health.
What is the posterior probability that the individual with poor health will have more doctor visits?
This goes beyond the posterior probabilities we have calculated comparing expected responses in previous lessons.
Here we will create Monte Carlo samples for the responses themselves.
This is done by taking the Monte Carlo samples of the model parameters, and for each of those, drawing a sample from the likelihood.
'
x1 = c(0, 35, 0) # good health
x2 = c(1, 35, 35) # bad health

head(mod_csim)

loglam1 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x1
loglam2 = mod_csim[,"int"] + mod_csim[,c(2,1,3)] %*% x2

lam1 = exp(loglam1)
lam2 = exp(loglam2)

(n_sim = length(lam1))

y1 = rpois(n=n_sim, lambda=lam1)
y2 = rpois(n=n_sim, lambda=lam2)

plot(table(factor(y1, levels=0:18))/n_sim, pch=2, ylab="posterior prob.", xlab="visits")
points(table(y2+0.1)/n_sim, col="red")

mean(y2 > y1)
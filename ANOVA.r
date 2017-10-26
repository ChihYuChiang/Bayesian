data("PlantGrowth")
?PlantGrowth
head(PlantGrowth)
boxplot(weight ~ group, data=PlantGrowth)

library("rjags")
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


#--Model checking
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
effectiveSize(mod_sim)

(pm_params = colMeans(mod_csim))

yhat = pm_params[1:3][data_jags$grp]
resid = data_jags$y - yhat
plot(resid)

plot(yhat, resid)


#--Results
summary(mod_sim)
HPDinterval(mod_csim) #The HPDinterval() function in the coda package calculates intervals of highest posterior density for each parameter.

#There is a high posterior probability that the mean yield for treatment 2 is greater than the mean yield for the control group.
mean(mod_csim[,3] > mod_csim[,1])

#We have about 50/50 odds that adopting treatment 2 would increase mean yield by at least 10%.
mean(mod_csim[,3] > 1.1*mod_csim[,1])



#-----------------
#Two way ANOVA

#--Data
data("warpbreaks")
?warpbreaks
head(warpbreaks)

table(warpbreaks$wool, warpbreaks$tension)
boxplot(breaks ~ wool + tension, data=warpbreaks)
boxplot(log(breaks) ~ wool + tension, data=warpbreaks)


#--Two way additive model (No interaction)
X = model.matrix( ~ wool + tension, data=warpbreaks)
head(X)
tail(X)

mod2_string = " model {
    for( i in 1:length(y)) {
y[i] ~ dnorm(mu[i], prec)
mu[i] = int + alpha*isWoolB[i] + beta[1]*isTensionM[i] + beta[2]*isTensionH[i]
}

int ~ dnorm(0.0, 1.0/1.0e6)
alpha ~ dnorm(0.0, 1.0/1.0e6)
for (j in 1:2) {
beta[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(3/2.0, 3*1.0/2.0)
sig = sqrt(1.0 / prec)
} "

data2_jags = list(y=log(warpbreaks$breaks), isWoolB=X[,"woolB"], isTensionM=X[,"tensionM"], isTensionH=X[,"tensionH"])

params2 = c("int", "alpha", "beta", "sig")

mod2 = jags.model(textConnection(mod2_string), data=data2_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5e3)

## convergene diagnostics
plot(mod2_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)

#Summary
summary(mod2_sim)
(dic2 = dic.samples(mod2, n.iter=1e3))
dic1

boxplot(log(breaks) ~ wool + tension, data=warpbreaks)


#--Two-way cell means model (with interaction effect)
mod3_string = " model {
    for( i in 1:length(y)) {
        y[i] ~ dnorm(mu[woolGrp[i], tensGrp[i]], prec)
    }
    
    for (j in 1:max(woolGrp)) {
        for (k in 1:max(tensGrp)) {
            mu[j,k] ~ dnorm(0.0, 1.0/1.0e6)
        }
    }
    
    prec ~ dgamma(3/2.0, 3*1.0/2.0)
    sig = sqrt(1.0 / prec)
} "

str(warpbreaks)

data3_jags = list(y=log(warpbreaks$breaks), woolGrp=as.numeric(warpbreaks$wool), tensGrp=as.numeric(warpbreaks$tension))

params3 = c("mu", "sig")

mod3 = jags.model(textConnection(mod3_string), data=data3_jags, n.chains=3)
update(mod3, 1e3)

mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=5e3)
mod3_csim = as.mcmc(do.call(rbind, mod3_sim))

plot(mod3_sim, ask=TRUE)

## convergence diagnostics
gelman.diag(mod3_sim)
autocorr.diag(mod3_sim)
effectiveSize(mod3_sim)
raftery.diag(mod3_sim)

#Summary
(dic3 = dic.samples(mod3, n.iter=1e3))
dic1

summary(mod3_sim)

HPDinterval(mod3_csim)

par(mfrow=c(3,2)) # arrange frame for plots
densplot(mod3_csim[,1:6], xlim=c(2.0, 4.5))

prop.table( table( apply(mod3_csim[,1:6], 1, which.min) ) )
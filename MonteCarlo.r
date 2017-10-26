#--Simulation
#Paras
set.seed(32) # Initializes the random number generator so we can replicate these results. To get different random numbers, change the seed.
m = 1e4
a = 2.0
b = 1.0 / 3.0

#Draw
theta = rgamma(n=m, shape=a, rate=b)


#--Observations (comparing theoretical result and sample result)
#Plot
hist(theta, freq=FALSE)
curve(dgamma(x=x, shape=a, rate=b), col="blue", add=TRUE)

#Mean
mean(theta) # sample mean
a / b # true expected value

#Var
var(theta) # sample variance
a / b^2 # true variance of Gamma(a,b)

#Prob that theta < 5
ind = theta < 5.0 # set of indicators, TRUE if theta_i < 5
mean(ind) # automatically converts FALSE/TRUE to 0/1

pgamma(q=5.0, shape=a, rate=b) # true value of Pr( theta < 5 )

#Quantile
quantile(x=theta, probs=0.9)
qgamma(p=0.9, shape=a, rate=b) # true value of 0.9 quantile


#--Monte Carlo error
#For entire dist
se = sd(theta) / sqrt(m)
2.0 * se # we are reasonably confident that the Monte Carlo estimate is no more than t his far from the truth

#For specific part of dist
ind = theta < 5.0
se = sd(ind) / sqrt(m)
2.0 * se # we are reasonably confident that the Monte Carlo estimate is no more than t his far from the truth




#--Simulation (hierarchical model)
m = 10e4
y = numeric(m) # create the vectors we will fill in with simulations
phi = numeric(m)
for (i in 1:m) {
  phi[i] = rbeta(n=1, shape1=2.0, shape2=2.0)
  y[i] = rbinom(n=1, size=10, prob=phi[i])
}
# which is equivalent to the following 'vectorized' code
phi = rbeta(n=m, shape1=2.0, shape2=2.0)
y = rbinom(n=m, size=10, prob=phi)

#Marginal dist y
mean(y)
plot(prop.table(table(y)), ylab="P(y)", main="Marginal distribution of y")

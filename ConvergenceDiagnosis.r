#--Log function to compute log likelihood (for the joint distribution)
lg = function(mu, n, ybar) {
  mu2 = mu^2
  n * (ybar * mu - mu2 / 2.0) - log(1 + mu2)
}

#--Function to execute the Random-Walk Metropolis-Hastings sampler with normal proposals.
mh = function(n, ybar, n_iter, mu_init, cand_sd) {
  ## Random-Walk Metropolis-Hastings algorithm
  
  ## step 1, initialize
  mu_out = numeric(n_iter)
  accpt = 0
  mu_now = mu_init
  lg_now = lg(mu=mu_now, n=n, ybar=ybar)
  
  ## step 2, iterate
  for (i in 1:n_iter) {
    ## step 2a
    mu_cand = rnorm(n=1, mean=mu_now, sd=cand_sd) # draw a candidate
    
    ## step 2b
    lg_cand = lg(mu=mu_cand, n=n, ybar=ybar) # evaluate log of g with the candidate
    lalpha = lg_cand - lg_now # log of acceptance ratio
    alpha = exp(lalpha)
    
    ## step 2c
    u = runif(1) # draw a uniform variable which will be less than alpha with probability min(1, alpha)
    if (u < alpha) { # then accept the candidate
      mu_now = mu_cand
      accpt = accpt + 1 # to keep track of acceptance
      lg_now = lg_cand
    }
    
    ## collect results
    mu_out[i] = mu_now # save this iteration's value of mu
  }
  
  ## return a list of output
  list(mu=mu_out, accpt=accpt/n_iter)
}

#Here is an example of a chain that has most likely converged.
set.seed(61)
post0 = mh(n=n, ybar=ybar, n_iter=10e3, mu_init=0.0, cand_sd=0.9)
coda::traceplot(as.mcmc(post0$mu[-c(1:500)]))

#If the chain is stationary, it should not be showing any long-term trends. The average value for the chain should be roughly flat. It should not be wandering as in this example:
set.seed(61)
post1 = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=0.04)
coda::traceplot(as.mcmc(post1$mu[-c(1:500)]))

#If this is the case, you need to run the chain many more iterations, as seen here:
set.seed(61)
post2 = mh(n=n, ybar=ybar, n_iter=100e3, mu_init=0.0, cand_sd=0.04)
coda::traceplot(as.mcmc(post2$mu))


#--Autocorrelation to decide sample size
coda::autocorr.plot(as.mcmc(post0$mu))
coda::autocorr.diag(as.mcmc(post0$mu))

coda::autocorr.plot(as.mcmc(post1$mu))
coda::autocorr.diag(as.mcmc(post1$mu))

str(post2) # contains 100,000 iterations
coda::effectiveSize(as.mcmc(post2$mu)) # effective sample size of ~350

## thin out the samples until autocorrelation is essentially 0. This will leave you with approximately independent samples. The number of samples remaining is similar to the effective sample size.
coda::autocorr.plot(as.mcmc(post2$mu), lag.max=500)

thin_interval = 400 # how far apart the iterations are for autocorrelation to be essentially 0.
thin_indx = seq(from=thin_interval, to=length(post2$mu), by=thin_interval)
head(thin_indx)

post2mu_thin = post2$mu[thin_indx]
traceplot(as.mcmc(post2$mu))
traceplot(as.mcmc(post2mu_thin))
coda::autocorr.plot(as.mcmc(post2mu_thin), lag.max=10)
effectiveSize(as.mcmc(post2mu_thin))
length(post2mu_thin)

str(post0) # contains 10,000 iterations
coda::effectiveSize(as.mcmc(post0$mu)) # effective sample size of ~2,500

#Raftery and Lewis diagnostic for effective sample size
raftery.diag(as.mcmc(post0$mu))
raftery.diag(as.mcmc(post0$mu), q=0.005, r=0.001, s=0.95)


#--Burn-in
set.seed(62)
post3 = mh(n=n, ybar=ybar, n_iter=500, mu_init=10.0, cand_sd=0.3)
coda::traceplot(as.mcmc(post3$mu))


#--Multiple chains, Gelman-Rubin
set.seed(61)

nsim = 500
post1 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=15.0, cand_sd=0.4)
post1$accpt
post2 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=-5.0, cand_sd=0.4)
post2$accpt
post3 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=7.0, cand_sd=0.1)
post3$accpt
post4 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=23.0, cand_sd=0.5)
post4$accpt
post5 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=-17.0, cand_sd=0.4)
post5$accpt
pmc = mcmc.list(as.mcmc(post1$mu), as.mcmc(post2$mu), 
                as.mcmc(post3$mu), as.mcmc(post4$mu), as.mcmc(post5$mu))
str(pmc)

coda::traceplot(pmc)
coda::gelman.diag(pmc)
coda::gelman.plot(pmc)
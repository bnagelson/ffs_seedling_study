
## FFS Seedling Study 
## JAGS model scripts


## dens.mod2 ####

dens.mod2 <- "
model { #### CHANGE THIS TO 'MODEL' BEFORE RUNNING!!!!!! also need to un-comment the parameter estimates

## Model ### ----------------------------------
for(plot in 1:N){
  observed[plot] ~ dnegbin(p[plot], spread)
  p[plot] <- spread/(spread+mu[plot])     


  #expected
  log(mu[plot]) <- intercept +
                   trt.eff[trt[plot]] +                           #FE for treatment (intercept)
                   stand.eff[stand[plot],trt[plot]] +              #RE for stand
                   timesince.eff[trt[plot]]*ztimesince[plot]+       #FE for time-since-treatment (varies by treatment)
                   #timesince.eff.sq[trt[plot]]*ztimesince[plot]^2 + #quadratic for time-since-treatment
                   inprod(beta[], z.model.matrix[plot,]) +                #the model matrix for continuous vars
                   log(plotsize[plot])                             #offset term for plot size (sampling area)
}

for(plot in 1:N){
# #For posterior predictive check
     predicted[plot] ~ dnegbin(p[plot], spread)
#   #E_obs[plot] <- (observed[plot]-mu[plot])   #observed-expected
#   #E_pred[plot] <- (predicted[plot]-mu[plot]) #predicted-expected
#   #SE_obs[plot] <- E_obs[plot]^2 #squared errors
#   #SE_pred[plot] <- E_pred[plot]^2
}

# Priors Model 1 ##
spread ~ dunif(0,5)


### treatment, stand, and plot effects - using post-sweeping as per Ogle & Barber (2020)
for(t in 1:ntrt){
  for(s in 1:nstands){
    stand.eff[s,t] ~ dnorm(0,stand.tau) #don't monitor
    stand.eff.star[s,t] <- stand.eff[s,t] - ave.stand.eff[t]  #monitor this!!!
  }
  #compute mean of all stands within treatments
  ave.stand.eff[t] <- mean(stand.eff[,t])
}
#average of all stand effects
ave.ave.stand.eff <- mean(ave.stand.eff[]) 

#treatment effects
for(t in 1:ntrt) {
  trt.eff[t] ~ dnorm(0, 0.01)
  trt.eff.star[t] <- trt.eff[t] + ave.stand.eff[t] - ave.trt.eff - ave.ave.stand.eff #monitor this!!!
}
#average of all treatment effects
ave.trt.eff <- mean(trt.eff[])

### overal intercept
intercept ~ dnorm(0, 1e-6)
intercept.star <- intercept + ave.trt.eff + ave.ave.stand.eff
 
#prior for variance of stand-level random effects
stand.tau <- 1/(stand.sigma^2)
stand.sigma ~ dunif(1e-3,1e3)


 ### for continuous, non-hierarchical vars
 for(k in 1:nbeta){
   beta[k] ~ dnorm(0, 0.01)
 }

 #time since treatment - slope varies by treatment
 for(x in 1:ntrt){
    timesince.eff[x] ~ dnorm(0, .01)
    #timesince.eff.sq[x] ~ dnorm(0, .01)
 }

#end model
}

# model {
# fake <- 0
# }
"
#save script
writeLines(dens.mod2, con = file.path("dens.mod2.txt"))


### ========================================
### Occupancy Models
### ========================================


occ.mod2 <- " #Line 1, for debugging
model { 
  for (plot in 1:N){
    for(sample in 2:nSample[plot]){
  #logit function for probability of colonization  
      logit(pcol[plot,sample])    <-       intercept.c + 
                                           trteff.c[trt[plot]]+                              #treatment effect 
                                           plot.type.c[type[plot]] +  #random effect for plot type
                                           locationeff.c[stand[plot],trt[plot],location[plot]] + #RE for plot location
                                           standeff.c[stand[plot],trt[plot]]+                 #random effect for stand within treatment
                                           timesinceeff.c[trt[plot]]*timesince[plot,sample]+ #time since the most recent treatment
    		                                   inprod(beta.c[],betavar.array[plot,sample,])
    		                                   
  #logit function for probability of extinction (inverse of peristence)                                   
      logit(pext[plot,sample])    <-        intercept.e +
                                            trteff.e[trt[plot]]+                              #treatment effect 
                                            plot.type.e[type[plot]]+
                                            locationeff.e[stand[plot],trt[plot],location[plot]] + #RE for plot location
                                            standeff.e[stand[plot],trt[plot]]+                #random effect for stand within treatment
                                            timesinceeff.e[trt[plot]]*timesince[plot,sample]+ #time since the most recent treatment
                                            inprod(beta.e[],betavar.array[plot,sample,])
                                          
  #on/off switch to determine whether to use pcol or pext, based on occupancy of last sample                               
      pocc[plot,sample] <- y[plot,sample-1]*pext[plot,sample]+(1-y[plot,sample-1])*pcol[plot,sample]   # prob of occupancy this sample
      y[plot,sample] ~ dbern(pocc[plot,sample])
      
      
  #goodness-of-fit with pearsons residuals
      Pears_resid[plot,sample] <- (y[plot,sample] - pocc[plot,sample])/sqrt(pocc[plot,sample]*(1-pocc[plot,sample]))
      predicted[plot,sample] ~ dbern(pocc[plot,sample]) #predictions
      Pears_resid_new[plot,sample] <- (predicted[plot,sample] - pocc[plot,sample])/sqrt(pocc[plot,sample]*(1-pocc[plot,sample]))
      D[plot,sample] <- pow(Pears_resid[plot,sample], 2)
      D_new[plot,sample] <- pow(Pears_resid_new[plot,sample], 2)
    
      
      ## WOuld like to do this, but there is no 'mu' value... hmmm
  #      #For posterior predictive check
  # E_obs[plot,sample] <- (y[plot,sample]-mu[plot,sample])   #observed-expected
  # E_pred[plot,sample] <- (yP[plot,sample]-mu[plot,sample]) #predicted-expected
  # SE_obs[plot,sample] <- E_obs[plot,sample]^2 #squared errors
  # SE_pred[plot,sample] <- E_pred[plot,sample]^2
    }
  }
  
  # SSE_obs = sum(SE_obs[,])
  # SSE_pred = sum(SE_pred[,])

  #sum within rows to avoid dealing with NAs
  for (plot in 1:N){
    D_rows[plot] <- sum(D[plot,2:nSample[plot]])
    Dnew_rows[plot] <- sum(D_new[plot,2:nSample[plot]])
  }
  
  fit <- sum(D_rows)
  fit_new <- sum(Dnew_rows)
  
  ###
  ### Priors
  ###
  
### plot type effect - post-sweeping
  for(p in 1:2){
    plot.type.c[p] ~ dnorm(0, pt.c.tau)
    plot.type.e[p] ~ dnorm(0, pt.e.tau)
    #monitor these
    plot.type.c.star[p] <- plot.type.c[p] - ave.plot.type.c
    plot.type.e.star[p] <- plot.type.e[p] - ave.plot.type.e
  }
  ave.plot.type.c <- mean(plot.type.c[])
  ave.plot.type.e <- mean(plot.type.e[])
  pt.c.tau <- 1/(pt.c.sd^2)
  pt.e.tau <- 1/(pt.e.sd^2)
  pt.c.sd ~ dgamma(2,.01)
  pt.e.sd ~ dgamma(2,.01)
   
### stand effect (random effect) - using post-sweeping as per Ogle & Barber (2020)
  for(t in 1:ntrt){
    trteff.c[t] ~ dnorm(0, .5)
    #indentifiable trteff = unidentifiable trteff - average trteff - average of all standeffs + average of all locationeffs[t] - average of all locationeffs
    trteff.c.star[t] <- trteff.c[t] + ave.standeff.c[t] - ave.trteff.c - ave.ave.standeff.c + ave.ave.locationeff.c[t] - ave.ave.ave.location.c
    trteff.e[t] ~ dnorm(0, 0.5)
    trteff.e.star[t] <- trteff.e[t] + ave.standeff.e[t] - ave.trteff.e - ave.ave.standeff.e + ave.ave.locationeff.e[t] - ave.ave.ave.location.e
    
    for(s in 1:(nstands)){
      standeff.c[s,t] ~ dnorm(trteff.c[t],stand.c.tau) #don't monitor
      standeff.e[s,t] ~ dnorm(trteff.e[t],stand.c.tau) #don't monitor
      standeff.c.star[s,t] <- standeff.e[s,t] - ave.standeff.e[t] + ave.locationeff.c[s,t] - ave.ave.locationeff.c[t] #monitor this!!!
      standeff.e.star[s,t] <- standeff.e[s,t] - ave.standeff.e[t] + ave.locationeff.e[s,t] - ave.ave.locationeff.e[t] #monitor this!!!
      
      for(p in 1:nplots[s,t]){
        locationeff.c[s,t,p] ~ dnorm(standeff.c[s,t], location.c.tau)
        locationeff.e[s,t,p] ~ dnorm(standeff.e[s,t], location.e.tau)
        locationeff.c.star[s,t,p] <- locationeff.c[s,t,p] - ave.locationeff.c[s,t]
        locationeff.e.star[s,t,p] <- locationeff.e[s,t,p] - ave.locationeff.e[s,t]
      }
      #compute stand-level means
      ave.locationeff.c[s,t] <- mean(locationeff.c[s,t,1:nplots[s,t]])
      ave.locationeff.e[s,t] <- mean(locationeff.e[s,t,1:nplots[s,t]])
    }
    #compute treatment-level means
    ave.standeff.c[t] <- mean(standeff.c[,t])
    ave.standeff.e[t] <- mean(standeff.e[,t])
    #compute treatment-level location mean
    ave.ave.locationeff.c[t] <- mean(ave.locationeff.c[,t])
    ave.ave.locationeff.e[t] <- mean(ave.locationeff.e[,t])
  }
  stand.c.tau <- 1/stand.c.sd^2
  stand.e.tau <- 1/stand.e.sd^2 
  stand.c.sd ~ dgamma(2,.01)
  stand.e.sd ~ dgamma(2,.01)
  location.c.tau <- 1/location.c.sd^2
  location.e.tau <- 1/location.e.sd^2 
  location.c.sd ~ dgamma(2,.01)
  location.e.sd ~ dgamma(2,.01)
  
  
  ave.trteff.c <- mean(trteff.c[])
  ave.trteff.e <- mean(trteff.e[])
  ave.ave.standeff.c <- mean(ave.standeff.c[])
  ave.ave.standeff.e <- mean(ave.standeff.e[])
  ave.ave.ave.location.c <- mean(ave.ave.locationeff.c[])
  ave.ave.ave.location.e <- mean(ave.ave.locationeff.e[])
  
### intercepts
intercept.c ~ dnorm(0, .5)
intercept.e ~ dnorm(0, .5)
# identifiable intercepts (add back in the averaged FE and REs)
intercept.c.star <- intercept.c + ave.trteff.c + ave.ave.standeff.c + ave.plot.type.c + ave.ave.ave.location.c
intercept.e.star <- intercept.e + ave.trteff.e + ave.ave.standeff.e + ave.plot.type.e + ave.ave.ave.location.e
  
### time since treatment - slope varies by treatment
  for(x in 1:ntrt){
      timesinceeff.c[x] ~ dnorm(0, .01)
      timesinceeff.e[x] ~ dnorm(0, .01)
  }


### regression coef priors - continuous variables
  for(n in 1:nbetavars){
    beta.c[n] ~ dnorm(0, .001)
  }
  
  for(n in 1:nbetavars){
    beta.e[n] ~ dnorm(0, .001)
  }
  
  # pre.occ.eff.c ~ dnorm(0, .001)
  # pre.occ.eff.e ~ dnorm(0, .001)

# #### un-transform coefficients
# beta.c[1:ncolvars] <- zbeta.c[1:ncolvars]/col.sd[1:ncolvars] 
# beta.e[1:npervars] <- zbeta.e[1:npervars]/per.sd[1:npervars] 
# timesinceeff.c[1:4]<- ztimesinceeff.c[1:4]/time.sd
# timesinceeff.e[1:4]<- ztimesinceeff.e[1:4]/time.sd
# intercept.c <- zintercept.c - sum(zbeta.c[1:ncolvars]*col.mean[1:ncolvars]/col.sd[1:ncolvars], ztimesinceeff.c[1:4]*time.mean/time.sd)
# intercept.e <- zintercept.e - sum(zbeta.e[1:npervars]*per.mean[1:npervars]/per.sd[1:npervars], ztimesinceeff.e[1:4]*time.mean/time.sd)

} # end model
"
writeLines(occ.mod2, con = file.path("occ.mod2.txt"))

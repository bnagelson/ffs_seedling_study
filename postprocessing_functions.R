###
### FFS Seedling Study
### Post-processing

packs = c('dplyr','tidyr','ggplot2','cowplot', 'DHARMa','mcmcplots', 'reshape')
lapply(packs, function(x) {
 if(!require(x, character.only = T)) install.packages(x); library(x, character.only = T) 
})

myspecies=c('BO','DF','IC','PP','SP','WF')
names(myspecies) = myspecies

#load model data
# load("C:/Users/bryan/Documents/UNR/FFS Seedling Study/bayes/RDatafiles/dens.mod3.output.rdata", verbose=T)
# load("C:/Users/bryan/Documents/UNR/FFS Seedling Study/bayes/RDatafiles/occ.mod2.output.rdata", verbose=T)
load("occ.mod.output.rdata", verbose=T)
load("dens.mod.output.rdata", verbose=T)


## function for diagnostic checks
diag.check = function(species, model, datalist) {
  mod = model[[species]]
  observed = datalist[[species]]$observed
  
  #outlist = lapply(seq_along(myspecies), function(x){
    predicted = mod$sims.list$predicted
    #observed = data$observed
    expected = mod$sims.list$mu
    
    SE_obs = (observed-expected)^2
    SE_sim = (predicted-expected)^2
    RMSE_obs = apply(SE_obs, 1, mean)
    RMSE_sim = apply(SE_sim, 1, mean)
    pval = length(which(RMSE_sim>RMSE_obs))/length(RMSE_sim)
    
    sim.pred = mod$sims.list$predicted
    sim.mu = mod$sims.list$mu
    sim=createDHARMa(simulatedResponse = t(sim.pred), 
                     observedResponse = observed,
                     fittedPredictedResponse = apply(sim.mu, 2, median),
                     integerResponse = T)
    #disp = testDispersion(sim, plot=F)
    #zi = testZeroInflation(sim, plot=F)$p.value
    #res = plotResiduals(sim)
    
    return(list(pval = pval,
                DHARMaSim=sim))
  #})
  #return(outlist)
}


bayes.r2 = function(mod, species = c(myspecies)) {
  #dev
  mod=dens.mod6.output[['DF']]
  y_pred = mod$sims.list$predicted
  var_fit = apply(y_pred, 1, var)
  var_res = apply(mod$sims.list$predicted - mod$sims.list$mu, 1, sum)^2
  r2 = var_fit/(var_fit+var_res)
  return(r2)
}

#### ---------------------
## Forest Plot function --
#### ---------------------
dens.fp.df <- function(
    model=dens.mod2.output,
    data=model.data,
    jagsdata = jags.data.list,
    modelmatrix = model.mats,
    species='PP',
    my.ci=c(.1,.9),
    coeff.type='scaled') {
  
  mod=model[[species]]
  mydata=data[[species]]
  mymod = model[[species]]
  jagsdata = jagsdata[[species]]
  mysims=mymod$sims.list
  mymeans=mymod$mean
  covariate.names=colnames(modelmatrix[[species]])
  #values for unscaling
  cont.sd = apply(modelmatrix[[species]], 2, sd) #SDs of continuous variables
  cont.mean=apply(modelmatrix[[species]], 2, mean)
  time.sd=sd(mydata$timesince)
  time.mean=mean(mydata$timesince)
  
  #scaled parameters
  zbeta.coeffs = mysims$beta
  zint = mysims$intercept
  ztime.coeffs=mysims$timesince.eff
  #ztime2.coeffs=mysims$timesince.eff.sq
  
  #terms to unscale categorical parameters
  beta.unscale=rowSums(zbeta.coeffs*(cont.mean/cont.sd))
  time.unscale=rowSums(ztime.coeffs*(time.mean/time.sd))
  #unscaled parameters
  beta.coeffs = mysims$beta %*% diag(1/cont.sd)
  time.coeffs=mysims$timesince.eff/time.sd
  int = zint - beta.unscale - time.unscale 
  trt.coeffs=mymod$sims.list$trt.eff
  
  #keep names of predictor variables with coefficient matrices
  colnames(beta.coeffs) = covariate.names
  colnames(zbeta.coeffs) = covariate.names
  colnames(trt.coeffs) = c('Control','Fire Only','Mechanical Only', 'Mechanical+Fire')
  #colnames(ztrt.coeffs) = c('Control','Fire Only','Mechanical Only', 'Mechanical+Fire')
  colnames(time.coeffs)=paste0("Time:",c('Control','Fire Only','Mechanical Only', 'Mechanical+Fire'))
  colnames(ztime.coeffs)=paste0("Time:",c('Control','Fire Only','Mechanical Only', 'Mechanical+Fire'))
  # colnames(time2.coeffs)=paste0("Time.sq:",c('Control','Fire Only','Mechanical Only', 'Mechanical+Fire'))
  # colnames(ztime2.coeffs)=paste0("Time.sq:",c('Control','Fire Only','Mechanical Only', 'Mechanical+Fire'))
  
  #package all scaled coefficients into a dataframe
  ztemp.df =data.frame(int=zint, 
                       trt.coeffs, 
                       ztime.coeffs, # the trt.coeffs don't require unscaling because they are in reference to the intercept
                       zbeta.coeffs)
  #unscaled coefficients
  temp.df =data.frame(int, 
                      trt.coeffs, 
                      time.coeffs, #time2.coeffs, 
                      beta.coeffs)
  #vector for labeling the type of coefficient
  data.type.labels=c(rep('categorical',5), rep('time',4), rep('continuous',jagsdata$nbeta))
  
  #calculate summary statistics
  zmean.df = apply(ztemp.df,2,mean)
  mean.df = apply(temp.df,2,mean)
  zmedian.df=apply(ztemp.df,2,median)
  median.df=apply(temp.df,2,median)
  zci.df = t(apply(ztemp.df, 2, function(x) quantile(x, my.ci)))
  ci.df = t(apply(temp.df, 2, function(x) quantile(x, my.ci)))
  both.ci=rbind(zci.df, ci.df)
  
  #make better term names
  term.names.clean = c('Intercept', 
  'Control', 'Fire Only','Mech Only','Mech+Fire',
                       'Time[Control]', 'Time since Fire Only', 'Time since Mech Only', 'Time since Mech+Fire',
                       colnames(modelmatrix[[species]]))
  
  #combine scaled and unscaled coeffs into single dataframe
  all.coeffs=rbind(ztemp.df, temp.df)
  coeff.type.labels=c(rep('scaled', length(zmean.df)), rep('unscaled',length(mean.df)))
  tempdf1 = data.frame(term =rep(term.names.clean,2), 
                     coeff.type.labels,  
                     rep(data.type.labels, 2), 
                     species, 
                     c(zmean.df,mean.df), 
                     c(zmedian.df,median.df), 
                     both.ci, 
                     row.names = 1:(length(mean.df)*2))
  #'estimate' is the mean
  colnames(tempdf1)=c('term','coeff.type', 'data.type', 'species','estimate','median','conf.low','conf.high')
  tempdf2=tempdf1 %>% mutate(overlap0=case_when(conf.low<0&conf.high>0~'yes',TRUE~'no'),
                             model='Density')
  
  #determine function outputs based on inputs
  if(coeff.type=='both'){return(tempdf2)} else {
    if(coeff.type=='scaled'){
      return(tempdf2 %>% filter(coeff.type=="scaled"))} else {
        return(tempdf2 %>% filter(coeff.type=='unscaled'))
      }
  }
}

## Plotting function for forest plots

plot.dens.fp <- function(x, foc.species='PP', center.val='estimate', title=foc.species){
  #development
  #x=fp.df(coeff.type = 'unscaled')
  
  plot.data = x %>% filter(species==foc.species) %>%
    mutate(row=row_number())
  
  fp <- ggplot(data=plot.data,
               aes(x = reorder(term, -row),y = plot.data[,center.val], ymin = conf.low, ymax = conf.high))+
    geom_pointrange(size=.4)+
    geom_hline(yintercept =0, linetype=2)+
    #xlab('Variable')+ ylab("Coefficient (mean and 95% CI)")+
    geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=term),width=0.5,cex=1)+ 
    #facet_wrap(~species,strip.position="left",nrow=9,scales = "free_y") +
    ggtitle(title) +
    theme_cowplot()+
    theme_minimal_hgrid() +
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=7),
          axis.title=element_blank(),
          #strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"), #facet labels
          legend.position = "none")+
    coord_flip()
  return(fp)
}

### ============================================================================
### Occupancy model
### ============================================================================

# generate data.frame of model results
occ.fp.df <- function(
    model=occ.mod1.list,
    longdata = occ.model.data, #long-format data before processing into jagged array
    group.array = y.list,
    jags.data = occ.jags.data.list,
    betavars.names = occ.betavars,
    species='BO',
    my.ci=c(.1,.9),
    coeff.type='unscaled') {
  
  
  
  mod=model[[species]]
  #data=mydata[[species]]
  df=longdata[[species]]
  mysims=mod$sims.list
  mymeans=mod$mean
  #values for unscaling
  
  # col.sd = apply(df[,col.vars], 2, function(x) sd(x, na.rm=T)) #SDs of continuous variables
  # col.mean=apply(df[,col.vars], 2, function(x) mean(x, na.rm=T))
  # per.sd = apply(df[,per.vars], 2, function(x) sd(x, na.rm=T)) #SDs of continuous variables
  # per.mean=apply(df[,per.vars], 2, function(x) mean(x, na.rm=T))
  beta.mean = apply(df[, betavars.names], 2, function(x) mean(x, na.rm=T)) 
  beta.sd = apply(df[, betavars.names], 2, function(x) sd(x, na.rm=T)) 
  time.mean=mean(df$timesince, na.rm=T)
  time.sd=sd(df$timesince, na.rm=T)
  
  #scaled parameters
  zint.c = mysims$intercept.c.star
  zint.e = mysims$intercept.e.star
  zbeta.coeffs.c = mysims$beta.c
  ztrt.coeffs.c=mysims$trteff.c.star
  ztime.coeffs.c=mysims$timesinceeff.c
  zbeta.coeffs.e = mysims$beta.e
  ztrt.coeffs.e=mysims$trteff.e.star
  ztime.coeffs.e=mysims$timesinceeff.e
  
  #terms to unscale categorical parameters
  #beta.unscale.c=sum(zbeta.coeffs.c/cont.sd + cont.mean)
  #time.unscale.c=ztime.coeffs.c/time.sd + time.mean
  
  #unscaled parameters
  beta.coeffs.c = mysims$beta.c %*% diag(1/beta.sd)
  trt.coeffs.c = ztrt.coeffs.c # - beta.unscale - time.unscale
  time.coeffs.c=mysims$timesinceeff.c/time.sd
  beta.coeffs.e = mysims$beta.e %*% diag(1/beta.sd)
  trt.coeffs.e = ztrt.coeffs.e # - beta.unscale - time.unscale
  time.coeffs.e=mysims$timesinceeff.e/time.sd
  
  term.names.clean = c('Intercept', 
                       'Control', 'Fire Only','Mech Only','Mech+Fire',
                       'Time[Control]', 'Time since Fire Only', 'Time since Mech Only', 'Time since Mech+Fire',
                       betavars.names)
  
  colnames(beta.coeffs.c) = betavars.names
  colnames(zbeta.coeffs.c) = betavars.names
  colnames(beta.coeffs.e) = betavars.names
  colnames(zbeta.coeffs.e) = betavars.names
  #colnames(trt.coeffs.c) = c('Control','Fire Only','Mechanical Only', 'Mechanical+Fire')
  #colnames(trt.coeffs.e) = c('Control','Fire Only','Mechanical Only', 'Mechanical+Fire')
  # colnames(ztrt.coeffs.c) = c('Intercept', 'Control','Fire Only','Mech Only', 'Mech+Fire')
  # colnames(ztrt.coeffs.e) = c('Intercept','Control','Fire Only','Mech Only', 'Mech+Fire')
  # colnames(trt.coeffs.c) = c('Intercept','Control','Fire Only','Mech Only', 'Mech+Fire')
  # colnames(trt.coeffs.e) = c('Intercept','Control','Fire Only','Mech Only', 'Mech+Fire')
  # colnames(time.coeffs.c)=paste0("Time:",c('Control','Fire Only','Mech Only', 'Mech+Fire'))
  # colnames(time.coeffs.e)=paste0("Time:",c('Control','Fire Only','Mech Only', 'Mech+Fire'))
  # colnames(ztime.coeffs.c)=paste0("Time:",c('Control','Fire Only','Mech Only', 'Mech+Fire'))
  # colnames(ztime.coeffs.e)=paste0("Time:",c('Control','Fire Only','Mech Only', 'Mech+Fire'))
  #### stoppping here - need to add .c to coeffs
  ztemp.df.c =data.frame(zint.c, ztrt.coeffs.c, ztime.coeffs.c, zbeta.coeffs.c)
  temp.df.c =data.frame(zint.c, trt.coeffs.c, time.coeffs.c, beta.coeffs.c)
  ztemp.df.e =data.frame(zint.e, ztrt.coeffs.e, ztime.coeffs.e, zbeta.coeffs.e)
  temp.df.e =data.frame(zint.e, trt.coeffs.c, time.coeffs.e, beta.coeffs.e)
  ### still need .e's  ztemp.df = rbind(ztemp.df.c, ztemp.df.e)
  #temp.df = rbind(temp.df.c, temp.df.e)
  #ztemp.df = rbind(ztemp.df.c, ztemp.df.e)
  
  ## take summary stats of the dataframes
  zmean.df.c = apply(ztemp.df.c,2,mean)
  mean.df.c = apply(temp.df.c,2,mean)
  zmean.df.e = apply(ztemp.df.e,2,mean)
  mean.df.e = apply(temp.df.e,2,mean)
  zmedian.df.c=apply(ztemp.df.c,2,median)
  median.df.c=apply(temp.df.c,2,median)
  zmedian.df.e=apply(ztemp.df.e,2,median)
  median.df.e=apply(temp.df.e,2,median)
  zci.df.c = t(apply(ztemp.df.c, 2, function(x) quantile(x, my.ci)))
  ci.df.c = t(apply(temp.df.c, 2, function(x) quantile(x, my.ci)))
  zci.df.e = t(apply(ztemp.df.e, 2, function(x) quantile(x, my.ci)))
  ci.df.e = t(apply(temp.df.e, 2, function(x) quantile(x, my.ci)))
  all.ci=rbind(zci.df.c,zci.df.e,ci.df.c,ci.df.e)
  #combine scaled and unscaled coeffs
  #all.coeffs=rbind(ztemp.df, temp.df)
  nrow(all.ci)
  coeff.type.labels=c(rep('scaled', length(zmean.df.c)+length(zmean.df.e)), rep('unscaled',length(mean.df.c)+length(mean.df.e)))
  zdata.type.labels.c=c('Intercept', rep('categorical',4), rep('time',4), rep('continuous',length(betavars.names)))
  data.type.labels.c=c('Intercept',rep('categorical',4),rep('time',4), rep('continuous',length(betavars.names)))
  zdata.type.labels.e=c('Intercept',rep('categorical',4), rep('time',4), rep('continuous',length(betavars.names)))
  data.type.labels.e=c('Intercept',rep('categorical',4),rep('time',4), rep('continuous',length(betavars.names)))
  zmod.type.labels = c(rep('colonization',length(zdata.type.labels.c)), rep('persistence',length(zdata.type.labels.e)))
  mod.type.labels = c(rep('colonization',length(data.type.labels.c)), 
                      rep('persistence',length(data.type.labels.e)))
  ## this is still wrong
  fp.df = data.frame(parameter =rep(term.names.clean, 4),#c(names(zmean.df.c), names(zmean.df.e), names(mean.df.c), names(mean.df.e)), 
                     coeff.type.labels,  
                     c(zdata.type.labels.c,zdata.type.labels.e,data.type.labels.c,data.type.labels.e), 
                     c(zmod.type.labels,mod.type.labels),
                     species, c(zmean.df.c,zmean.df.e, mean.df.c,mean.df.e), 
                     c(zmedian.df.c,zmedian.df.e,median.df.c,median.df.e),
                     all.ci,
                     row.names = 1:length(coeff.type.labels))
  colnames(fp.df)=c('term','coeff.type', 'data.type','model', 'species','estimate','median','conf.low','conf.high')
  fp.df=fp.df %>% mutate(overlap0=case_when(conf.low<0&conf.high>0~'yes',TRUE~'no'))
  
  if(coeff.type=='both'){return(fp.df)} else {
    if(coeff.type=='scaled'){
      return(fp.df %>% filter(coeff.type=="scaled"))} else {
        return(fp.df %>% filter(coeff.type=='unscaled'))
      }
  }
}
#end function

#df=occ.fp.df(coeff.type = 'both')

#generate plot of model results
plot.occ.fp <- function(df, 
                       myspecies='BO', 
                       title=unique(df$species)){
  funpacks = c('ggplot2','cowplot')
  lapply(funpacks, require, character.only=T)
  df$row = 1:nrow(df)
  mymod.df=df[which(df$species==myspecies),]
  mymod.df$term = factor(mymod.df$term, levels = rev(unique(mymod.df$term)))
  fp <- ggplot(data=mymod.df,
               aes(x = term,
                   y = estimate, 
                   ymin = conf.low, 
                   ymax = conf.high))+
    #geom_pointrange(size=.4)+
    #these are the nuts and bolts
    geom_errorbar(aes(ymin=conf.low, ymax=conf.high,col=term,linetype=model, group=term, alpha=coeff.type),width=0.6,cex=1,position=position_dodge2(width = 1))+ 
    geom_point(aes(col=term, alpha=coeff.type), position = position_dodge2(width=.6))+
    #specify colors, linetypes, etc.
    scale_linetype_manual(values = c('solid','dotted'), 
                          labels=c('colonization','persistence'),
                          name="process")+
    #scale_color_discrete(guide='none')+
    scale_alpha_manual(values = c(1, .3), guide='none')+
    geom_hline(yintercept =0, linetype=2)+
    xlab('Variable')+ ylab("Coefficient (mean and 90% CI)")+
    #facet_wrap(~species,strip.position="left",nrow=9,scales = "free_y") +
    ggtitle(title) +
    theme_cowplot()+
    theme_minimal_hgrid() +
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_text(size=8),
          axis.text.x=element_text(size=7),
          axis.title=element_text(size=12,face="bold"),
          #strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"), #facet labels
          legend.position = "top")+
    coord_flip()
  return(fp)
}

##------------------------------------------------------------------
## logistic curves for select variables ####
##------------------------------------------------------------------
# for plotting univariate relationships with pcol or pext
#function inputs
curves.fun = function(model=occ.mod6.list,
                      datalist=occ.jags.data.list,
                      df = occ.model.data,
                      species='IC',
                      myvar='species.ba',
                      ndraws=30,
                      n.xvalues=100)
  #end inputS
{
  
  #development
  # model=occ.mod6.list
  # datalist=occ.jags.data.list
  # df = occ.model.data
  # species='PP'
  # myvar='timesince'
  # ndraws=30
  # n.xvalues=100
  
  #Need these for later
  logistic <- function(x) 1/(1+exp(-x))
  inv.logit <- function(x) exp(x)/(1+exp(x))
  
  #subset input data
  mod=model[[species]]
  mydraws = sample(mod$mcmc.info$n.samples, ndraws)
  mydf=df[[species]]
  data=datalist[[species]]
  
  #housekeeping
  nplots=data$nplots
  nsamples=mod$mcmc.info$n.samples
  # beta.means.c=apply(mydf[,col.vars],2, function(x) mean(x, na.rm=T))
  # beta.sd.c=apply(mydf[,col.vars],2, function(x) sd(x, na.rm=T))
  # beta.means.e=apply(mydf[,per.vars],2, function(x) mean(x, na.rm=T))
  # beta.sd.e=apply(mydf[,per.vars],2, function(x) sd(x, na.rm=T))
  beta.mean = apply(mydf[, occ.betavars], 2, function(x) mean(x, na.rm=T)) 
  beta.sd = apply(mydf[, occ.betavars], 2, function(x) sd(x, na.rm=T)) 
  time.mean=mean(mydf$timesince, na.rm=T)
  time.sd=sd(mydf$timesince, na.rm=T)

  #unscale coefficient samples
  beta.unscale.c= apply(mod$sims.list$beta.c * beta.mean/beta.sd,1,sum) #+ mod$sims.list$timesinceeff.c*time.mean/time.sd#Bx*mean(x)/sd(x) then sum it all up for each model iteration
  beta.unscale.e= apply(mod$sims.list$beta.e * beta.mean/beta.sd,1,sum) #+ mod$sims.list$timesinceeff.e*time.mean/time.sd#Bx*mean(x)/sd(x) then sum it all up for each model iteration
  
  #full x values used in model
  x.var = unlist(mydf[,myvar])
  #zx.var = scale(x.var)
  
  #x values to be used for plotting
  x.values =  seq(min(x.var,na.rm=T),max(x.var,na.rm=T),length.out=n.xvalues) #zx.values*sd(x.var, na.rm = T) + mean(x.var, na.rm=T)
  zx.values = (x.values-mean(x.var, na.rm=T))/sd(x.var, na.rm=T) #seq(min(zx.var,na.rm=T),max(zx.var,na.rm=T),length.out=n.xvalues)
  
  #put them in a df
  x.df = data.frame(Var1=1:n.xvalues, zx.values=zx.values, x.values = x.values)
  
  
  if(myvar=='timesince'){
    ztimesince=zx.values
    timesince=x.values
    other.values = beta.mean
    
    zvar.coeff.c= rep(0,mod$mcmc.info$n.samples) #don't need this, since myvar==timesince
    var.coeff.c= rep(0,mod$mcmc.info$n.samples) #don't need this, since myvar==timesince
    ztimesinceeff.c = mod$sims.list$timesinceeff.c#[mydraws,]
    timesinceeff.c = ztimesinceeff.c/time.sd
    zother.coeffs.c=mod$sims.list$beta.c#[mydraws,]
    other.coeffs.c = zother.coeffs.c/beta.sd
    
    zvar.coeff.e= rep(0,mod$mcmc.info$n.samples)
    var.coeff.e= rep(0,mod$mcmc.info$n.samples)
    ztimesinceeff.e = mod$sims.list$timesinceeff.e#[mydraws,]
    timesinceeff.e = ztimesinceeff.e/time.sd
    zother.coeffs.e=mod$sims.list$beta.e#[mydraws,]
    other.coeffs.e = zother.coeffs.e/beta.sd
    
    ztrteff.c =   mod$sims.list$trteff.c.star#[mydraws,]
    trteff.c = ztrteff.c - timesinceeff.c*time.mean/time.sd
    ztrteff.e =   mod$sims.list$trteff.e.star#[mydraws,]
    trteff.e = ztrteff.e - timesinceeff.e*time.mean/time.sd
    
  } else {
    coeff.ind = which(occ.betavars==myvar)
    
    timesince = time.mean #if myvar != timesince, this should be zero #mean(mydf$timesince, na.rm = T)
    ztimesince = 0
    other.values = beta.mean[-coeff.ind]
    
    ztimesinceeff.c = mod$sims.list$timesinceeff.c
    ztimesinceeff.e = mod$sims.list$timesinceeff.e
    timesinceeff.c = mod$sims.list$timesinceeff.c/time.sd #[mydraws,]#/sd(mydf$timesince, na.rm=T)
    timesinceeff.e = mod$sims.list$timesinceeff.e/time.sd
    
    other.vars.mean = apply(mydf[,occ.betavars[-coeff.ind]],2,function(x)mean(x,na.rm=T))
    other.vars.sd=apply(mydf[,occ.betavars[-coeff.ind]],2,function(x)sd(x,na.rm=T))
    
    zother.coeffs.c=mod$sims.list$beta.c[,-coeff.ind]
    other.coeffs.c = zother.coeffs.c/other.vars.sd
    zvar.coeff.c = mod$sims.list$beta.c[,coeff.ind]
    var.coeff.c = zvar.coeff.c/beta.sd[coeff.ind]
    zother.coeffs.e=mod$sims.list$beta.e[,-coeff.ind]
    other.coeffs.e = zother.coeffs.e/other.vars.sd
    zvar.coeff.e = mod$sims.list$beta.e[,coeff.ind]
    var.coeff.e = zvar.coeff.e/beta.sd[coeff.ind]
    
    # ztrteff.c =   mod$sims.list$trteff.c.star
    # ztrteff.e =   mod$sims.list$trteff.e.star
    # trteff.c = ztrteff.c - timesinceeff.c*time.mean/time.sd
    # trteff.e = ztrteff.e - timesinceeff.e*time.mean/time.sd
    
    trteff.c=matrix(0,nrow=nsamples, ncol=4)
    trteff.e=matrix(0,nrow=nsamples, ncol=4)
    ztrteff.c=matrix(0,nrow=nsamples, ncol=4)
    ztrteff.e=matrix(0,nrow=nsamples, ncol=4)
    
  }
  #categorical terms - only need intercept (unless myvar==timesince, but trt.effs are included above)
  # plot.type.c= mod$sims.list$plot.type.c.star#[mydraws,]
  # trteff.c =   mod$sims.list$trteff.c.star#[mydraws,]
  zintercept.c= mod$sims.list$intercept.c.star#[mydraws]# - beta.unscale.c
  intercept.c = zintercept.c - beta.unscale.c
  # standeff.c = mod$sims.list$standeff.c.star#[mydraws,,]
  # locationeff.c=mod$sims.list$locationeff.c.star#[mydraws,,,]
  # plot.type.e=  mod$sims.list$plot.type.e.star#[mydraws,]
  # trteff.e =    mod$sims.list$trteff.e.star#[mydraws,]
  zintercept.e=  mod$sims.list$intercept.e.star#[mydraws]# - beta.unscale.c
  intercept.e = zintercept.e - beta.unscale.e
  # standeff.e =  mod$sims.list$standeff.e.star#[mydraws,,]
  # locationeff.e=mod$sims.list$locationeff.e.star#[mydraws,,,]
  
  y.c=array(NA, dim=c(n.xvalues, ndraws, data$ntrt))  #rows=f(x), cols=draws, slices=trts
  y.e=array(NA, dim=c(n.xvalues, ndraws, data$ntrt)) #rows=f(x), cols=draws, slices=trts
  zy.c=array(NA, dim=c(n.xvalues, ndraws, data$ntrt))
  zy.e=array(NA, dim=c(n.xvalues, ndraws, data$ntrt))
  dimnames(y.c)=list(c(1:n.xvalues),c(1:ndraws),c('trt1','trt2','trt3','trt4'))
  dimnames(y.e)=list(c(1:n.xvalues),c(1:ndraws),c('trt1','trt2','trt3','trt4'))
  dimnames(zy.c)=list(c(1:n.xvalues),c(1:ndraws),c('trt1','trt2','trt3','trt4'))
  dimnames(zy.e)=list(c(1:n.xvalues),c(1:ndraws),c('trt1','trt2','trt3','trt4'))
  
  #trteff.c.mat=matrix(rep(trteff.c,4), ndraws, 4)
  #trteff.e.mat=matrix(rep(trteff.e,4), ndraws, 4)
  for(i in 1:ndraws){
    for(t in 1:4){
      for(s in 1:3){
        for(p in 1:nplots[s,t]){
          # y.c[,i,t] = logistic(intercept.c[mydraws[i]] + 
          #                        # plot.type.c[mydraws[i]]+
          #                        trteff.c[mydraws[i],t]  +
          #                        # standeff.c[mydraws[i],s,t] +
          #                        # locationeff.c[mydraws[i],s,t,p] +
          #                        timesinceeff.c[mydraws[i],t]*timesince + #timesince = 0 if myvar!=timesince
          #                        sum(other.coeffs.c[mydraws[i],] *  other.values)+
          #                        var.coeff.c[mydraws[i]]*x.values)
          # 
          # y.e[,i,t] = logistic(intercept.e[mydraws[i]] +
          #                        # plot.type.e[mydraws[i]]+
          #                        trteff.e[mydraws[i],t]  +
          #                        #standeff.e[mydraws[i],s,t] +
          #                        #locationeff.e[mydraws[i],s,t,p] +
          #                        timesinceeff.e[mydraws[i],t]*ztimesince +
          #                        sum(other.coeffs.e[mydraws[i],] * other.values)+
          #                        var.coeff.e[mydraws[i]]*x.values)
          zy.c[,i,t] = logistic(zintercept.c[mydraws[i]] + 
                                 # plot.type.c[mydraws[i]]+
                                 ztrteff.c[mydraws[i],t]  +
                                 # standeff.c[mydraws[i],s,t] +
                                 # locationeff.c[mydraws[i],s,t,p] +
                                 ztimesinceeff.c[mydraws[i],t]*ztimesince + #timesince = 0 if myvar!=timesince
                                 sum(zother.coeffs.c[mydraws[i],] *  0)+
                                 zvar.coeff.c[mydraws[i]]*zx.values)
          
          zy.e[,i,t] = logistic(zintercept.e[mydraws[i]] +
                                 # plot.type.e[mydraws[i]]+
                                 ztrteff.e[mydraws[i],t]  +
                                 #standeff.e[mydraws[i],s,t] +
                                 #locationeff.e[mydraws[i],s,t,p] +
                                 ztimesinceeff.e[mydraws[i],t]*ztimesince +
                                 sum(zother.coeffs.e[mydraws[i],] * 0)+
                                 zvar.coeff.e[mydraws[i]]*zx.values)
          
          
        }
      }
    }
  }
  #y.e.mean=apply(y.e, c(1,3), mean)
  
  y.c.median=matrix(NA, n.xvalues, 4);dimnames(y.c.median)=list(c(1:n.xvalues), c('trt1','trt2','trt3','trt4'))
  y.e.median=matrix(NA, n.xvalues, 4);dimnames(y.e.median)=list(c(1:n.xvalues), c('trt1','trt2','trt3','trt4'))
  zy.c.median=matrix(NA, n.xvalues, 4);dimnames(zy.c.median)=list(c(1:n.xvalues), c('trt1','trt2','trt3','trt4'))
  zy.e.median=matrix(NA, n.xvalues, 4);dimnames(zy.e.median)=list(c(1:n.xvalues), c('trt1','trt2','trt3','trt4'))
  for(t in 1:4){
    # y.c.median[,t]=logistic(median(intercept.c)+
    #                           #mean(plot.type.c) +
    #                           apply(trteff.c, 2, median)[t] +
    #                           #median(standeff.c) +
    #                           #median(locationeff.c,na.rm=T)[t] +
    #                           apply(timesinceeff.c, 2, median)[t]*timesince +
    #                           sum(apply(other.coeffs.c, 2, median)*other.values)+
    #                           median(var.coeff.c)*x.values)
    # y.e.median[,t]=logistic(median(intercept.e)+
    #                           #median(plot.type.e) +
    #                           apply(trteff.e, 2, median)[t] +
    #                           #median(standeff.e) +
    #                           #median(locationeff.e,na.rm=T)[t] +
    #                           apply(timesinceeff.e, 2, median)[t]*timesince +
    #                           sum(apply(other.coeffs.e, 2, median)*other.values) +
    #                           median(var.coeff.e)*x.values)
    zy.c.median[,t]=logistic(median(zintercept.c)+
                              #mean(plot.type.c) +
                              apply(ztrteff.c, 2, median)[t] +
                              #median(standeff.c) +
                              #median(locationeff.c,na.rm=T)[t] +
                              apply(ztimesinceeff.c, 2, median)[t]*ztimesince +
                              sum(apply(zother.coeffs.c, 2, median)*0)+
                              median(zvar.coeff.c)*zx.values)
    zy.e.median[,t]=logistic(median(zintercept.e)+
                              #median(plot.type.e) +
                              apply(ztrteff.e, 2, median)[t] +
                              #median(standeff.e) +
                              #median(locationeff.e,na.rm=T)[t] +
                              apply(ztimesinceeff.e, 2, median)[t]*ztimesince +
                              sum(apply(other.coeffs.e, 2, median)*0) +
                              median(zvar.coeff.e)*zx.values)
    
  }
  # plot(x.var, mydf$occ)
  # for(i in 1:ndraws){
  #   for(t in 1:4)
  #     lines(x.values,y.c[,i,t],col=t)
  # }
  # for(t in 1:4){
  #   lines(x.values, y.c.median[,t], lwd=3, col=5)
  # }
  
  #using scaled coefficients and data
  y.c = zy.c
  y.e = zy.e
  y.c.median = zy.c.median
  y.e.median = zy.e.median
  
  y.c.long = suppressWarnings(melt(y.c, varnames = c('Var1', 'draw', 'trt')) %>% 
    left_join(.,x.df, by='Var1') %>%
    mutate(species=species,
           myvar=myvar,
           model='colonization')) %>%
    filter(x.values>0)
  y.e.long = suppressWarnings(melt(y.e, varnames = c('Var1', 'draw', 'trt')) %>% 
    left_join(.,x.df, by='Var1') %>%
    mutate(species=species,
           myvar=myvar,
           model='persistence')) %>%
    filter(x.values>0)
  colnames(y.c.long) = c('n.x','n.draw','trt','y','zx.values','x.values','species','myvar','model')
  colnames(y.e.long) = c('n.x','n.draw','trt','y','zx.values','x.values','species','myvar','model')
  y.c.median.long=as.data.frame(y.c.median) %>%
    mutate(x.values=x.values,
           zx.values=zx.values) %>%
    pivot_longer(trt1:trt4, names_to = 'trt',values_to = 'y.median') %>%
    mutate(species=species,
           myvar=myvar,
           model='colonization') %>%
    filter(x.values>0)
  y.e.median.long=as.data.frame(y.e.median) %>%
    mutate(x.values=x.values,
           zx.values=zx.values) %>%
    pivot_longer(trt1:trt4, names_to = 'trt',values_to = 'y.median') %>%
    mutate(species=species,
           myvar=myvar,
           model='persistence') %>%
    filter(x.values>0)
  
  # #trying something else for medians
  # y.c.median.long = y.c.long %>%
  #   group_by(n.x, zx.values, x.values, trt, species, myvar, model) %>%
  #   summarise(y.median = median(y))
  # y.e.median.long = y.e.long %>%
  #   group_by(n.x, zx.values, x.values, trt, species, myvar, model) %>%
  #   summarise(y.median = median(y))
  
  return.list = list(y.c=y.c.long,
                     y.e=y.e.long,
                     c.median=y.c.median.long,
                     e.median = y.e.median.long)
  
  
  
  return(return.list)
}

#TROUBLESHOOTING
#inputs for curveplotfun
# my.species = species
# var=myvar
# listdf=occ.model.data
# #dfdata=datadf
# alldf = y.e.long
# med.df=y.e.median.long
# xlab="var"
# #end inputs
# ################ start function ################
# #establish set colors for treatments
# treatcolors = cal_palettes$superbloom3[c(7,1,5,2)]
# 
# #establish set colors for species 
# spp.colors = cal_palettes$kelp1[c(3,5,2,1,4,6)]
# 
# #change figure labeling based on variable
# if(var!='timesince') {
#   trt.filter <- 'trt1'
#   trt.labels <- my.species
#   spp.filter <- myspecies
#   mycolors <- spp.colors[which(myspecies==my.species)]} else {
#   #var==timesince  
#   trt.filter <- c('trt1','trt2','trt3','trt4')
#   trt.labels <- c('Control','Fire Only','Mech Only','Mech+Fire')
#   spp.filter <- my.species
#   mycolors <- treatcolors
#  }
# 
# #filter input data
#   datadf = listdf[[spp.filter]]
#   #var.vals=dfdata[,var]
#   alldf=alldf %>% filter(species%in%spp.filter, myvar==var, trt%in%trt.filter)
#   meddf=med.df %>% filter(species%in%spp.filter, myvar==var, trt%in%trt.filter)
#   #meddf$trt = factor(meddf$trt, labels=c('Control','Fire.Only','Mechanical.Only','Mechanical.Fire'))
#   
#   #now plot
#   ggplot(data=listdf[[my.species]], aes_string(x=var,y='occupancy'))+
#     #geom_point(position = position_jitter(height =.2,width = 1), size=.1) +
#     scale_color_manual(values = mycolors, labels=trt.labels) +
#     labs(color='Treatment', y='probability', x=xlab)+
#     #individual mcmc samples
#     geom_line(data=alldf, aes(x=zx.values, 
#                               y=y,
#                               group=interaction(n.draw,trt),
#                               color=trt), alpha=.2) +
#     #medians
#     geom_line(data=meddf, aes(x=zx.values, y=y.median, group=trt, color=trt), lwd=2) +
#     theme_cowplot()+
#     facet_wrap(~model)+
#     ggtitle(paste(my.species, var, sep = " ")) +
#     geom_hline(yintercept = .5) +
#     theme(legend.position = 'bottom', 
#           legend.margin = margin(t=-.5, unit = 'cm'),
#           legend.title = element_blank(),
#           legend.text = element_text(size=10),
#           plot.margin = unit(c(0,0.01,-.1,-.5),'cm'),
#           plot.title = element_text(size=12, margin =margin(.05, 0, .05, 0, "cm")) ,
#           strip.text.x = element_text(size=10, margin = margin(.05, 0, .05, 0, "cm")), 
#           axis.text.x = element_text(size=9,vjust = 3),
#           axis.text.y=element_text(size=10))





## ==========================================================
## Load model objects, generate model results, and export functions and model results ####
## ==========================================================

dens.mod.results = bind_rows(lapply(myspecies, function(x){dens.fp.df(model = dens.mod6.output,
                                                   species=x, 
                                                   coeff.type = 'both')}))
occ.mod.results = bind_rows(lapply(myspecies, function(x){occ.fp.df(model = occ.mod6.list,
                                                                 species=x, 
                                                                 coeff.type = 'both')}))

#univariate relationships - individual runs
occ.curves.y = bind_rows(lapply(seq_along(myspecies), function(x) {
  lapply(seq_along(c(occ.betavars, 'timesince')), function(j) {
    output.list = curves.fun(model = occ.mod6.list,
               datalist = occ.jags.data.list,
               df = occ.model.data, 
               species = myspecies[x],
               myvar = c(occ.betavars, 'timesince')[j],
               ndraws=50, 
               n.xvalues = 100)
    output.df = bind_rows(output.list$y.c, output.list$y.e)
    return(output.df)
  })
})) %>%
  #align naming conventions with other products
  mutate(treatment = factor(trt, levels=c('trt1','trt2','trt3','trt4'), 
                            labels = c('Control','Fire Only','Mech Only','Mech+Fire')),
         term = case_when(myvar=='trtnum'~'Treatment applications',
                           myvar=='species.ba'~'Conspecific BA',
                           myvar=='total.ba'~'Total BA',
                           myvar=='timesince'&treatment=='Control' ~ 'Time[Control]',
                           myvar=='timesince'&treatment=='Fire Only' ~ 'Time since Fire Only',
                           myvar=='timesince'&treatment=='Mech Only' ~ 'Time since Mech Only',
                           myvar=='timesince'&treatment=='Mech+Fire' ~ 'Time since Mech+Fire',
                           TRUE~myvar))



#univariate relationships - medians
occ.curves.median = bind_rows(lapply(seq_along(myspecies), function(x) {
  lapply(seq_along(c(occ.betavars, 'timesince')), function(j) {
    output.list = curves.fun(model = occ.mod6.list,
                             datalist = occ.jags.data.list,
                             df = occ.model.data, 
                             species = myspecies[x],
                             myvar = c(occ.betavars, 'timesince')[j],
                             ndraws=50, 
                             n.xvalues = 100)
    output.df = bind_rows(output.list$c.median, output.list$e.median)
  })
})) %>%
  mutate(treatment = factor(trt, levels=c('trt1','trt2','trt3','trt4'), 
                            labels = c('Control','Fire Only','Mech Only','Mech+Fire')),
         #align naming convention with other products
  term = case_when(myvar=='trtnum'~'Treatment applications',
                    myvar=='species.ba'~'Conspecific BA',
                    myvar=='total.ba'~'Total BA',
                    myvar=='timesince'&treatment=='Control' ~ 'Time[Control]',
                    myvar=='timesince'&treatment=='Fire Only' ~ 'Time since Fire Only',
                    myvar=='timesince'&treatment=='Mech Only' ~ 'Time since Mech Only',
                    myvar=='timesince'&treatment=='Mech+Fire' ~ 'Time since Mech+Fire',
                    TRUE~myvar))


#save functions and data into .rdata file
save(dens.fp.df, 
     plot.dens.fp, 
     occ.fp.df,
     plot.occ.fp,
     dens.mod.results,
     occ.mod.results,
     occ.curves.y,
     occ.curves.median,
     file = "Cresults_and_functions.rdata")

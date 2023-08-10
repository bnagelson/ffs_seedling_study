
###
### FFS Seedling Study
### Abundance model - data preparation and model fitting


packs <- c('dplyr','tidyr', 'MCMCvis', 'mcmcplots','bayesplot', 'forcats', 'tictoc','jagshelper')
lapply(packs, require, character.only=T)
#current.comp = "E:/Other computers/My PC"
#my.wd=file.path(current.comp, "UNR/FFS Seedling Study")

#sometimes JagsUI::jags produces an error that it cannot find the jags installation
# this seems to be the fix: 
# Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.0")

### prep data ####

#load data
load("data_processing_foster.rdata",
     verbose=T)
load("cool_warm_precip_v2.rdata",
     verbose=T) 

myspecies = c('BO','DF','IC','PP','SP','WF')

#climate data
clim.wide = williams2021$abund.data 

#join
all = left_join(abund.all, clim.wide,
                by=c('PlotID','year')) 

#prep a clean data frame for each species - only necessary predictors, response data, and grouping vars
betavars = c('trtnum','species.ba','total.ba',
               'cool.ppt','warm.ppt','warm.cwd') #predictors (except for timesince)
grouping = c('comp','treatment', 'trt','plotsize') #grouping variables

model.data = lapply(seq_along(myspecies), function(x) {
  all %>% 
    filter(species==myspecies[x],
           timestep!='pre_treatment') %>%
    select(grouping, betavars, timesince, seedlings.count) %>%
    arrange(trt, comp) %>%
    drop_na(c(betavars, grouping))
})
names(model.data) = myspecies

#construct model matrices for continuous predictors (timesince comes later)
model.mats = lapply(seq_along(myspecies), function(x) {
  model.matrix(~ trtnum +
                 species.ba +
                 total.ba +
                 cool.ppt +
                 warm.ppt + 
                 warm.cwd, #mean annual precip for the three years prior to sampling, with annual precip calculated as june-may
               data = model.data[[x]])[,-1] #remove intercept column
})
names(model.mats) = myspecies

#timseince 
timesince.list = lapply(seq_along(myspecies), function(x) {
  model.data[[x]] %>%
    select(timesince)
})
names(timesince.list) = myspecies

#package data together for JAGS
jags.data.list = lapply(seq_along(myspecies), function(x) {
  df = model.data[[x]] %>%
  #manipulate compartment labeling into a nested factor (1:3 within each treatment)
    mutate(standfac = factor(comp, levels = c(40,240,590,60,340,400,190,350,490,180,380,570),
                             labels = rep(1:3,4)))
  list(N = nrow(df),
       observed = df$seedlings.count,
       plotsize = df$plotsize,
       stand=as.numeric(df$standfac),
       nstands = 3,
       trt=df$trt,
       ntrt=4,
       nbeta = length(betavars),
       z.model.matrix = apply(model.mats[[x]],2,scale),
       ztimesince = as.numeric(scale(timesince.list[[x]])))
})
names(jags.data.list) = myspecies

### JAGS ###################################################

#define parameters to monitor
params =  c('intercept.star','trt.eff.star',
             'timesince.eff', 
             'beta',
             'stand.eff.star',
            'predicted',
            'mu') #monitoring 'spread' slows down the model considerably

### Run models
dens.mod6.output = lapply(seq_along(myspecies), function(x) {
  print(paste0('Working on the ', myspecies[x], ' model'))
  #tic()
  mod = jagsUI::jags(data = jags.data.list[[x]], 
             model.file = "dens.mod2.txt",
             parameters.to.save = params,
             #inits = list(initsfun(),initsfun()),
             n.iter=15000,
             n.adapt=1500,
             n.thin = 10,
             n.chains=4,
             parallel = T)
  #toc()
  #beepr::beep()
  return(mod)
  #42min with 15,000 iters, 1500 adapt, 2 chains, parallel=F
  #20min with 10,000 iters, 1500 adapt, 3 chains, parallel=T (good mixing)
})
names(dens.mod6.output) = myspecies

# save(jags.data.list, 
#      model.mats,
#      model.data,
#      dens.mod6.output,
#      file = "dens.mod.output.rdata")





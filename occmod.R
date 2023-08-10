### FFS Seedling Study
### Occupacy modeling

#load packages and set working directory
packs <- c('dplyr','tidyr', 'MCMCvis', 'mcmcplots','bayesplot', 'forcats', 'tictoc')
lapply(packs, require, character.only=T)

#sometimes JagsUI::jags produces an error that it cannot find the jags installation
# this seems to be the fix: 
# Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.1")

### prep data ####

#load data
load("data_processing_foster.rdata",
     verbose=T)
load("cool_warm_precip_v2.rdata",
     verbose=T) #window.calc has 3- and 5-year means for winter, spring, and dry-season precip,
# wateryear cwd, and january min temps

myspecies <- c('BO','DF','IC','PP','SP','WF')


clim.wide = williams2021$occ.data

all = left_join(occ.joined, clim.wide, 
                by=c('ID', 'year'))

#prep a clean data frame for each species - only necessary predictors, response data, and grouping vars
occ.betavars = c('trtnum','species.ba','total.ba',
                 'cool.ppt','warm.ppt','warm.cwd')
grouping = c('ID','PlotID','comp','treatment', 'trt', 'type', 'timestep.num', 'year')

occ.model.data = lapply(seq_along(myspecies), function(x) {
  df = all %>% 
    replace_na(., list(trtnum=0)) %>%
    filter(species==myspecies[x],
           #timestep!='pre_treatment'
           ) %>%
    select(grouping, occ.betavars, timesince, occupancy) %>%
    arrange(ID, PlotID, trt, comp) %>%
    drop_na(c(occ.betavars, grouping)) %>%
    group_by(ID) %>%
    mutate(samp1 = min(timestep.num),
           samp2 = sort(timestep.num, decreasing = F)[2],
           samp3 = sort(timestep.num, decreasing = F)[3],
           samp4 = sort(timestep.num, decreasing = F)[4],
           samp5 = sort(timestep.num, decreasing = F)[5],
           samp6 = sort(timestep.num, decreasing = F)[6],
           samp7 = sort(timestep.num, decreasing = F)[7],
           samp8 = sort(timestep.num, decreasing = F)[8],
           samp9 = sort(timestep.num, decreasing = F)[9]) %>% ungroup() %>% 
    pivot_longer(cols = c(samp1:samp9), values_to = 'newtime1', 
                 names_to = 'sample', names_prefix = "samp") %>%
    #'sample' is the nth time the plot was sampled. 
    distinct(ID, year, sample, .keep_all = T) %>%
    filter(timestep.num==newtime1) %>%
    select(-newtime1) %>%
    group_by(ID) %>%
    mutate(nsamples=n()) %>%
    relocate(nsamples, .after=sample) %>% ungroup() %>%
    filter(nsamples>1)
  return(df)
})
names(occ.model.data) = myspecies

y.list = lapply(seq_along(myspecies), function(x){
  occ.model.data[[x]] %>%
  mutate(standfac = factor(comp, levels = c(40,240,590,60,340,400,190,350,490,180,380,570),
                           labels = rep(1:3,4)),
         typefac = factor(type, levels = c('CSP','FI'),
                          labels = 1:2)) %>%
    group_by(comp) %>% 
    mutate(location=as.integer(fct_inorder(PlotID))) %>% 
    ungroup() %>%
    select(ID, PlotID, location, comp, standfac, trt, type, typefac, nsamples, sample, occupancy) %>%
    #select(ID, sample, occupancy) %>%
    pivot_wider(#id_cols = c('ID', 'PlotID', 'comp', 'type', 'nsamples'), 
                names_from = sample, values_from = occupancy) #%>%
    #select(-ID)
})
names(y.list) = myspecies

##function to convert individual vars into into N x nSample matrices
matrix.indvars = function(df, value, remove.IDs=TRUE, scale=FALSE){
  value.scale = scale(df[,value])
  if(scale==TRUE){df[,value] = value.scale}
  temp = df %>% #filter(!is.na(sample)) %>%
    select(ID, PlotID, comp, trt, sample, nsamples, type, value)
  
  temp.wide = temp %>%  pivot_wider(id_cols = c(ID, PlotID, comp, trt,type, nsamples), 
                                     names_from = sample, 
                                     values_from = value) %>%
    ungroup()
  if(remove.IDs==TRUE){return(temp.wide %>% select(-c(ID, PlotID, comp, trt, type, nsamples)))} else 
    return(temp.wide)
}

#colonization occupancy array
betavar.array.list = lapply(seq_along(myspecies), function(x) {
  array(unlist(lapply(occ.betavars, 
                      FUN = function(j) matrix.indvars(occ.model.data[[x]], j, scale = T))),
        dim = c(nrow(y.list[[x]]),max(occ.model.data[[x]]$nsamples),length(occ.betavars)),
        dimnames=list(NULL,NULL,occ.betavars))
})
names(betavar.array.list) = myspecies
dim(betavar.array.list$BO)

# package data for JAGS
occ.jags.data.list = lapply(seq_along(myspecies), function(x) {
  # df = occ.model.data[[x]] %>%
  #   mutate(standfac = factor(comp, levels = c(40,240,590,60,340,400,190,350,490,180,380,570),
  #                            labels = rep(1:3,4)))
  nplots = y.list[[x]] %>%
    group_by(comp) %>%
    mutate(nplots=n_distinct(PlotID)) %>% 
    ungroup() %>%
    distinct(trt, standfac,nplots) %>%
    arrange(trt) %>%
    pivot_wider(names_from = trt, values_from = nplots)
  
  list(N=nrow(y.list[[x]]),
       y=as.matrix(select(y.list[[x]], `1`:`9`)),
       nSample=y.list[[x]]$nsamples,
       trt=y.list[[x]]$trt,
       ntrt=4,
       stand=as.numeric(y.list[[x]]$standfac),
       nstands=3,
       nplots=as.matrix(nplots[,2:5]),
       location=y.list[[x]]$location,
       type=as.numeric(y.list[[x]]$typefac),
       timesince=matrix.indvars(occ.model.data[[x]], 'timesince', scale = T),
       betavar.array = betavar.array.list[[x]],
       nbetavars = length(occ.betavars))
})
names(occ.jags.data.list)=myspecies

######################################

### Run JAGS

params = c('intercept.c.star','intercept.e.star',
           'plot.type.c.star','plot.type.e.star',
           #'locationeff.c.star','locationeff.e.star',
           "trteff.c.star", "trteff.e.star",
           "standeff.c.star", "standeff.e.star",
           "timesinceeff.c",
           "timesinceeff.e",
           'beta.c','beta.e',
           'predicted',
           #'Pears_resid'
           'fit',
           'fit_new'
           )
occ.mod9.list = lapply(seq_along(myspecies), function(x) {
  print(paste0('Working on the ', myspecies[x], ' model'))
  print(paste0("Start time: ", Sys.time()))
  tic()
  mod = jagsUI::jags(data = occ.jags.data.list[[x]],
                     model.file = "occ.mod2.txt",
                     parameters.to.save = params,
                     #inits = list(initsfun(),initsfun()),
                     n.iter=15000,
                     n.adapt=1500,
                     n.thin = 10,
                     n.chains=4,
                     parallel = T)
  toc()
  beepr::beep()
  return(mod)
  #42min with 15,000 iters, 1500 adapt, 2 chains, parallel=F
  #20min with 10,000 iters, 1500 adapt, 3 chains, parallel=T (good mixing)
})
names(occ.mod9.list) = myspecies

save(occ.jags.data.list, 
     occ.model.data,
     y.list,
     betavar.array.list,
     occ.mod9.list,
     occ.betavars,
     file = "occ.mod9.output.rdata")


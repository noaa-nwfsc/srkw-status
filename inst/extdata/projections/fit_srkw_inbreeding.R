library(mgcv)
library(R2jags)
library(dplyr)

mcmc_iter = 50000
# Load in the whale data
data(orca)

# Expand the data to large format, to include survival/fecundity
# over time
whaleData = kwdemog::expand(orca)
whaleData$region = ifelse(whaleData$pod %in%
    c("J001","K001","L001"), "SRKW", "NRKW")

# A handful of unknown sexes need to be randomly filled in
# Doesn't really matter because these are all juveniles, and survival rates are the same
whaleData$sexF1M2[which(whaleData$sexF1M2==0)] =
  sample(c(1,2), size=length(which(whaleData$sexF1M2==0)),
    replace=T)

# To simplify a few things from the workshop,
# I'm using a simpler approach and just allowing survival / fecundity to be age based and
whales_since76 = as.character(whaleData$animal[whaleData$birth > 1975])

# NRKW data doesn't exactly line up, last year missing
whaleData$alive[which(whaleData$region== "NRKW" & whaleData$year>2018)] = NA

# convert ages to stages
ages2stages = read.csv("inst/extdata/ages2stages.csv", stringsAsFactors = FALSE)
ages2stages$sex = ages2stages$sexF1M2
whaleData = dplyr::left_join(whaleData, ages2stages)

# add in inbreeding data
genetics = read.csv("inst/extdata/projections/inbreeding coefs with 95 ci and pedigree and wang r values and homozygosity Feb10 2017.csv")
genetics = dplyr::select(genetics, het, Offspring) %>%
  dplyr::mutate(first_char = substr(Offspring,1,1)) %>%
  rename(animal = Offspring) %>%
  dplyr::mutate(nchar = nchar(as.character(animal)))
genetics$animal = as.character(genetics$animal)
genetics$animal[which(genetics$nchar==2)] = paste0(substr(genetics$animal[which(genetics$nchar==2)],1,1), "00",substr(genetics$animal[which(genetics$nchar==2)],2,2))
genetics$animal[which(genetics$nchar==3)] = paste0(substr(genetics$animal[which(genetics$nchar==3)],1,1), "0",substr(genetics$animal[which(genetics$nchar==3)],2,3))
genetics = dplyr::select(genetics, -nchar,-first_char)
whaleData = dplyr::left_join(whaleData, genetics)

# add in new genetics data
genetics = read.csv("inst/extdata/projections/indHet.autosomes.kw151.snp.final.minGQ20.maxmDP17.recode.012.csv",
  stringsAsFactors = FALSE)
genetics = genetics[,c("inds","indHet")]
genetics = genetics[7:107,]
genetics = dplyr::mutate(genetics, first_char = substr(inds,1,1)) %>%
  rename(animal = inds) %>%
  dplyr::mutate(nchar = nchar(as.character(animal)))
genetics$animal = as.character(genetics$animal)
genetics$animal[which(genetics$nchar==2)] = paste0(substr(genetics$animal[which(genetics$nchar==2)],1,1), "00",substr(genetics$animal[which(genetics$nchar==2)],2,2))
genetics$animal[which(genetics$nchar==3)] = paste0(substr(genetics$animal[which(genetics$nchar==3)],1,1), "0",substr(genetics$animal[which(genetics$nchar==3)],2,3))
genetics = dplyr::select(genetics, -nchar,-first_char)
whaleData = dplyr::left_join(whaleData, genetics)

whaleData$old_het = whaleData$het
whaleData$new_het = whaleData$indHet

dplyr::group_by(whaleData,animal) %>%
  summarize(m = new_het[1], pop = population[1]) %>%
  group_by(pop) %>%
  summarize(m = max(m,na.rm=T))

# fit model to animals without the het data
nohet = dplyr::filter(whaleData,is.na(het),
  includeFec==1,
  sexF1M2==1,
  alive==1,
  age%in%seq(10,42),
  !is.na(gave_birth))
nohet$region[which(nohet$region=="SRKW")]="aSRKW"

fec_srkw = jagam(gave_birth ~ s(year) + s(age,k=3) + region,
  family = "binomial",
  data = nohet,
  file="fec_nohet.jags")
require(rjags)
load.module("glm") ## improved samplers for GLMs often worth loading
jm <-jags.model("fec_nohet.jags",data=fec_srkw$jags.data,inits=fec_srkw$jags.ini,n.chains=1)
sam <- jags.samples(jm,c("b"),n.iter=30000,thin=10)
jam <- sim2jam(sam,fec_srkw$pregam)
b = as.matrix(as.mcmc.list(sam$b))

# 13 betas, 1 = intercept, 2 = region(NRKW), 3:11 are year, 12:13 are age
b[,2] = 0
bmu = apply(b,2,mean)
btau = 1/apply(b,2,var)
btau[2] = 0.0071

K1 = solve(cov(b[,3:11]))
K2 = solve(cov(b[,12:13]))

het = dplyr::filter(whaleData,!is.na(het),
  includeFec==1,
  sexF1M2==1,
  alive==1,
  age%in%seq(10,42),
  !is.na(gave_birth))
fec_het = jagam(gave_birth ~ s(year) +
    s(age,k=3) + het, family = "binomial",
  data = het,
  file="fec_het.jags")

data2 = list("y"=fec_het$jags.data$y,
  "n"=fec_het$jags.data$n,
  "X"=fec_het$jags.data$X,
  "w"=fec_het$jags.data$w,
  "K1"=K1,
  "K2"=K2,
  "bmu"=bmu,
  "btau"=btau)
jm_het <-jags.model("fec_het_prior.jags",
  data=data2, inits=NULL,n.chains=3)

b_srkw <- jags.samples(jm_het,c("b"),
  n.iter=mcmc_iter,thin=10)
b2 = as.matrix(as.mcmc.list(b_srkw$b))



# read in data that only includes srkw
srkw = dplyr::filter(whaleData, region=="SRKW")

fems_srkw = dplyr::filter(srkw,
    includeFec==1,
    sexF1M2==1,
    alive==1,
    age%in%seq(10,42),
    !is.na(gave_birth))

fec_gam = mgcv::gam(gave_birth ~ age + I(age^2)
  + s(year,k=3) + het,
  family = "binomial",
  data = fems_srkw)

library(glmmTMB)
fec_glm = glmmTMB::glmmTMB(gave_birth ~ age + I(age^2) + het,
  family = "binomial",
  data = fems_srkw)

fec_srkw = jagam(gave_birth ~ s(year) +
    s(age,k=3), family = "binomial",
  data = fems_srkw, file="inst/extdata/projections/temp.jags")
#saveRDS(fec_srkw, "inst/extdata/projections/pregam.RDS")
#fec_srkw = readRDS("inst/extdata/projections/pregam.RDS")

data2 = list("y"=fec_srkw$jags.data$y,
    "n"=fec_srkw$jags.data$n,
    "X"=fec_srkw$jags.data$X,
    "w"=fec_srkw$jags.data$w,
    "K1mean"=as.numeric(K1mean),
    "K1cov"=K1cov,
    "K2mean"=as.numeric(K2mean),
    "K2cov"=K2cov)
jm_srkw <-jags.model("inst/extdata/projections/fec_srkw.jags",
  data=data2, inits=NULL,n.chains=3)

b_srkw <- jags.samples(jm_srkw,c("b"),
    n.iter=mcmc_iter,thin=10)

sim_srkw = sim2jam(sam = b_srkw,
    pregam = fec_srkw$pregam)

# save model for making predictions
saveRDS(sim_srkw, "inst/extdata/projections/srkw_fec_age-year.rds")
saveRDS(b_srkw, "inst/extdata/projections/srkw_fec_age-year_b.rds")


# Re-run the same model, using only the age prior this time -- no year effects
data2 = list("y"=fec_srkw$jags.data$y,
  "n"=fec_srkw$jags.data$n,
  "X"=fec_srkw$jags.data$X,
  "w"=fec_srkw$jags.data$w,
  "K1mean"=rep(0, length(as.numeric(K1mean))),
  "K1cov"=K1cov,
  "K2mean"=K2mean,
  "K2cov"=K2cov)
jm_srkw <-jags.model("inst/extdata/projections/fec_srkw.jags",
  data=data2, inits=NULL,n.chains=3)

b_srkw <- jags.samples(jm_srkw,c("b"),
  n.iter=mcmc_iter,thin=10)

sim_srkw = sim2jam(sam = b_srkw,
  pregam = fec_srkw$pregam)

# save model for making predictions
saveRDS(sim_srkw, "inst/extdata/projections/srkw_fec_age.rds")
saveRDS(b_srkw, "inst/extdata/projections/srkw_fec_age_b.rds")

###########################################################
# Survival models here
###########################################################

# extract coefficients as priors
b = readRDS("inst/extdata/projections/nrkw_survival_betas.rds")

# extract coefficients as priors
bmean = apply(b[,1:6],2,mean)
btau = solve(cov(b[,1:6]))
K1mean = apply(b[,7:15],2,mean)
K1cov = solve(cov(b[,7:15]))

# read in data that only includes srkw
srkw = dplyr::filter(whaleData, region=="SRKW")

all_srkw = dplyr::filter(srkw,
  includeSurv==1,
  !is.na(alive))

fit = mgcv::gam(alive ~ s(year) + stage * het,
  family = "binomial",data = all_srkw)

fit = glmmTMB::glmmTMB(alive ~ stage + het + (1|year),
  family = "binomial",data = all_srkw)

surv_srkw = jagam(alive ~ s(year) + stage,
 family = "binomial",
 data = all_srkw,
 file="inst/extdata/projections/temp.jags")
#saveRDS(surv_srkw,"pregam_surv.rds")
#surv_srkw = readRDS("inst/extdata/projections/pregam_surv.rds")

data2 = list("y"=surv_srkw$jags.data$y,
  "n"=surv_srkw$jags.data$n,
  "X"=surv_srkw$jags.data$X,
  "w"=surv_srkw$jags.data$w,
  "bmean"=as.numeric(bmean),
  "bcov"=btau,
  "K1mean"=as.numeric(K1mean),
  "K1cov"=K1cov)
jm_srkw <-jags.model("inst/extdata/projections/surv_srkw.jags",data=data2,
  inits=NULL,n.chains=3)
b_srkw <- jags.samples(jm_srkw,c("b"),
  n.iter=mcmc_iter,thin=10)

sim_srkw = sim2jam(sam = b_srkw,
  pregam = surv_srkw$pregam)

# save model for making predictions
saveRDS(sim_srkw, "inst/extdata/projections/srkw_surv_age-year.rds")
saveRDS(b_srkw, "inst/extdata/projections/srkw_surv_age-year_b.rds")

# Re-run the same model, using only the age prior this time -- no year effects
data2 = list("y"=surv_srkw$jags.data$y,
  "n"=surv_srkw$jags.data$n,
  "X"=surv_srkw$jags.data$X,
  "w"=surv_srkw$jags.data$w,
  "bmean"=as.numeric(bmean),
  "bcov"=btau,
  "K1mean"=rep(0, length(K1mean)),
  "K1cov"=K1cov)
jm_srkw <-jags.model("inst/extdata/projections/surv_srkw.jags",
  data=data2, inits=NULL,n.chains=3)

b_srkw <- jags.samples(jm_srkw,c("b"),
  n.iter=mcmc_iter,thin=10)

sim_srkw = sim2jam(sam = b_srkw,
  pregam = surv_srkw$pregam)

# save model for making predictions
saveRDS(sim_srkw, "inst/extdata/projections/srkw_surv_age.rds")
saveRDS(b_srkw, "inst/extdata/projections/srkw_surv_age_b.rds")

# fecundity model in stan
dat = dplyr::filter(whaleData,
  includeFec==1,
  sexF1M2==1,
  alive==1,
  age%in%seq(10,42), !is.na(gave_birth))
dat$year = dat$year - min(dat$year) + 1
dat$srkw = ifelse(dat$region=="SRKW",2,1)
nyear = max(dat$year)

dat$het = dat$new_het
dat$het[which(is.na(dat$het))] = 0

index_nodat = which(dat$het==0)
index_animal = as.numeric(as.factor(dat$animal[index_nodat]))
index_n = max(index_animal)
n_nodat = length(index_nodat)

library(rstan)
standat = list(N = nrow(dat), birth = dat$gave_birth,
  age = dat$age, age2 = dat$age^2, age3 = dat$age^3,
  year = dat$year, nyear = nyear, srkw = dat$srkw,
  index_animal = index_animal, index_n = index_n,
  index_nodat = index_nodat, n_nodat = n_nodat,
  het = dat$het)

fit = stan(file="fecstan.stan",
  data=standat, warmup=3000,iter=5000, chains=3,
  pars = c("age1p","age2p","bhet","sigma","intercept",
    "effects"))
pars = extract(fit)


# survival model in stan
dat = dplyr::filter(whaleData,
  includeSurv==1,
  !is.na(alive),
  !is.na(stage))
dat$year = dat$year - min(dat$year) + 1
dat$srkw = ifelse(dat$region=="SRKW",2,1)

#dat$stage[which(is.na(dat$stage))] = "calf"
dat$numstate = as.numeric(as.factor(dat$stage))
nyear = max(dat$year)

dat$het = dat$new_het
dat$het[which(is.na(dat$het))] = 0

index_nodat = which(dat$het==0)
index_animal = as.numeric(as.factor(dat$animal[index_nodat]))
index_n = max(index_animal)
n_nodat = length(index_nodat)

library(rstan)
standat = list(N = nrow(dat), alive = dat$alive,
  stages = dat$numstate,
  year = dat$year, nyear = nyear, srkw = dat$srkw,
  index_animal = index_animal, index_n = index_n,
  index_nodat = index_nodat, n_nodat = n_nodat,
  het = dat$het)

fit = stan(file="survstan.stan",
  data=standat, warmup=3000,iter=5000, chains=3,
  pars = c("bstage","bhet","sigma","intercept",
    "effects"), control = list(adapt_delta=0.99,
      max_treedepth=20))
pars = extract(fit)

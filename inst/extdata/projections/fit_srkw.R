library(mgcv)
library(R2jags)

mcmc_iter = 50000
# Load in the whale data
data(orca)

# Expand the data to large format, to include survival/fecundity
# over time
whaleData = kwdemog::expand(orca)
whaleData$region = ifelse(whaleData$pod %in% c("J001","K001","L001"), "SRKW", "NRKW")

# A handful of unknown sexes need to be randomly filled in
# Doesn't really matter because these are all juveniles, and survival rates are the same
whaleData$sexF1M2[which(whaleData$sexF1M2==0)] = sample(c(1,2), size=length(which(whaleData$sexF1M2==0)), replace=T)

# To simplify a few things from the workshop,
# I'm using a simpler approach and just allowing survival / fecundity to be age based and
whales_since76 = as.character(whaleData$animal[whaleData$birth > 1975])

# NRKW data doesn't exactly line up, last year missing
whaleData$alive[which(whaleData$region== "NRKW" & whaleData$year>2018)] = NA

# convert ages to stages
ages2stages = read.csv("inst/extdata/ages2stages.csv", stringsAsFactors = FALSE)
ages2stages$sex = ages2stages$sexF1M2
whaleData = dplyr::left_join(whaleData, ages2stages)

# extract coefficients as priors
b = readRDS("inst/extdata/projections/nrkw_fecundity_betas.rds")

# note that the first b coef not included -- that's the intercept/level
K1mean = apply(b[,2:10],2,mean)
K1cov = solve(cov(b[,2:10]))
K2mean = apply(b[,11:12],2,mean)
K2cov = solve(cov(b[,11:12]))

# read in data that only includes srkw
srkw = dplyr::filter(whaleData, region=="SRKW")

fems_srkw = dplyr::filter(srkw,
    includeFec==1,
    sexF1M2==1,
    alive==1,
    age%in%seq(10,42),
    !is.na(gave_birth))

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



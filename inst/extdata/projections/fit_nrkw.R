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

# run survival models for NRKW
# filter northern residents
all_nrkw = dplyr::filter(whaleData,
  includeSurv==1,
  population=="NRKW",
  !is.na(alive))

surv_nrkw = jagam(alive ~ s(year) + stage,
  family = "binomial",
  data = all_nrkw,
  file="inst/extdata/projections/surv_nrkw.jags")

jm_nrkw <-jags.model("inst/extdata/projections/surv_nrkw.jags",
  data=surv_nrkw$jags.data,
  inits=surv_nrkw$jags.ini,
  n.chains=3)
b_nrkw <- jags.samples(jm_nrkw,c("b"),
  n.iter=mcmc_iter, thin=10)

# save raw posteriors
b = as.matrix(as.mcmc.list(b_nrkw$b))
saveRDS(b, file="inst/extdata/projections/nrkw_survival_betas.rds")

# run fecundity models for NRKW
# filter northern residents
all_nrkw = dplyr::filter(whaleData,
  includeFec==1,
  population=="NRKW",
  age >= 10,
  age <= 42,
  sexF1M2==1,
  !is.na(alive))

fec_nrkw = jagam(gave_birth ~ s(year) + s(age,k=3),
  family = "binomial",
  data = all_nrkw,
  file="inst/extdata/projections/fec_nrkw.jags")

jm_nrkw <-jags.model("inst/extdata/projections/fec_nrkw.jags",
  data=fec_nrkw$jags.data,
  inits=fec_nrkw$jags.ini,
  n.chains=3)
b_nrkw <- jags.samples(jm_nrkw,c("b"),
  n.iter=mcmc_iter, thin=10)

b = as.matrix(as.mcmc.list(b_nrkw$b))
saveRDS(b, file="inst/extdata/projections/nrkw_fecundity_betas.rds")

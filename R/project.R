#' Function to project SRKW population in the future using individual based models
#'
#' @param whale_data the expanded data frame, returned by `expand()`
#' @param seed the random seed used to initialize simulations, defaults to `123`
#' @param n_years the length of the time series projections, defaults to 30 time steps
#' @param n_iter the number of iterations per scenario, defaults to 200
#' @param scenarios a character string, or vector of strings describing which scenario
#' to be run. Should be entered as a range of years separated by a colon,
#' e.g. "1981:1985" to use rates from 1981 to 1985 or ("1981:1985","1981:2021")
#' @param verbose whether to print status updates to screen, defaults to FALSE
#' @param p_female The probability of a female birth, defaults to the empirical value since 1980
#' @param n_births The number of births used to calculate `p_female`. By default, this defaults to the empirical
#' value observed since 1980 the probability of female birth. Larger values result in more precise estimates of `p_female`
#' 
#' @import dplyr
#' @import mgcv
#' @importFrom stats rnorm runif rbeta
#' @importFrom stringr str_locate
#'
#' @export
project = function(whale_data, seed = 123, 
                      n_years = 30,
                      scenarios = c("1981:2020"),
                      n_iter = 200, 
                      verbose = FALSE,
                      p_female = NA,
                      n_births = NA) {
  set.seed(seed)
  
  data = whale_data
  
  # First fit models with survival
  whale_data[which(whale_data$animal %in% c("L095","L098","L112") & whale_data$alive==0),] = NA
  # filter out animals born > 1960
  whale_data = dplyr::filter(whale_data, birth > 1960)
  # fill in missing sexes randomly
  whale_data$sexF1M2[which(whale_data$sexF1M2==0)] = sample(1:2, size=length(which(whale_data$sexF1M2==0)), replace=T)
  whale_data = dplyr::rename(whale_data, sex = sexF1M2)
  surv_fit <- gam(alive ~ s(year) + age+I(age^2)+sex,data=dplyr::filter(whale_data,includeSurv==1),
             family="binomial")
  
  # Second fit models to fecundity
  whale_data = data
  whale_data = dplyr::filter(whale_data, sexF1M2==1, 
                            includeFec==1, 
                            age >= 10,
                            age<=43,
                            alive == 1,
                            !is.na(gave_birth))
  fec_fit <- gam(gave_birth ~ s(year) + age+I(age^2),data=whale_data,
             family="binomial")
  
  total_popsize = list()
  repro_females = list()
  
  # do the projections here
  for(s in 1:length(scenarios)) {
    indx = str_locate(scenarios[s],":")
    if(!is.na(indx[1])) {
      year_start = as.numeric(substr(scenarios[s],1,indx[1]-1))
      year_end = as.numeric(substr(scenarios[s],indx[2]+1, nchar(scenarios[s])))
      years_to_sample = seq(year_start, year_end)
    }
    
    # These were the sex ratios we'd used during the workshops -- this is
    # a random variable, rather than fixed, so the uncertainty is propagated
    srb = dplyr::filter(data,birth>1980) %>%
      dplyr::group_by(animal) %>% 
      dplyr::summarize(sex=sexF1M2[1]) %>% 
      dplyr::filter(sex!=0)
    if(is.na(p_female)) {
      p_female = table(srb$sex)[1]/nrow(srb)
    }
    if(is.na(n_births)) {
      n_births = nrow(srb)
    }
    beta_1 = p_female*(n_births+1)
    beta_2 = (1-p_female)*(n_births+1)

    popSize = array(0, dim=c(n_iter, n_years))
    popFems = array(0, dim=c(n_iter, n_years))
    currentPop = dplyr::filter(data, 
                               year == max(whale_data$year),
                               alive==1)
    for(i in 1:n_iter) {
      # get the current age / sex / pod structure of the population
      # to make this more general, I called "populationSegment" = pod, and "breedingGroup" = matriline
      initPopCurrent = currentPop[,c("animal", "pod", "matriline", "sexF1M2", "age")] # animal, pod, matriline, sex, age
      initPopCurrent = dplyr::rename(initPopCurrent, sex=sexF1M2)
      initPopCurrent$animal = as.character(initPopCurrent$animal)
      initPopCurrent$matriline = as.numeric(as.factor(initPopCurrent$matriline)) # , -29 because of NRs
      initPopCurrent$dad = NA # keep track of dads
      
      numWhale = dim(initPopCurrent)[1]
      initPopCurrent$hasCalf = 0
      initPopCurrent$hasCalf[which(as.character(currentPop$animal) %in% as.character(currentPop$mom[which(currentPop$age==1)]))] = 1
      initPopCurrent$mom = as.character(currentPop$mom)
      
      initPopCurrent$sex = as.numeric(initPopCurrent$sex)
      # assign sexes to unsexed individuals [0s]
      initPopCurrent$sex[which(initPopCurrent$Sex==0)] = ifelse(runif(length(which(initPopCurrent$sex==0))) < rbeta(length(which(initPopCurrent$sex==0)), beta_1, beta_2), 1, 2)
      #initPopCurrent$sex[which(initPopCurrent$Sex==0)] = ifelse(runif(length(which(initPopCurrent$sex==0))) < rnorm(length(which(initPopCurrent$sex==0)),p_female, sd_female), 1, 2)
      newID = 9999
      
      for(yrs in 1:dim(popSize)[2]) {
        # first step is find females available to give birth
        indx = which(initPopCurrent$sex == 1 & as.numeric(initPopCurrent$age) >= 10 & as.numeric(initPopCurrent$age) < 43 & initPopCurrent$hasCalf == 0)
        # second step is to calculate predicted fecundity rates and make fecundity stochastic - every individual's pregnancy is a coin flip
        
        # if a range of years is used to draw demographic rates from sample those randomly
        initPopCurrent$year = as.numeric(sample(c(as.character(years_to_sample)), size=1))
        
        if(length(indx) > 0) {
          # bind together the current moms with the matching fecundity data
          p_fec = mgcv::predict.gam(fec_fit, initPopCurrent[indx,], type="response")
          
          pregMoms = indx[which(runif(length(p_fec)) < p_fec)]
          # third step is to make moms that aren't mate limited give birth to calves of known sexes
          if(length(pregMoms) > 0) {
            
            for(ll in 1:length(pregMoms)) {
              # loop over moms and only let them breed if there's a mature male in a different matriline
              dads = which(initPopCurrent$sex == 2 & as.numeric(initPopCurrent$age) > 12)
              #dads = which(initPopCurrent$Sex == 2 & initPopCurrent$breedingGroup != initPopCurrent$breedingGroup[pregMoms[ll]] & as.numeric(initPopCurrent$age1) > 12)
              if(length(dads) > 0) {
                # assign the pod / matriline to be the same of the mom
                newpod = initPopCurrent$pod[pregMoms[ll]]
                newmat = initPopCurrent$matriline[pregMoms[ll]]
                # sex is stochastic
                newsex = ifelse(runif(1) < rbeta(1, beta_1, beta_2), 1, 2)
                # bookkeeping
                newage = 0
                newcalf = 0
                newmom = initPopCurrent$animal[pregMoms[ll]]
                
                # sample from potential dads in proprtion to their estimated relative reproductive output
                # their ages are initPopCurrent$age1[dads], and the fecundity model defined above outside of loops
                # sample from potential dads in proprtion to their estimated relative reproductive output
                # their ages are initPopCurrent$age1[dads], and the fecundity model defined above outside of loops
                
                #inflation_factor = 5
                #pred.male.fecundity[32:200] = pred.male.fecundity[31]
                #probs = pred.male.fecundity[as.numeric(initPopCurrent$age[dads])]
                #probs[which.max(initPopCurrent$age1[dads])] = probs[which.max(initPopCurrent$age[dads])] * inflation_factor
                #newdad = initPopCurrent$animal[sample(dads,1, prob=probs)]
                newdad = initPopCurrent$animal[sample(dads,1)]
                newID = newID + 1
                # add calves to the population
                newdf = data.frame(animal=newID,
                                   pod=newpod,matriline=newmat,sex=newsex,
                                   age=newage,dad=newdad,hasCalf=newcalf,
                                   mom=newmom,year=initPopCurrent$year[1])
                initPopCurrent = rbind(initPopCurrent, newdf)
                
              }# end if(length(dads) > 0) {
            } # end ll loop
          } # end if(length(pregMoms)
          
        }
        # bookkeeping: update whales that have calves
        initPopCurrent$hasCalf[which(initPopCurrent$hasCalf==1)] = 0
        initPopCurrent$hasCalf[pregMoms] = 1
        
        # step 4 is calculate predicted survival at age
        initPopCurrent$sexF1M2 = as.numeric(initPopCurrent$sex)
        p_surv = mgcv::predict.gam(surv_fit, initPopCurrent, type="response")

        initPopCurrent = dplyr::select(initPopCurrent, -sexF1M2)
        
        # step 5: stochastic survival to kill whales
        liveOrDie = rep(1,length(p_surv))
        dead = which(runif(length(p_surv)) > p_surv)
        liveOrDie[dead] = 0
        
        # step 6: see if any of these dead animals has a calf - if so, kill the calf too
        if(length(which(liveOrDie == 0)) > 0) {
          for(ll in 1:length(dead)) {
            # kill the calf
            if(is.na(initPopCurrent$hasCalf[dead[ll]]) == FALSE & initPopCurrent$hasCalf[dead[ll]] == 1) {
              liveOrDie[which(initPopCurrent$mom == dead[ll])] = 0
            }
          }
        }
        
        # step 7: bookkeeping at the end of the time step
        # first remove dead animals from the population
        if(length(dead) > 0 ) deadWhales = initPopCurrent[which(liveOrDie==0),]  ## MIKE ADDITION - list of who is dead
        if(length(dead) > 0) initPopCurrent = initPopCurrent[-which(liveOrDie==0),]
        # second age remaining whales
        initPopCurrent$age = as.numeric(initPopCurrent$age) + 1
        # third record pop size to later calculate recovery targets
        popSize[i,yrs] = dim(initPopCurrent)[1]
        popFems[i,yrs] = dim(dplyr::filter(initPopCurrent, sex==1,
                                           age <= 42, age >= 10))[1]
        if(verbose==TRUE) print(paste0("Iteration ",i, " : year ", yrs))
      } # this is end of yrs loop
      
    } # end i loop 
    
    total_popsize[[s]] = popSize
    repro_females[[s]] = popFems
  }
  
  return(list(fec_fit = fec_fit, surv_fit = surv_fit,
              total_popsize = total_popsize, 
              repro_females = repro_females))
}


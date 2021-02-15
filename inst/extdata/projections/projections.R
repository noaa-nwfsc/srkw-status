  # scenario used affects which rates are used for projections
  years_to_sample = unique(whaleData$year)
  if(scenario=="last5") years_to_sample = seq(year.end-4,year.end)
  if(scenario=="best") years_to_sample = 1985:1989

  # These were the sex ratios we'd used during the workshops -- this is
  # a random variable, rather than fixed, so the uncertainty is propagated
  pFemale = 45/100 # this is the number of female births/total up through 2018
  sdFemale = sqrt(pFemale*(1-pFemale)/100)

  nSims = 1000 # number of simulations
  NYEARS = 25 # number of years to project
  #years_to_sample= sort(unique(whaleData$year), decreasing=T)[1:5]
  popSize = array(0, dim=c(nSims, NYEARS))
  currentPop = whaleData[which(whaleData$year==year.end & whaleData$region =="SRKW" &
      is.na(whaleData$alive) == FALSE &
      whaleData$alive != 0 & whaleData$region =="SRKW"),]

  # These arrays / matrices are to keep track of things we might be interested in storing
  popSize.L = popSize # L pod population
  popSize.J = popSize # J pod population
  popSize.K = popSize # K pod population
  nMales = popSize # total males
  nFemales = popSize # reproductive females

  for(i in 1:nSims) {
#print(i)
      # get the current age / sex / pod structure of the population
      # to make this more general, I called "populationSegment" = pod, and "breedingGroup" = matriline
      initPopCurrent = currentPop[,c("animal", "pod", "matriline", "sexF1M2", "age")] # animal, pod, matriline, sex, age
      initPopCurrent = dplyr::rename(initPopCurrent, sex=sexF1M2)
      initPopCurrent$populationSegment = as.numeric(as.factor(initPopCurrent$pod)) # 3 becomes L pod, -14 because of NRs
      initPopCurrent$animal = as.character(initPopCurrent$animal)
      initPopCurrent$matriline = as.numeric(as.factor(initPopCurrent$matriline)) # , -29 because of NRs
      initPopCurrent$dad = NA # keep track of dads

      numWhale = dim(initPopCurrent)[1]
      initPopCurrent$hasCalf = 0
      initPopCurrent$hasCalf[which(as.character(currentPop$animal) %in% as.character(currentPop$mom[which(currentPop$age==1)]))] = 1
      initPopCurrent$mom = as.character(currentPop$mom)

      initPopCurrent$sex = as.numeric(initPopCurrent$sex)
      # assign sexes to unsexed individuals [0s]
      initPopCurrent$sex[which(initPopCurrent$Sex==0)] = ifelse(runif(length(which(initPopCurrent$sex==0))) < rnorm(length(which(initPopCurrent$sex==0)),pFemale, sdFemale), 1, 2)
      newID = 9999

      for(yrs in 1:dim(popSize)[2]) {
        # first step is find females available to give birth
        indx = which(initPopCurrent$sex == 1 & as.numeric(initPopCurrent$age) >= 10 & as.numeric(initPopCurrent$age) < 43 & initPopCurrent$hasCalf == 0)
        # second step is to calculate predcited fecundty rates and make fecundity stochastic - every individual's pregnancy is a coin flip

        if(length(indx) > 0) {
        predicted_fecundityLogit = predict(fecundity.model, newdata = data.frame("age"=as.numeric(initPopCurrent$age[indx]), "region"="SRKW", year = sample(years_to_sample, size=length(indx), replace=T)), se.fit=TRUE)
        # Make the individuals have correlated rates by drawing from quantile
        #predicted_fecundityAtAge = plogis(rnorm(length(indx), predicted_fecundityLogit$fit, predicted_fecundityLogit$se.fit))
        predicted_fecundityAtAge = plogis(qnorm(runif(1), mean = predicted_fecundityLogit$fit, sd = predicted_fecundityLogit$se.fit, lower.tail = TRUE, log.p = FALSE))

        pregMoms = indx[which(runif(length(indx)) < predicted_fecundityAtAge[as.numeric(initPopCurrent$age[indx])-9])]
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
              newsex = ifelse(runif(1) < rnorm(1, pFemale, sdFemale), 1, 2)
              # bookkeeping
              newage = 0
              newcalf = 0
              newmom = initPopCurrent$animal[pregMoms[ll]]

              # sample from potential dads in proprtion to their estimated relative reproductive output
              # their ages are initPopCurrent$age1[dads], and the fecundity model defined above outside of loops
              #inflation_factor = 3
              #pred.male.fecundity[32:200] = pred.male.fecundity[31]
              #probs = pred.male.fecundity[as.numeric(initPopCurrent$age1[dads])]
              #probs[which.max(initPopCurrent$age1[dads])] = probs[which.max(initPopCurrent$age1[dads])] * inflation_factor
              #newdad = initPopCurrent$animal[sample(dads,1, prob=pred.male.fecundity[as.numeric(initPopCurrent$age1[dads])])]

              # Dads are just assigned randomly -- this is slightly different from the paternity modeling we've done with Ford et al.
              newdad = sample(dads,1)
              newID = newID + 1
              # add calves to the population
              initPopCurrent = rbind(initPopCurrent, c(newID, newpod,newmat, newsex,newage,newdad,newcalf,newmom))

            }# end if(length(dads) > 0) {
          } # end ll loop
        } # end if(length(pregMoms)

        }
        # bookkeeping: update whales that have calves
        initPopCurrent$hasCalf[which(initPopCurrent$hasCalf==1)] = 0
        initPopCurrent$hasCalf[pregMoms] = 1

        # step 4 is calculate predicted survival at age
        ages2stages$sex = ages2stages$sexF1M2
        initPopCurrent$sex = as.numeric(initPopCurrent$sex)
        initPopCurrent$age = as.numeric(initPopCurrent$age)
        newStage = dplyr::left_join(initPopCurrent, ages2stages)
        predicted_survivalLogit = predict(survival.model, newdata = data.frame("stage"=newStage$stage,"region"="SRKW", year = sample(years_to_sample, size=nrow(newStage), replace=T)), se.fit=TRUE)
        #predSurv = plogis(rnorm(length(indx), predicted_survivalLogit$fit, predicted_survivalLogit$se.fit))
        predSurv = plogis(qnorm(runif(1), mean = predicted_survivalLogit$fit, sd = predicted_survivalLogit$se.fit, lower.tail = TRUE, log.p = FALSE))

        # step 5: stochastic survival to kill whales
        liveOrDie = rep(1,length(predSurv))
        dead = which(runif(length(predSurv)) > predSurv)
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

        popSize.L[i,yrs] = length(which(initPopCurrent$pod=="L001"))
        popSize.K[i,yrs] = length(which(initPopCurrent$pod=="K001"))
        popSize.J[i,yrs] = length(which(initPopCurrent$pod=="J001"))
        nFemales[i,yrs] = length(which(initPopCurrent$sex==1 & initPopCurrent$age < 43  & initPopCurrent$age >= 10))
        nMales[i,yrs] = length(which(initPopCurrent$sexF1M2==2))
      } # this is end of yrs loop

  }

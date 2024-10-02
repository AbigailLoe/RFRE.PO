pacman::p_load(matrixStats, dplyr, tidyr, LaplacesDemon)

#Generate recurrent events with gap times ~ Exp(lambda_i)
#lambda_i can be a vector or scalar
#s is span of study
# n is the number of observations/people in subgroup of study.
# n_j is the number of recurrent events simulated per subject

# need some way of getting after the final event...

gen_xia_lambda = function(study_size, numcovs) {
  sigma_x = 4
  numlevels = 3
  X = sample(1:numlevels, size = study_size, replace = TRUE)
  mean_model = function(x) {
    beta  = matrix(rep(c(3 / 4, 1 / 3, 1 / 2), each = study_size),
                   nrow = study_size,
                   ncol = 3)
    lambda = NULL
    for (i in 1:study_size) {
      lambda = c(lambda, beta[i, X[i]])
    }
    return(lambda)
  }
  
  dat = as.data.frame(cbind(mean_model(X), X))
  colnames(dat) = c("lambda", paste("X", 1:numcovs, sep = "_"))
  return(dat)
}

#x ~ exp(lambda) where fx = lambda e ^-lambda x
unifLambda = function(study_size, numcovs, lowerbound, upperbound){
  X = runif(study_size, min = lowerbound, max = upperbound)
  lambda = X
  return(cbind.data.frame(lambda, X))
}



##### TestGen3 #####
#this one also works better than TestGen2

testGen3 = function(study_size,numcovs) {
  X =matrix(runif(numcovs * 2*study_size, min = .02, max = 4),
            nrow = 2*study_size,
            ncol = numcovs)
  mean_model <- function(x) {
    beta  <- c(rep(2, numcovs))
    sqbeta = c( rep(-1.5, numcovs))
    interbeta = c( rep(.4, numcovs))
    sum(x * beta, x * x * sqbeta, x[1]*x *interbeta)
  }
  
  dat <- as.data.frame(cbind(apply(X, 1, mean_model), X))
  colnames(dat) <- c("Y", paste("X", 1:numcovs, sep = "_"))
  indices = c(which(dat$Y>=quantile(dat$Y, .5)))
  #gets the two extreme ends.
  dat = dat[indices,]
  dat$lambda = 1.5*exp(dat[["Y"]])
  dat2 = dat[ , -which(names(dat) %in% c("Y"))]
  return(dat2)
}

#12


restrictLambda3= function(desiredSize, upperBound, lowerBound){
  dat = testGen3(4* desiredSize, 4)
  dat = dat[which(dat$lambda< upperBound), ]
  dat = dat[which(dat$lambda > lowerBound), ]
  
  while(nrow(dat)<desiredSize){
    dat2 = testGen3(4*desiredSize, 4)
    dat2 = dat2[which(dat2$lambda< upperBound), ]
    dat2 = dat2[which(dat2$lambda > lowerBound), ]

    
    dat = rbind(dat, dat2) 
  }
  dat = dat[1:desiredSize, ]
  return(dat)
}

gen6mod = function(study_size,numcovs = 7, lower = 8/15, upper = 15) {
  X1 = rnorm(2*study_size, mean = 0, sd = 2)
  X2 = rnorm(2*study_size, 2, .8)
  X3 = rpois(2*study_size, 4)
  X4 = rbeta(2*study_size, 7, 1)+.1
  X5 = sample(c(2, 1, 0), size = 2*study_size, 
              replace = TRUE, prob = c(1/6, 1/6, 2/3))
  X6 = sample(c(-5, -2, 3, 2), size = 2*study_size, 
              replace = TRUE, prob = c(1/10, 1/3, 1/5, 11/30))
  X7 = sample(c(0:4), size = 2*study_size,
              replace = T)
  # X8 = rbeta(2*study_size, 4, 45)
  dat = cbind.data.frame(X1, X2, X3, X4, X5, X6, X7,  
                         exp(X2*sin(X1/X6))+ ifelse(X2>2, -X3, X3)+ 3*X1*X6 + X2^2*X4)
  names(dat) = c(paste("X", 1:7, sep = "_"), "lambda")
  dat = dat %>% filter(lambda <=upper & lambda >= lower)
  return(dat)
}

gen6= function(study_size){
  data = gen6mod(10*study_size)
  while(nrow(data)<study_size){
    data = rbind(data, gen6mod(10*study_size))
  }
  return(data[1:study_size, ])
}



updatedS1 = function(study_size, numcovs = 50, lower = 1/7, upper = 10) {
  X1 = rnorm(study_size, 2, 1.5)
  X2 = runif(study_size, min = 1,  max= 4)
  X3 = rgamma(study_size, 5)
  
  S = toeplitz((study_size:1)/study_size)
  R = rinvwishart(study_size, S)
  X4 = apply(R, 2, function(x) median(x))
  
  X5 = 1/(2.5+ rpois(study_size, 10))
  X6 = 3*rbeta(2*study_size, 4, 1)
  X7 = runif(study_size, min = -4, max = 1)
  x.sig = cbind(X1, X2, X3, X4, X5, X6, X7)
  
  x.noise =matrix(rnorm((numcovs-7) * study_size, mean = .02, sd = 6),
            nrow = study_size,
            ncol = numcovs-7)
  
  Y = X1*X2/X3 + sqrt(X5)- X4*sin(X4) + ifelse(X6 < X1, X1/X6, X6/X1)  + X7*X1^2 - X2*X4^2
 
  dat <- cbind.data.frame(Y, x.sig, x.noise)
  
  colnames(dat) <- c("Y", paste("X", 1:numcovs, sep = "_"))
  dat$lambda = exp(dat[["Y"]])
  #range(dat$lambda)
  dat2 = dat %>% dplyr::select(-Y) %>% filter(lambda<= upper & lambda >= lower)
  #hist(dat2$lambda)
  return(dat2)
}

restrictS1 = function(study_size, numcovs, low, up){
  print("started first iteration")
  data = updatedS1(2*study_size, numcovs, lower = low, upper = up)
  my.iter = 1
  while(nrow(data)<study_size){
    data = rbind(data, updatedS1(2*study_size, numcovs, lower = low, upper = up))
    print(my.iter)
    my.iter = my.iter +1
  }
  return(data[1:study_size, ])
}

# test = restrictS1(500, 7, 1/8, 4)
# test= test %>% mutate(trueRelationship =
#                         X_1*X_2/X_3 + sqrt(X_5)- X_4*sin(X_4) -
#                         ifelse(X_6 < X_1, X_1/X_6, X_6/X_1) + X_7*X_1^2 - X_2*X_4^2,
#                       fooVar = ifelse(X_6 < X_1, X_1/X_6, X_6/X_1))
# hist(test$trueRelationship)
# hist(test$lambda)
# plot(test$lambda~test$trueRelationship)
# View(test[which(test$trueRelationship< -10),])


empiricalCorrelation= function(tT){
  corMat = matrix(NA, nrow= ncol(tT), ncol = ncol(tT))
  for(i in 1:ncol(tT)){
    for(j in 1:i){
      corMat[i,j] = corMat[j,i] = cor(tT[,i], tT[,j])
    }
  }
  return(corMat)
}


gen_Sliding = function(tList = T_t_format,
                       covMat = covs,
                       censorTime = L,
                       checks = checks,
                       numcovs = numcovs) {
  # print(paste0("empirical average correlation: ", )
  
  
  T_t = data.frame(cbind(1:nrow(tList$T_t),covMat,censorTime, tList$T_t))
  colnames(T_t) = c("ID", colnames(covMat), "L", 1:ncol(tList$T_t))
  
  longCal = T_t %>% 
    pivot_longer(!c(ID, paste0("X_", 1:numcovs), L), values_to = "timeToEvent") %>% 
    na.omit()
  
  longCal$marker = as.numeric(as.character(
    factor(longCal$name, levels=c(1:length(checks)), 
                          labels=checks)))
  longCal$calTime = as.numeric(as.character(longCal$marker)) + longCal$timeToEvent
  
  #refactor to represent calendar time
  
  longCal$delta = ifelse(longCal$calTime< longCal$L, 1, 0)
  longCal = subset(longCal, select = -c(name) )

  return(longCal)
}


gen_R_star = function(lambdas, correl, study_size, 
                      study_length = study_length,
                      discard_excess = F, covs = covs, space) {
  indic = 0
  foo = 1
  nJparam = 1000
  if(correl >.2){
    nJparam = 1600
  }
  if(correl >.6){
    nJparam = 2000
  }
  while (indic == 0 ) {
    nJparam = nJparam * 2 ^ (foo - 1)
    group = generate_data_woloss(
      n = study_size,
      lambda_i = lambdas,
      rho = correl,
      s = (study_length), #think about why is s larger than the study size??
      # I don't think it should be necessarily huge!!
      J = nJparam,
      covs = covs,
      discard_excess = discard_excess,
      space = space
    )
    # print(paste0("Tries: ", foo, ". nJ: ", nJparam))
    foo = foo + 1
    indic = group$use_dataset
  }
  # View(cbind(rowMeans2(group$R), lambdas))
  R_star = group$R_star
  covs = group$covs
  lambdas = group$lambdas
  prev_tte = group$prev.window.tte

  R_star_long = R_star
  # R_star_long[R_star_long > 100*study_length] = NA  
  R_star[R_star > study_length] = NA
  R = group$R

  
  return(list(R_star= R_star, covs = covs, lambdas = lambdas, R = R, R_extra = R_star_long, prev_tte = prev_tte))
}


generate_data_woloss = function(n, lambda_i, rho, s, J,
                                discard_excess = T, covs =covs, space){
  use_dataset = 1
  # J is the number of gaptimes simulated per subject
  sigma = diag(rep(1, J)) + (1 - diag(rep(1, J))) * rho #set up covariance matrix
  # working correctly as of 3/8

  Y = mvrnorm(n, mu = rep(0, J), sigma) #simulate multivariate normal distribution
  #View(rowMeans(Y))
  #this command is taking a lot more time than expected.
  U = pnorm(Y) #transform to correlated uniform distributions
  # mean(rowMeans(U))
  R = - 1/lambda_i*log(U) # exponentials
  R_star_orig <- rowCumsums(R)
  
  #View(cbind(1/lambda_i, rowMeans(R)))
  badIndex = NULL
  
  burn_in_time =3* 1/(min(lambda_i))
  # burn_in_time = 0
  print(burn_in_time)
  
  
  R_star = R_star_orig - burn_in_time
  
  R_prev = R_star_orig - burn_in_time+space
  
  badIndex = which( R_star[, J]< s) # finds rows that have the final
  # observation not making it to end of study AFTER accounting for burn-in
  
  #need these for filtering for the lambdas AND R_star
  prop_discard = length(badIndex)/nrow(R_star)
  print(paste0("discarded ", prop_discard, "of observations"))
  if(length(badIndex) > 0 ){
    R_star = R_star[-badIndex,]
    lambda_i = lambda_i[-badIndex]
    covs = covs[-badIndex, ]
  }
  if (sum(R_star[, J] < s) > 0) {
    use_dataset = 0
  }
  R_star[R_star <= 0] = NA
  R_prev[R_prev <= 0] = NA
  R_prev[R_prev >space] = NA
  TT_prev = rowMins(R_prev, na.rm = TRUE)
  TT_prev[is.infinite(TT_prev)] = space
  # View(cbind(TT_prev, R_prev))
  return(list(
    # corresponding covariates.
    lambdas = lambda_i,
    covs = covs,
    R = R_star_orig,
    R_star = R_star,
    use_dataset = use_dataset,
    U = U,
    prev.window.tte = TT_prev
  )
  )
}



#Format data for simulation with Loss of follow up
# NEED TO CHANGE TO FIX THE PROBLEM OF ESTIMAND 2 NOT RELYING ON TAU

get_min_T_tau_wloss <- function(R_star, L, n, b, window , t) {
  # L is the censoring time
  #get T(t) for each window
  # as of 11/29, t is prespecified and equal for all people.
  # EVENTUALLY get t to be non unique across individuals.
  T_t = array(NA, c(n, b))
  delta = array(NA, c(n, b))
  for (j in 1:b) {

    R_star_j = R_star - t[j] #increments down a window.
    #look at cumulative time, remove a specified observation window.
    R_star_j[R_star_j <= 0] = NA #recurrent events observed before time t[j]
    L_j = L - t[j] #distance to censoring time
    L_j[L_j <= 0] = NA
    #take the minimum of the distance to censoring time, and the distance to event time.
    R_star_j_prime = ifelse(is.na(L_j), NA, rowMins(cbind(L_j, R_star_j), na.rm =
                                                      TRUE))
    #can be greater than tau.
    
    foo <- cbind(L_j, R_star_j_prime, R_star_j)
    #row-wise summary of finding the minimum of the window.
    T_t[, j] = R_star_j_prime
    L_j[is.na(L_j)] <- NA
    #if they are equal, then we either don't observe an event in that window, or we have already been
    # censored.
    delta[, j] = ifelse(all.equal(R_star_j_prime, L_j), 0, 1)
    foo2 <- cbind(R_star_j_prime, L_j, delta[, j])
    #if the minimum of the calendar time for a j window is equal to the censoring time,
    # then we set the delta value to 0
    # else, it gets set to 1.
  }
  # min_T_tau <- apply(T_t, c(1,2), function(x) min(x, tau))
  min_T_tau <- T_t
  min_T_tau[min_T_tau > window] = window
  return(list(
    T_t = T_t,
    min_T_tau = min_T_tau,
    delta = delta
  ))
}


formatTt <- function(R_star, L, n , t) {
  b = length(t)
  T_t = array(NA, c(n, b))
  delta = array(NA, c(n, b))
  for (j in 1:b) {
    
    R_star_j = R_star - t[j] #increments down a window.
    #look at cumulative time, remove a specified observation window.
    R_star_j[R_star_j <= 0] = NA #recurrent events observed before time t[j]
    L_j = L - t[j] #distance to censoring time
    L_j[L_j <= 0] = NA
    #take the minimum of the distance to censoring time, and the distance to event time.
    R_star_j_prime = ifelse(is.na(L_j), NA, rowMins(cbind(L_j, R_star_j), na.rm =
                                                      TRUE))
    #can be greater than tau.
    
    foo <- cbind(L_j, R_star_j_prime, R_star_j)
    #row-wise summary of finding the minimum of the window.
    T_t[, j] = R_star_j_prime
    L_j[is.na(L_j)] <- NA
    #if they are equal, then we either don't observe an event in that window, or we have already been
    # censored.
    # delta[, j] = if_else(isTRUE(all.equal(R_star_j_prime, L_j, tolerance = (.Machine$double.eps)*3)), 0, 1)
    delta[, j] = as.numeric(!near(L_j, R_star_j_prime))
    foo2 <- cbind(R_star_j_prime, L_j, delta[, j])
    #if the minimum of the calendar time for a j window is equal to the censoring time,
    # then we set the delta value to 0
    # else, it gets set to 1.
  }
  # min_T_tau <- apply(T_t, c(1,2), function(x) min(x, tau))
  return(list(
    T_t = T_t,
    delta = delta
  ))
}






genRMST = function(min_T_tau = min_T_tau,
                   study_size,
                   covariateMat = covs,
                   censorTime = L,
                   status = delta,
                   numWindows = b,
                   checks = checks,
                   tau = tau) {
  data_wide = data.frame(cbind(1:study_size, min_T_tau, covariateMat))
  colnames(data_wide) = c("ID",
                          paste0("min_T_tau_", 1:numWindows),
                          paste0("X_", 1:ncol(covariateMat)))
  data_wide$censorTime = censorTime
  data_long = gather(data_wide,
                     tj,
                     min_T_tau,
                     paste0("min_T_tau_", 1:numWindows),
                     factor_key = TRUE)
  data_long$timing = paste0("t_", substr(as.character(data_long$tj), 11, 11))
  data_long = data_long[!is.na(data_long$min_T_tau), ]
  delta_long = array(status)
  data_long$delta = delta_long[!is.na(delta_long)]
  pseudoval = NULL
  meas = NULL
  j = NULL
  for (j in 1:numWindows) {
    subdata = data_long[data_long$tj == paste0("min_T_tau_", j), ]
    meas_time = rep(checks[j], nrow(subdata))
    pseudoval_each = pseudomean(subdata$min_T_tau, subdata$delta, tmax =
                                  tau)
    pseudoval = c(pseudoval, pseudoval_each)
    meas = c(meas, meas_time)
  }
  data_long$pseudoval = pseudoval
  data_long$meas = meas
  return(data_long)
  
}






rmse <- function(a, b) {
  sqrt(mean((a - b) ^ 2))
}

myPseudo = function(dat = data_long, tau ){
  rowVecCheck = NULL
  dat$pseudoEst = rep(NA, nrow(dat)) #initialize the empty column 
  s = unique(dat$marker) #get the landmarks that we want to observe probabilities at
  for( si in s){
    #print(si)
     siDat= dat %>% filter(marker == si) %>%
      as.data.frame() 
      kmWithAll = survfit(Surv(siDat$timeToEvent, siDat$delta)~1, 
              se.fit =FALSE, type ="kaplan")
    hatS= summary(kmWithAll, tau, extend = TRUE)$surv # gets hatS^tau for all people
    n = nrow(siDat)
    for( i in 1:nrow(siDat)){
      # print(paste("person I: ", i))
      missingIDi = siDat[-i, ]
      sfit_noi = survfit(Surv(missingIDi$timeToEvent, missingIDi$delta)~1,
                         se.fit = FALSE, type = "kaplan")
      hatSnoi = summary(sfit_noi, tau, extend = TRUE)$surv
      IDi = siDat$ID[i]
      pseudoval = n * hatS - (n-1)*hatSnoi
      correctRow = which(dat$ID == IDi & dat$marker == si)
      #print(correctRow)
      rowVecCheck = c(rowVecCheck, correctRow)
      dat$pseudoEst[correctRow] = pseudoval #I think this needs to be correct!
    }
  }
  return(list(returnedData = dat, debugRow = rowVecCheck))
}






## need to write a function that determines the bias based on the lambdas.
# our hazard is exponential. That means that survival time is 
# exp(-lambda * (currTime - window) )

# debugging:
# 
# View(generatedDF)
# survivalWindow = window.set
# data = generatedDF
# # to get the prediction, move to Sim_setting file.
# predictionDF = finalPredictions

survivalBias = function(lambdas, generatedDF, survivalWindow, predictionDF){
  # survivalWindow: tau in P(T>tau).
  # lambdas: from the gen.data function
  # generatedDF: from gen.data
  # predictionDF is from the Sim_Setting file and what happens when we run
  # predictions.
  
  # need to get the individual level bias on average
  # average across individuals in a dataset.
  # average across all replicates (to be done in the testing file)

  myDF = cbind.data.frame(generatedDF, predictionDF)
  expectedSurvival = exp(-lambdas * survivalWindow)
  head(expectedSurvival )
  
  ID = 1:length(expectedSurvival)
  expectedDF = cbind.data.frame(expectedSurvival,
                                ID)
  head(expectedDF)
  useDF = myDF %>% 
    left_join(expectedDF, by= c("ID"= "ID")) %>% 
    group_by(ID) %>% 
    # filter(delta == 1) %>%  #gets only the observed events. Is this something we
    # want to consider?
    add_count(ID) %>% 
    dplyr::select(ID, magic, naive, hrf, expectedSurvival) %>% 
    group_by(ID) %>% 
    mutate(
           magicError = mean(magic-expectedSurvival),
           naiveError = mean(naive-expectedSurvival),
           hrfError = mean(hrf-expectedSurvival)
           
           ) %>% filter(row_number()==1)
  biasRes = colMeans(useDF[, c("magicError",
                               
                               "hrfError",
                               
                               "naiveError"
                               )])
  
  #"magicAUC", "rangerAUC",  "hrfAUC", "merfAUC","naiveAUC"
    
  #how should we deal with a different number of observations per person?
  
  return(biasRes)
}



c_stat = function(predicted, actual){
  if(length(predicted) != length( actual)){
    return("OOPS")
  }else {
    n = length(predicted)
    index.mat = t(combn(n , 2))
    concord.pairs = 0
    discord.pairs = 0
    comp.pairs = 0
    for( i in 1:nrow(index.mat)){
      pos.i = index.mat[i, 1]
      pos.j = index.mat[i ,2]
      if((predicted[pos.i] == predicted[pos.j]) || (actual[pos.i] == actual[pos.j])){
        comp.pairs = comp.pairs + 1
      }
      else{
        if((predicted[pos.i] >= predicted[pos.j]) && 
           (actual[pos.i] >= actual[pos.j])){
          concord.pairs = concord.pairs + 1
        }else if ((predicted[pos.i] < predicted[pos.j]) && 
                  (actual[pos.i] < actual[pos.j])){
          concord.pairs = concord.pairs + 1
        } else{
          discord.pairs = discord.pairs + 1
        }
        }
    }
    return( (concord.pairs ) / (concord.pairs + discord.pairs))}
}


# Pseduo code: check to see if the length of the vectors is the same
# if it is: then we sample from the predicted, the calTime, and the delta
# if both observations are censored, skip it
# if one observation is censored



c_censored = function(predicted, calTime, observed){
  # predicted will be the pseudo-observation
  # calTime will be the T_i, not the actual time of survival
  # observed is the delta_i which indicates when things happen
  if((length(predicted) != length(calTime))|| (length(calTime) != length(observed))){
    return("OOPS")
  }else {
    n = length(predicted)
    index.mat = t(combn(n , 2))
    concord.pairs = 0
    discord.pairs = 0
    comp.pairs = 0
    for( i in 1:nrow(index.mat)){
      pos.i = index.mat[i, 1]
      pos.j = index.mat[i ,2]
      # calculate times that we don't want to do anything or count pairs
      if((predicted[pos.i] == predicted[pos.j]) ||  # if predicted times are the same
         (calTime[pos.i] == calTime[pos.j]) || # if the calendar times are the same
         sum(observed[c(pos.i, pos.j)]) == 0||  #if they are both censored
         (calTime[pos.i]< calTime[pos.j] & observed[pos.i] == 0)|| # if the earlier event is censored first
         (calTime[pos.j]< calTime[pos.i] & observed[pos.j] == 0) # if the earlier event is censored first
         ){
        comp.pairs = comp.pairs + 1
      }
      else{
        if((predicted[pos.i] >= predicted[pos.j]) && 
           (calTime[pos.i] >= calTime[pos.j])){
          concord.pairs = concord.pairs + 1
        }else if ((predicted[pos.i] < predicted[pos.j]) && 
                  (calTime[pos.i] < calTime[pos.j])){
          concord.pairs = concord.pairs + 1
        } else{
          discord.pairs = discord.pairs + 1
        }
      }
    }
    return( (concord.pairs ) / (concord.pairs + discord.pairs))}
}


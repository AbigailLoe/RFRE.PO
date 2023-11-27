####### Small useful meta functions ########
# in general, the only particularly illuminating function is gen6mod
# as it contains the data generation mechanism.

pacman::p_load(matrixStats, dplyr, tidyr, LaplacesDemon)


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


# Gets the prespecified cehck ins and formats the data correctly.
# returns a data frame with the calendar times, not yet in PO format.
gen_Sliding = function(tList = T_t_format,
                       covMat = covs,
                       censorTime = L,
                       checks = checks,
                       numcovs = numcovs) {
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


# creates the correlated exponential gap times and corresponding calendar times.

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
  R_star = group$R_star
  covs = group$covs
  lambdas = group$lambdas
  prev_tte = group$prev.window.tte

  R_star_long = R_star 
  R_star[R_star > study_length] = NA
  R = group$R
  return(list(R_star= R_star, covs = covs, lambdas = lambdas, R = R, R_extra = R_star_long, prev_tte = prev_tte))
}

# generates data without loss to follow up.
generate_data_woloss = function(n, lambda_i, rho, s, J,
                                discard_excess = T, covs =covs, space){
  use_dataset = 1
  # J is the number of gaptimes simulated per subject
  sigma = diag(rep(1, J)) + (1 - diag(rep(1, J))) * rho #set up covariance matrix

  Y = mvrnorm(n, mu = rep(0, J), sigma) #simulate multivariate normal distribution
  U = pnorm(Y) #transform to correlated uniform distributions
  # mean(rowMeans(U))
  R = - 1/lambda_i*log(U) # exponentials
  R_star_orig <- rowCumsums(R)
  
  #View(cbind(1/lambda_i, rowMeans(R)))
  badIndex = NULL
  burn_in_time =3* 1/(min(lambda_i))
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


# formats to get residual survival time conditional on survival up to
# current time point or check-in.

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

# create pseudo-values!

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



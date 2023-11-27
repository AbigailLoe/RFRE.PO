#########################################################
# A file for generating the data necessary for simulation
#########################################################

pacman:: p_load(foreach, doParallel,
                tidyr, dplyr,
                survival, geeM, gdata,
                MASS, pseudo, matrixStats,
                LaplacesDemon)

source("internalFxns.R")

# a function to generate the recurrent event data necessary for simulation.
gen.data = function(study_length = 5,
                    space = 1/2, window = 1,
                    num_covs = 7,
                    study_size = 10 ^ 3
                    ,p.censoring = 1 / 14,
                    correl = 0,
                    estimand = 1, 
                    discard_excess = T){
  prelimN = 2*study_size
  checks = seq(from = 0, to = study_length-window, by = space)
  
  b = length(checks)
  l_m = gen6(prelimN) # data generation mechanism.
  
  lambdas = l_m[,"lambda"]
  length(lambdas)
  range(lambdas)
  range(1/lambdas)
  covs = l_m[, - which(names(l_m)=="lambda")] #assumes lambdas at end!!!!

  # get the recurrent event calendar times.
  foo = gen_R_star(lambdas = lambdas, 
                   correl= correl, 
                   study_size = prelimN, 
                   study_length = study_length,
                   discard_excess = discard_excess,
                   covs = covs,
                   space = space)
  # now we have roughly 2 times as many R_star as we need; this is necessary for 
  # highly correlated cases.
  
  R_star=foo$R_star
  covs = foo$covs
  lambdas = foo$lambdas
  R = foo$R
  prev_event = foo$prev_tte
  extendedTime  = foo$R_extra
  dim(R)
  length(lambdas)
  
  R_star = R_star[1:study_size, ]
  covs = covs[1:study_size, ]
  extendedTime = extendedTime[1:study_size, ]
  lambdas = lambdas[1:study_size]
  prev_event = prev_event[1:study_size]

  
  if(!all.equal(length(lambdas), study_size, nrow(covs),
     nrow(R_star))){
    return("something went wrong here. Throwing a string as an error")
  }
  
  if (p.censoring == 0) {
    L = rep_len(x = study_length, length.out = study_size)
  } else{
    L = rexp(study_size, p.censoring) 
    L[L > study_length] = study_length
  }

  
  extendedTime[extendedTime<L] = NA

  test = cbind.data.frame( nextEvent = rowMins(extendedTime, na.rm = TRUE),ID = 1:study_size)
  
  print("number of events before second window:")
  print((table(rowSums2(R_star<space, na.rm = TRUE))))
  maxEvent = cbind.data.frame(ID= 1:study_size, nextEvent=rowSums2(R_star<space, na.rm = TRUE)+1)


 {
    T_t_format = formatTt(R_star = R_star, L= L, 
                          n = study_size, t = checks)
    
    tT = T_t_format$T_t
    for(i in 1:length(checks)){
      tT[,i] = tT[,i]+checks[i]
    }

   # get the window start times.
   
    data_long=gen_Sliding(tList = T_t_format,
                          covMat= covs,
                          censorTime = L,
                          checks = checks,
                          numcovs = num_covs)
    
   #final check of delta values.
    data_long = data_long %>% 
        mutate(delta = as.numeric(!near(calTime, L))) 


    # get pseudo observations.
    interim = myPseudo(dat = data_long, tau = window)
    
    data_long = interim[["returnedData"]]
    # Do a final check of censored vs not:
    
      
    data_long = data_long %>% na.omit()
    
    avgLambda = 1/mean(lambdas)
    
    # add a new row so that we can have a sense of "history"
        dat2 = data_long %>% 
      as_tibble() %>%
      group_by(ID) %>% 
      group_modify(~ add_row(.x,.before=0))%>% 
      mutate(timeToEventPermissible = if_else(timeToEvent <= space, timeToEvent,
                                              space)) %>%
      mutate(timeToEventImpute = if_else(is.na(timeToEventPermissible), 
                                        rexp(n = 1, mean(lambdas)),
                                         timeToEventPermissible),
             timeToEventImpute = if_else(timeToEventImpute <= space, 
                                         timeToEventImpute,
                                         space), #take the min of a and the imputed variable
        avg_exacerb_impute =lag(cummean(timeToEventImpute))) %>%
      na.omit()
    
    prev_event_df = cbind.data.frame(ID = 1:study_size, prev_event)
    
    tester = dat2 %>% 
      as_tibble() %>%
      group_by(ID) %>% 
      group_modify(~ add_row(.x,.before=0))
    
    # now need to add in the past event time.
    for(i in 1:length(unique(tester$ID))){
      temp.df= tester %>% filter(ID == i)
      temp.df$timeToEventPermissible[1] = prev_event_df[i, "prev_event"]
      tester[which(tester$ID == i), ] = temp.df
    }
     
    tester = tester %>% 
      group_by(ID) %>% 
      mutate(avg_exacerb_true =lag(cummean(timeToEventPermissible))) %>% 
      na.omit()
    
    dat3 = tester %>% left_join(test, by = "ID")
    censoredRows = which(dat3$delta==0)
    dat3$timeToEvent[censoredRows]=dat3$nextEvent[censoredRows]-dat3$marker[censoredRows]
    
  }
  return(list(data = dat3, 
              generatingLambdas = lambdas, 
              Tt = T_t_format$T_t,
              nextEvent = maxEvent
              ))
}


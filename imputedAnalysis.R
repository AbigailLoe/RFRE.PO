pacman::p_load(dplyr, tidyr, htree, geepack, stats, miceadds, xtable, Hmisc)

getwd()
# set working directory to where the appropriate files are kept
# setwd("~/Dropbox (University of Michigan)/LLZ/Simulation/dataApplication/imputedCat")
m = 10
set.seed(16)


# need to do a 70/30 split for all of the datasets
hrf.param2 = hrf.param = hrf.se = hrf.se2 =
  summer.param = summer.vcov =wald.param = wald.vcov=
  c.param = c.vcov =  c.hrf1 = c.hrf2 = c.wald= c.summer = 
  c.c =c.hrf1.se = c.hrf2.se = c.wald.se = c.summer.se = c.c.se = list()
for(l in 1:10){
  print(l)
  # do the permutation test with 50 iterates each
  runData = read.csv(paste0("runData", l, ".csv")) %>% 
    mutate(rateExacerb = rateExacerb*100,
           trtgroup = if_else(trtgroup == 2, 0, 1), # re-level treatment.
           # create indicators that may be of use.
           gold2 = if_else(goldclass == 2, 1, 0),
           sleep1 = if_else(sleep_dysfunction >= 1, 1, 0),
           panic4= if_else(runData$sudden_feelings_panic>3, 1, 0)
           )
  set.seed(16)
  folds = if_else(cut(seq(1, 1035), breaks = 10, labels = F)<=7,
                  1, 0) # randomly sample teh training adn testing data.
  idpool = runData$sub.id
  id.sel=idpool[folds==1] # takes only 70% of data
  
  # create training and validation sets
  train.data = runData[runData$sub.id  %in% id.sel, ]
  test.data = runData[!(runData$sub.id  %in% id.sel), ]
  
  selnHRF2 = unique(c( setdiff(names(runData), 
                               c("timeToEvent", "delta", "censor",  
                                 "pseudoEst", "calTime", "toDelete", "missingVal",
                                 "numpills", "ID", "sub.id", "marker", "numExacerb",
                                 "POFEV1_","visit", "firstTime", "contTimeSince", 
                                 "contTimeSince2", "timeSince"))
  )); selnHRF2
  
  #########  HTREE #########
  trainDf2 = as.data.frame(train.data[,c("pseudoEst", selnHRF2)])
  testDF2 = as.data.frame(test.data[, c("pseudoEst", selnHRF2)])
  set.seed(16)
  htree.rf2=hrf(x=trainDf2,time=train.data$checkIn,
                id=as.numeric(train.data$sub.id),
                yindx="pseudoEst",ntrees=500,
                historical=FALSE, se=FALSE,
                control = list(nodesize = 40))
  
  set.seed(16)
  pred.hrf2 =  predict_hrf(htree.rf2,x= testDF2,time=test.data$checkIn,
                           id=as.numeric(test.data$sub.id),
                           all.trees=FALSE,se=FALSE)
  
  
  # variable important data frames are computationally intensive to calculate
  # I have attached the datasframes as they have been created with 100 permutations
  varImportDF2 = read.csv(paste0("varImportDF2_", l, ".csv"))
  
  # need a variance covariance matrix
  # need a parameter estimate.
  varImportDF2 = varImportDF2 %>%
    mutate(estimate = `Marginalized.error`-`Model.error`,
           se.2 = abs(Relative.change/ `Z.value`) ,
           se.w = abs(estimate/Z.value))%>%
    filter( `Z.value` >= 1.96 & estimate >0 & se.2>0)
  # print(varImportDF2)
  hrf.param2[[l]] = t(data.frame(estimate = varImportDF2$Relative.change,
                                 row.names = varImportDF2$Predictor))
  
  
  hrf.se2[[l]] = t(data.frame(se = varImportDF2$se.2,
                              row.names = varImportDF2$Predictor))
  c.hrf2[[l]] = rcorr.cens(pred.hrf2, Surv( test.data$timeToEvent, test.data$delta))[["C Index"]]
  
  # calculate bootstrap standard errors.
  blankSE = NULL
  for(i in 1:100){
    print(i/100)
    myboot.ind = sample(1:length(pred.hrf2), replace = T)
    blankSE =c(blankSE, rcorr.cens(pred.hrf2[myboot.ind], 
                                   Surv( test.data$timeToEvent[myboot.ind], 
                                         test.data$delta[myboot.ind]))[["C Index"]])
  }
  c.hrf2.se[[l]] = sd(blankSE)
  
  print(c.hrf2)
  
  # wald forward selection model
  currString = "pseudoEst~oneMonth+bronchitis+gender+trtgroup+
      beta_blocker_non_selective+fev1fvcpct_+SGRQ_symptoms+oxygen+hosp1yr+
      rateExacerb+icslama+icslaba+beta_blocker_use"
  formula = as.formula(currString)
  curr.geebin <- geese(formula = formula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  summary(curr.geebin)$mean
  
  wald.param[[l]]= t(summary(curr.geebin)$mean["estimate"])
  wald.vcov[[l]]= curr.geebin$vbeta.ajs
  
  cov.test.gee=as.matrix(cbind(rep(1,nrow(test.data)),
                               test.data[, c("oneMonth",
                                             "bronchitis",
                                             "gender",
                                             "trtgroup",
                                             "beta_blocker_non_selective",
                                             "fev1fvcpct_",
                                             "SGRQ_symptoms",
                                             "oxygen",
                                             "hosp1yr","rateExacerb",
                                             "icslama",
                                             "icslaba",
                                             "beta_blocker_use")]))
  
  betaMat=matrix(curr.geebin$beta)
  currPred = 1/(1+exp(-cov.test.gee %*% betaMat))
  c.wald[[l]] = rcorr.cens(currPred, Surv( test.data$timeToEvent, test.data$delta))[["C Index"]]
  
  blankSE = NULL
  for(i in 1:100){
    print(i/100)
    myboot.ind = sample(1:nrow(test.data), replace = T)
    cov.test.gee=as.matrix(cbind(rep(1,nrow(test.data)),
                                 test.data[myboot.ind, c("oneMonth",
                                                         "bronchitis",
                                                         "gender",
                                                         "trtgroup",
                                                         "beta_blocker_non_selective",
                                                         "fev1fvcpct_",
                                                         "SGRQ_symptoms",
                                                         "oxygen",
                                                         "hosp1yr","rateExacerb",
                                                         "icslama",
                                                         "icslaba",
                                                         "beta_blocker_use")]))
    currPred = 1/(1+exp(-cov.test.gee %*% betaMat))
    blankSE = c(blankSE, rcorr.cens(currPred, 
                                    Surv( test.data$timeToEvent[myboot.ind], 
                                          test.data$delta[myboot.ind]))[["C Index"]])
  }
  
  c.wald.se[[l]] = sd(blankSE)
  
  
  
  train.data.baseline = train.data %>%
    group_by(sub.id) %>%
    mutate(fev1_baseline = first(fev1_liters),
           age10 = age/10)
  
  test.data = test.data %>%
    group_by(sub.id) %>%
    mutate(fev1_baseline = first(fev1_liters),
           age10 = age/10)
  summer.geebin <- geese(pseudoEst ~ trtgroup +nowsmk+fev1_baseline + gender +age10,
                         data=train.data.baseline,scale.fix=TRUE,family=gaussian,
                         id=train.data.baseline$sub.id,
                         mean.link="logit",corstr="unstructured",jack=TRUE)
  #
  summer.param[[l]]= t(summary(summer.geebin)$mean["estimate"])
  summer.vcov[[l]]= summer.geebin$vbeta.ajs
  
  
  cov.test.gee=as.matrix(cbind(rep(1,nrow(test.data)),
                               test.data[, c("trtgroup",
                                             "nowsmk",
                                             "fev1_baseline",
                                             "gender",
                                             "age10")]))
  
  betaMat=matrix(summer.geebin$beta)
  summerPred = 1/(1+exp(-cov.test.gee %*% betaMat))
  
  c.summer[[l]] = rcorr.cens(summerPred, Surv( test.data$timeToEvent, test.data$delta))[["C Index"]]
  blankSE = NULL
  for(i in 1:100){
    print(i/100)
    myboot.ind = sample(1:nrow(test.data), replace = T)
    cov.test.gee=as.matrix(cbind(rep(1,nrow(test.data)),
                                 test.data[myboot.ind, c("trtgroup",
                                                         "nowsmk",
                                                         "fev1_baseline",
                                                         "gender",
                                                         "age10")]))
    currPred = 1/(1+exp(-cov.test.gee %*% betaMat))
    blankSE = c(blankSE, rcorr.cens(currPred, 
                                    Surv( test.data$timeToEvent[myboot.ind], 
                                          test.data$delta[myboot.ind]))[["C Index"]])
  }
  c.summer.se[[l]] = sd(blankSE)
  
}


wald.res = pool_mi(qhat= wald.param, u = wald.vcov)
summer.res =(pool_mi(qhat= summer.param, u = summer.vcov))
# c.res =(pool_mi(qhat= c.param, u = c.vcov))

wald.df = summary(wald.res) %>% 
  as.data.frame %>% 
  mutate(var.name = c("Intercept", "oneMonth", "bronchitis",
                      "gender", "trtgroup", "beta_blocker_non_selective",
                      "fev1fvcpct_", "SGRQ_symptoms", "oxygen",
                      "hosp1yr","rateExacerb", "icslama", "icslaba",
                      "beta_blocker_use"),
         model.type = "wald") %>% 
  relocate(var.name, .before = results)

summer.df = summary(summer.res) %>% 
  as.data.frame %>% 
  mutate(var.name = c("Intercept", "trtgroup",
                      "nowsmk",
                      "fev1_baseline", 
                      "gender",
                      "age10"),
         model.type = "summer") %>% 
  relocate(var.name, .before = results)




# gets all of the possibly significant variables
allPossVars =
  Reduce(union, list(colnames(hrf.param2[[1]]),
                     colnames(hrf.param2[[2]]),
                     colnames(hrf.param2[[3]]),
                     colnames(hrf.param2[[4]]),
                     colnames(hrf.param2[[5]]),
                     colnames(hrf.param2[[6]]),
                     colnames(hrf.param2[[7]]),
                     colnames(hrf.param2[[8]]),
                     colnames(hrf.param2[[9]]),
                     colnames(hrf.param2[[10]])))

# go back into the variable importance dataframes and access them.
for(l in 1:10){
  varImportDF2 = read.csv(paste0("varImportDF2_", l, ".csv"))
  
  # need a variance covariance matrix
  # need a parameter estimate.
  varImportDF2 = varImportDF2 %>%
    mutate(estimate = `Marginalized.error`-`Model.error`,
           se.2 = abs(Relative.change/ `Z.value`) ,
           se.w = abs(estimate/Z.value))%>%
    filter(Predictor %in% allPossVars)
  # print(varImportDF2)
  hrf.param[[l]] = t(data.frame(estimate = varImportDF2$Relative.change,
                                row.names = varImportDF2$Predictor))
  
  
  hrf.se[[l]] = t(data.frame(se = varImportDF2$se.2,
                             row.names = varImportDF2$Predictor))
}


vcov.hrf = list((hrf.se[[1]])^2 %*% diag(1, nrow=(length(hrf.se[[1]]))),
                (hrf.se[[2]])^2%*%diag(1, nrow=(length(hrf.se[[2]]))),
                (hrf.se[[3]])^2%*%diag(1, nrow=(length(hrf.se[[3]]))),
                (hrf.se[[4]])^2%*%diag(1, nrow=(length(hrf.se[[4]]))),
                (hrf.se[[5]])^2%*%diag(1, nrow=(length(hrf.se[[5]]))),
                (hrf.se[[6]])^2%*%diag(1, nrow=(length(hrf.se[[6]]))),
                (hrf.se[[7]])^2%*%diag(1, nrow=(length(hrf.se[[7]]))),
                (hrf.se[[8]])^2%*%diag(1, nrow=(length(hrf.se[[8]]))),
                (hrf.se[[9]])^2%*%diag(1, nrow=(length(hrf.se[[9]]))),
                (hrf.se[[10]])^2%*%diag(1, nrow=(length(hrf.se[[10]]))))



hrf.same.var.df = summary(pool_mi(qhat= hrf.param, u = vcov.hrf)) %>% 
  as.data.frame() %>% 
  mutate(var.name = allPossVars,
         model.type = "hrfSameVars")


total.imputed.df = rbind(
  summer.df,
  wald.df,
  # c.df,
  hrf.same.var.df
) %>% as.data.frame(.) %>% 
  rename(lower = "(lower",
         upper = "upper)")

row.names(total.imputed.df) = NULL





c.wald.combine =(pool_mi(qhat= c.wald, u = (unlist(c.wald.se))^2%*%diag(1, nrow=10)))
c.hrf2.combine =(pool_mi(qhat= c.hrf2, u = (unlist(c.hrf2.se))^2%*%diag(1, nrow=10)))
c.summer.combine =(pool_mi(qhat= c.summer, u = (unlist(c.summer.se))^2%*%diag(1, nrow=10)))
c.results.df = rbind.data.frame(wald.c.value = summary(c.wald.combine),
                                hrf.imputed.c.value = summary(c.hrf2.combine),
                                summer.c.value = summary(c.summer.combine))




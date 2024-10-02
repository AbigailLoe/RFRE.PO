rm(list = ls())
pacman::p_load(dplyr, tidyr, htree, geepack, stats, miceadds, xtable, Hmisc, 
               viridisLite, survival)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("internalFxns.R")
m = 10
set.seed(16)


source("~/University of Michigan Dropbox/Abby Loe/LLZ/Simulation/internalFxns.R")

# need to do a 70/30 split for all of the datasets
hrf.res = wald.res= summer.res = as.data.frame(matrix(NA, nrow = 0, ncol = 6))
c.hrf2.se = c.wald.se = c.summer.se = c.c.se = 
  c.vcov = summer.vcov = wald.vcov = summer.param = wald.param = hrf.param2 = hrf.param=list()

hrf.iterator = wald.iterator = summer.iterator = 1
names(hrf.res) = names(wald.res) = names(summer.res) = c("impute", "time.point", "c.stat", "c.se", "po.mse", "mse.se")

t.i = c(0, 30, 60, 90, 120, 150, 180)
for(l in 1:10){
  
  print(l)
  # do the permutation test with 50 iterates each
  runData = read.csv(paste0("runData", l, ".csv")) %>% 
    mutate(rateExacerb = rateExacerb*100,
           trtgroup = if_else(trtgroup == 2, 0, 1), # re-level treatment.
           # create indicators that may be of use.
           gold2 = if_else(goldclass == 2, 1, 0),
           sleep1 = if_else(sleep_dysfunction >= 1, 1, 0),
           panic4= if_else(sudden_feelings_panic>3, 1, 0)
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
                                 "contTimeSince2", "timeSince", "rateExacerb", 
                                 "catTimeSince", "noRecentMem", "lastYear", 
                                 "sixMonths", "threeMonths", "oneMonth", "X", 
                                 "experiencedSevere", "hosp1yr"))
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
  test.data$hrf.predictions = as.vector(pred.hrf2)
  # gets the predictions for hrf for a single imputed data set.
  # now need to analyze c-statistic and po.mse at each itme point and get a measure of
  # variablility (bootstraph)
  for(time.point in t.i){
    subdata = test.data %>% filter(checkIn == time.point)
    c.stat = rcorr.cens(subdata$hrf.predictions, Surv( subdata$timeToEvent, subdata$delta))[["C Index"]]
    po.mse = mean((subdata$hrf.predictions-subdata$pseudoEst)^2)
    # c.hrf2[wald.iterator, ] = 
    cSE = mseSE= NULL
    for(i in 1:100){
      print(i/100)
      myboot.ind = sample(1:nrow(subdata), replace = T)
      cSE =c(cSE, rcorr.cens(subdata$hrf.predictions[myboot.ind], 
                             Surv( subdata$timeToEvent[myboot.ind], 
                                   subdata$delta[myboot.ind]))[["C Index"]])
      mseSE = c(mseSE, mean(
        (subdata$hrf.predictions[myboot.ind]-subdata$pseudoEst[myboot.ind])^2)
      )
    }
    c.se = sd(cSE)
    mse.se = sd(mseSE)
    hrf.res[hrf.iterator, ] = c(l, time.point, c.stat, c.se, po.mse, mse.se)
    hrf.iterator = hrf.iterator+1
  }
  
  # wald forward selection model
  currString = "pseudoEst~bronchitis+gender+trtgroup+beta_blocker_non_selective+
  fev1fvcpct_+SGRQ_symptoms+oxygen+icslama+icslaba+beta_blocker_use"
  formula = as.formula(currString)
  curr.geebin <- geese(formula = formula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  summary(curr.geebin)$mean
  
  cov.test.gee=as.matrix(cbind(rep(1,nrow(test.data)),
                               test.data[, c(
                                             "bronchitis",
                                             "gender",
                                             "trtgroup",
                                             "beta_blocker_non_selective",
                                             "fev1fvcpct_",
                                             "SGRQ_symptoms",
                                             "oxygen",
                                             "icslama",
                                             "icslaba",
                                             "beta_blocker_use")]))
  
  betaMat=matrix(curr.geebin$beta)
  currPred = 1/(1+exp(-cov.test.gee %*% betaMat))
  test.data$wald.predictions = currPred
  for(time.point in t.i){
    subdata = test.data %>% filter(checkIn == time.point)
    c.wald = rcorr.cens(subdata$wald.predictions, Surv( subdata$timeToEvent, subdata$delta))[["C Index"]]
    po.mse = mean((subdata$wald.predictions-subdata$pseudoEst)^2)
    cSE = mseSE= NULL
    for(i in 1:100){
      print(i/100)
      myboot.ind = sample(1:nrow(subdata), replace = T)
      cSE =c(cSE, rcorr.cens(subdata$wald.predictions[myboot.ind], 
                             Surv( subdata$timeToEvent[myboot.ind], 
                                   subdata$delta[myboot.ind]))[["C Index"]])
      mseSE = c(mseSE, mean(
        (subdata$wald.predictions[myboot.ind]-subdata$pseudoEst[myboot.ind])^2)
      )
    }
    c.se = sd(cSE)
    mse.se = sd(mseSE)
    wald.res[wald.iterator, ] = c(l, time.point, c.wald, c.se, po.mse, mse.se)
    wald.iterator = wald.iterator+1
  }
  
  
  
  
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
  
  
  cov.test.gee=as.matrix(cbind(rep(1,nrow(test.data)),
                               test.data[, c("trtgroup",
                                             "nowsmk",
                                             "fev1_baseline",
                                             "gender",
                                             "age10")]))
  
  betaMat=matrix(summer.geebin$beta)
  summerPred = 1/(1+exp(-cov.test.gee %*% betaMat))
  
  
  test.data$summer.pred = summerPred
  for(time.point in t.i){
    subdata = test.data %>% filter(checkIn == time.point)
    c.summer = rcorr.cens(subdata$summer.pred, Surv( subdata$timeToEvent, subdata$delta))[["C Index"]]
    po.mse = mean((subdata$summer.pred-subdata$pseudoEst)^2)
    cSE = mseSE= NULL
    for(i in 1:100){
      print(i/100)
      myboot.ind = sample(1:nrow(subdata), replace = T)
      cSE =c(cSE, rcorr.cens(subdata$summer.pred[myboot.ind], 
                             Surv( subdata$timeToEvent[myboot.ind], 
                                   subdata$delta[myboot.ind]))[["C Index"]])
      mseSE = c(mseSE, mean(
        (subdata$summer.pred[myboot.ind]-subdata$pseudoEst[myboot.ind])^2)
      )
    }
    c.se = sd(cSE)
    mse.se = sd(mseSE)
    summer.res[summer.iterator, ] = c(l, time.point, c.summer, c.se, po.mse, mse.se)
    summer.iterator = summer.iterator+1
  }
  
  
}


# now we want to take a dataframe of the 
summer.res$model = "summer"
wald.res$model = "wald"
hrf.res$model = "hrf"

total.res = rbind.data.frame(summer.res, wald.res, hrf.res) %>% 
  group_by(model, time.point) %>% 
  mutate(average.c = mean(c.stat),
         v.t.c = mean((c.se)^2)+ 11/(10*9)*sum((c.stat - average.c)^2),
         average.mse = mean(po.mse),
         v.t.mse = mean((mse.se)^2)+ 11/(10*9)*sum((po.mse - average.mse)^2)
  ) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>% 
  mutate(c.lower.bound = average.c- 1.96*sqrt(v.t.c),
         c.upper.bound = average.c+ 1.96*sqrt(v.t.c),
         mse.lower.bound = average.mse -1.96 * sqrt(v.t.mse),
         mse.upper.bound = average.mse + 1.96 * sqrt(v.t.mse),
         model = as.factor(model),
         time.point2 = time.point+rnorm(nrow(.), 0, .5))

color_palette_mods <- c("wald" = "#482576FF",
                        "hrf" = "#1E9C89FF",
                        "summer" = "#FDE725FF")


c.plot = ggplot(total.res)+  
  geom_errorbar(aes(x = time.point2, ymin = c.lower.bound, 
                    ymax = c.upper.bound, color = model, alpha = .1), 
                show.legend = FALSE, width = 3,
                linewidth = 1.5)+
  geom_line(aes(x = time.point, y = average.c, color = model),
            linewidth = 1.5)+
  geom_point(aes(x = time.point, y = average.c, shape = model),
             show.legend = FALSE,
             size = 2)+
  scale_color_manual(values = color_palette_mods,
                     labels=c("hrf" = "RFRE.PO",
                              "wald" = "Wald Forward Selection Model",
                              "summer" = "Clinical Input Model"),
                     name = "Model"
  )+
  labs(title = "C-Statistic for Each Follow-up Window Without History Variables",
       x = "Follow-up Window Start-Time",
        y = "C-Index")+
  theme(
    legend.position = "bottom"
  )+theme_bw()+ ylim(.45, .68); c.plot

ggsave("~/University of Michigan Dropbox/Abby Loe/LLZ/Simulation/Graphics_for_Paper/azith_c_time_vary_NO_HISTORY.png",
       plot = c.plot, width = 10, height = 8)
ggsave("~/University of Michigan Dropbox/Abby Loe/LLZ/Simulation/Graphics_for_Paper/azith_c_time_vary_NO_HISTORY.pdf",
       plot = c.plot, width = 10, height = 8)
ggsave("~/University of Michigan Dropbox/Abby Loe/LLZ/Simulation/Graphics_for_Paper/azith_c_time_vary_NO_HISTORY.eps",
       plot = c.plot, width = 10, height = 8)

mse.plot = ggplot(total.res)+
  geom_errorbar(aes(x = time.point2, ymin = mse.lower.bound, 
                    ymax = mse.upper.bound, color = model, alpha = .1,
                    width = 3),show.legend = FALSE,
                linewidth = 1.5)+
  geom_line(aes(x = time.point, y = average.mse, color = model),
            linewidth = 1.5)+
  geom_point(aes(x = time.point, y = average.mse, shape = model),
             size = 2,
             show.legend = FALSE)+
  scale_color_manual(values = color_palette_mods,
                     labels=c("hrf" = "RFRE.PO",
                              "wald" = "Wald Forward Selection Model",
                              "summer" = "Clinical Input Model"),
                     name = "Model"
  )+
labs(title = "MSE for Each Follow-up Window without History Variables",
     x = "Follow-up Window Start-Time",
     y = "MSE")+
  theme(
    legend.position = "bottom"
  )+theme_bw()+ ylim(.22, .3); mse.plot

ggsave("~/University of Michigan Dropbox/Abby Loe/LLZ/Simulation/Graphics_for_Paper/azith_mse_time_vary_NO_HISTORY.png",
       plot = mse.plot, width = 10, height = 8)
ggsave("~/University of Michigan Dropbox/Abby Loe/LLZ/Simulation/Graphics_for_Paper/azith_mse_time_vary_NO_HISTORY.pdf",
       plot = mse.plot,width = 10, height = 8)
ggsave("~/University of Michigan Dropbox/Abby Loe/LLZ/Simulation/Graphics_for_Paper/azith_mse_time_vary_NO_HISTORY.eps",
       plot = mse.plot,width = 10, height = 8)


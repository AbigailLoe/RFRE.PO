pacman::p_load(dplyr, tidyr, ggplot2, htree, geepack, 
               Hmisc, utils, foreach, doParallel, stringr)

#load in data and useful functions file.
source("~/Dropbox (University of Michigan)/LLZ/Simulation/broken_LM/internalFxns.R")
runData = read.csv("~/Dropbox (University of Michigan)/LLZ/Simulation/dataApplication/imputedCat/runData1.csv")


runData = runData %>% mutate(rateExacerb = rateExacerb *100)
n=length(unique(runData$sub.id));n
unique(runData$sub.id)
# tail(unique(runData$sub.id), 40)

idpool=unique(runData$sub.id)
n == length(idpool)

set.seed(16)




# all of the longitudinal variables should not be updated
# ie 

# need to format some of the data.
runData$clinic = as.numeric(as.factor(runData$clinic))
runData$age10 = runData$age/10
runData$smokeStat = if_else(runData$nowsmk>0, 1, 0)
runData$trtgroup = runData$trtgroup - 1


selnHRF = unique(c( setdiff(names(runData), 
                 c("timeToEvent", "delta", "censor",  
                   "pseudoEst", "calTime", "toDelete", "missingVal",
                   "numpills", "ID", "sub.id", "marker", "numExacerb",
                   "POFEV1_","visit", "firstTime", "contTimeSince", 
                   "contTimeSince2"))
                 )); selnHRF
# selnHRf contains all of the permissible variables and variables that are not 
# a type of response.

set.seed(16)
folds = if_else(cut(seq(1, length(unique(runData$sub.id))), breaks = 10, labels = F)<=7,
                1, 0)
idpool = runData$sub.id
id.sel=idpool[folds==1] # takes only 70% of data
runData$gold2 = if_else(runData$goldclass == 2, 1, 0)
runData$gold1 = if_else(runData$goldclass == 2, 1, 0)
train.data = runData[runData$sub.id  %in% id.sel, ]
test.data = runData[!(runData$sub.id  %in% id.sel), ]



###### FORWARD SELECTION LOOPS - Wald Statistic #####

forward.seln.vars = setdiff(c(selnHRF, "age10"), c("age10", "checkIn", "marker", "visit" , 
                                                   "clinic", "numExacerb", "timeSince",
                                                   "contTimeSince"))

pList= matrix(data = NA, nrow = length(forward.seln.vars), ncol = 2)

for(i in 1:length(forward.seln.vars)){
  var.name = forward.seln.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0("pseudoEst~", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian, 
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList = pList %>% 
  as.data.frame()

names(pList) = c("variable", "pVal")
pList %>% arrange(as.numeric(pVal)) %>% View()


# numExacerb is number of exacerbations on study time

currString = "pseudoEst~oneMonth"
formula = as.formula(currString)
curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(forward.seln.vars, c("numExacerb", "oneMonth"))
pList2= matrix(data = NA, nrow = length(updated.vars), ncol = 2)
for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian, 
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList2[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList2 = pList2 %>% 
  as.data.frame()
names(pList2) = c("variable", "pVal")
pList2 %>% arrange(as.numeric(pVal)) %>% View()





currString = "pseudoEst~oneMonth+lastYear"
# experiencedSevere, goldclass, hosp1yr, nowsmk, race,
# trtgroup, variables all of the table introduce
# too much model instability!!


# instable.variables = pList2$variable[which(pList2$pVal== 0)]

formula = as.formula(currString)
curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean


currString = "pseudoEst~oneMonth+lastYear"
formula = as.formula(currString)
updated.vars = setdiff(updated.vars, c( "lastYear"))
pList3= matrix(data = NA, nrow = length(updated.vars), ncol = 2)
for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian, 
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList3[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList3 = pList3 %>% 
  as.data.frame()

names(pList3) = c("variable", "pVal")
pList3 %>% arrange(as.numeric(pVal)) %>% View()


currString = "pseudoEst~oneMonth+lastYear+bronchitis"
formula = as.formula(currString)

curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean





updated.vars = setdiff(updated.vars, "bronchitis")
pList4= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

# currString = paste0(currString, "")
for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList4[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList4 = pList4 %>%
  as.data.frame()

names(pList4) = c("variable", "pVal")
pList4 %>% arrange(as.numeric(pVal)) %>% View()

# ok, so

instable.variables4 = pList4$variable[which(pList4$pVal == 0 )]



# table(train.data$pain)
currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, c(instable.variables4,
                                       "gender"))


pList5= matrix(data = NA, nrow = length(updated.vars), ncol = 2)
for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList5[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList5 = pList5 %>%
  as.data.frame()

names(pList5) = c("variable", "pVal")
pList5 %>% arrange(as.numeric(pVal)) %>% View()


table(train.data$trtgroup)
train.data$age10 = train.data$age/10
# final model is!!!
currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "trtgroup")

pList6= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList6[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList6 = pList6 %>%
  as.data.frame()

names(pList6) = c("variable", "pVal")
pList6 %>% arrange(as.numeric(pVal)) %>% View()




currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+beta_blocker_non_selective"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "beta_blocker_non_selective")

pList7= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList7[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList7 = pList7 %>%
  as.data.frame()

names(pList7) = c("variable", "pVal")
pList7 %>% arrange(as.numeric(pVal)) %>% View()




currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+
beta_blocker_non_selective+fev1fvcpct_"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "fev1fvcpct_")

pList8= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList8[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList8 = pList8 %>%
  as.data.frame()

names(pList8) = c("variable", "pVal")
pList8 %>% arrange(as.numeric(pVal)) %>% View()


currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+
beta_blocker_non_selective+fev1fvcpct_+SGRQ_symptoms"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "SGRQ_symptoms")

pList9= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList9[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList9 = pList9 %>%
  as.data.frame()

names(pList9) = c("variable", "pVal")
pList9 %>% arrange(as.numeric(pVal)) %>% View()



currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+
beta_blocker_non_selective+fev1fvcpct_+SGRQ_symptoms+oxygen"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "oxygen")

pList10= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList10[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList10 = pList10 %>%
  as.data.frame()

names(pList10) = c("variable", "pVal")
pList10 %>% arrange(as.numeric(pVal)) %>% View()




currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+
beta_blocker_non_selective+fev1fvcpct_+SGRQ_symptoms+oxygen+hosp1yr"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "hosp1yr")

pList11= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList11[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList11 = pList11 %>%
  as.data.frame()

names(pList11) = c("variable", "pVal")
pList11 %>% arrange(as.numeric(pVal)) %>% View()


currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+
beta_blocker_non_selective+fev1fvcpct_+SGRQ_symptoms+oxygen+hosp1yr+
rateExacerb"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "rateExacerb")

pList12= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList12[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList12 = pList12 %>%
  as.data.frame()

names(pList12) = c("variable", "pVal")
pList12 %>% arrange(as.numeric(pVal)) %>% View()





currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+
beta_blocker_non_selective+fev1fvcpct_+SGRQ_symptoms+oxygen+hosp1yr+
rateExacerb+icslama"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "icslama")

pList13= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList13[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList13 = pList13 %>%
  as.data.frame()

names(pList13) = c("variable", "pVal")
pList13 %>% arrange(as.numeric(pVal)) %>% View()



train.data$panic1 = if_else(train.data$sudden_feelings_panic <4, 0, 1)
currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+
beta_blocker_non_selective+fev1fvcpct_+SGRQ_symptoms+oxygen+hosp1yr+
rateExacerb+icslama+icslaba"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "icslaba")

pList14= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList14[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList14 = pList14 %>%
  as.data.frame()

names(pList14) = c("variable", "pVal")
pList14 %>% arrange(as.numeric(pVal)) %>% View()




currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+
beta_blocker_non_selective+fev1fvcpct_+SGRQ_symptoms+oxygen+hosp1yr+
rateExacerb+icslama+icslaba+beta_blocker_use"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
updated.vars = setdiff(updated.vars, "beta_blocker_use")

pList15= matrix(data = NA, nrow = length(updated.vars), ncol = 2)

for(i in 1:length(updated.vars)){
  var.name = updated.vars[i]
  print(var.name)
  currentFormula = as.formula (paste0(currString, "+", var.name))
  curr.geebin <- geese(formula = currentFormula,
                       data=train.data,scale.fix=TRUE,family=gaussian,
                       id=train.data$sub.id,
                       mean.link="logit",corstr="unstructured",jack=TRUE)
  pList15[i , ] =c(var.name, summary(curr.geebin)$mean[[var.name, "p"]])
}

pList15 = pList15 %>%
  as.data.frame()

names(pList15) = c("variable", "pVal")
pList15 %>% arrange(as.numeric(pVal)) %>% View()



currString = "pseudoEst~oneMonth+lastYear+bronchitis+gender+trtgroup+
beta_blocker_non_selective+fev1fvcpct_+SGRQ_symptoms+oxygen+hosp1yr+
rateExacerb+icslama+icslaba+beta_blocker_use"
formula = as.formula(currString)


curr.geebin <- geese(formula = formula,
                     data=train.data,scale.fix=TRUE,family=gaussian, 
                     id=train.data$sub.id,
                     mean.link="logit",corstr="unstructured",jack=TRUE)
summary(curr.geebin)$mean
cov.test.gee=as.matrix(cbind(rep(1,nrow(test.data)),
                             test.data[, c("oneMonth",
                                           "lastYear",
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
# ask about htis.... because if we did the full formula, we wouldn't actually
# keep them all

betaMat=matrix(curr.geebin$beta)
# need to get into how to fit by hand...
currPred = 1/(1+exp(-cov.test.gee %*% betaMat))

c_censored(predicted = currPred, calTime = test.data$timeToEvent,
           observed = test.data$delta)




####### OK THE FINAL MODEL IS #########
###
# currString = "pseudoEst~timeSince+anxiety_HADS+sleep1+nowsmk"
#########


test.data$sleep1 = if_else(test.data$sleep_dysfunction >= 1, 1, 0)
cov.test.gee= cbind(rep(1, nrow(test.data)),
  test.data[,c("timeSince", "anxiety_HADS",
               "sleep1", "nowsmk")])

dim(cov.test.gee)
cov.test.gee=matrix(as.numeric(unlist(cov.test.gee)), nrow= nrow(test.data))

betaMat=matrix(curr.geebin$beta)
naivePred = 1/(1+exp(-cov.test.gee %*% betaMat))
# sepPredList[["naive"]]<-naivePred
head(naivePred)

# c_stat(predicted = naivePred, actual = test.data$pseudoEst)
c_censored(predicted = naivePred, calTime = test.data$timeToEvent, observed = test.data$delta)
# 0.4011788

cov.train.gee= cbind(rep(1, nrow(train.data)),
                    train.data[,c("timeSince", "anxiety_HADS",
                                  "sleep1", "nowsmk")])

dim(cov.train.gee)
cov.train.gee=matrix(as.numeric(unlist(cov.train.gee)), nrow= nrow(train.data))

betaMat=matrix(curr.geebin$beta)
naivePred = exp(cov.train.gee %*% betaMat)/(1+exp(cov.train.gee %*% betaMat))
# sepPredList[["naive"]]<-naivePred
head(naivePred)

# c_stat(predicted = naivePred, actual = test.data$pseudoEst)
c_censored(predicted = naivePred, calTime = train.data$timeToEvent, observed = train.data$delta)
# 0.387121

################### FORWARD SELECTION LOOPS - C-statistic ###################

permissible.names = setdiff(c(selnHRF), 
                            c("checkIn", "timeSince", "marker", "visit", "X" , 
                              "age10"))


chooseNewCModel= function(numVars){
  n = nrow(train.data)
  foo = t(combn(permissible.names , numVars))
  resFrame = matrix(data = NA, nrow = nrow(foo), ncol = 2)
  
  for(i in 1:nrow(foo)){
    myString = "pseudoEst~"
    # get the varaibles and as a formula
    for(j in 1:numVars){
      myString = paste0(myString, foo[i, j], "+")
    }
    myFormula = as.formula(str_sub(myString, start = 1, end = -2))
    print(myFormula)
    tryCatch({
      curr.geebin <- geese(formula = myFormula,
                           data = train.data, scale.fix = TRUE, family = gaussian,
                           id = train.data$sub.id,
                           mean.link = "logit", corstr = "unstructured", jack = TRUE)
      cov.test.gee = as.matrix(cbind(rep(1, nrow(train.data)), train.data[, foo[i, ]]))
      betaMat = matrix(curr.geebin$beta)
      currPred = 1 / (1 + exp(-cov.test.gee %*% betaMat))
      new.c = rcorr.cens(currPred, Surv(train.data$timeToEvent, train.data$delta))
      resFrame[i, ] = c(str_sub(myString, start = 1, end = -2), new.c[["C Index"]])
    }, error = function(e) {
      # Handle the error here, you can print an error message or perform other actions
      cat("Error occurred: ", conditionMessage(e), "\n")
      # Continue to the next iteration
    })
  }
  
  return(resFrame)
 
}





#####


hrf.z.df = viDF %>% filter(Z.value>=1.96) %>% 
  dplyr::select(Predictor, Z.value) %>% 
  mutate(model.appears = "hrf")

rownames(wald.z.df)= rownames(c.z.df)=rownames(summer.z.df) = NULL
names(wald.z.df)==names(summer.z.df)
names(hrf.z.df)=names(summer.z.df)

z.graphic.df = rbind(summer.z.df, hrf.z.df, wald.z.df, c.z.df)

write.csv(z.graphic.df,"z.graphic.df.csv")

  

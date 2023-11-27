#########################################################
# for analyzing the simulation performance
#########################################################

pacman::p_load(dplyr, tidyr, pseudo, geepack, survival,
                htree, ggplot2, tidymodels, parsnip,
                rsample, xtable, ranger, ggsurvfit,
                parallelly, Hmisc)


#setwd("~/Dropbox (University of Michigan)/LLZ/Simulation/broken_LM")
Sys.setenv(R_PARALLEL_MAKENODEPSOCK_SETUP_STRATEGY = "sequential")
Sys.setenv(R_PARALLEL_MAKENODEPSOCK_TRIES=10)
Sys.setenv(R_PARALLEL_RANDOM_PORTS = "10000:39999")
Sys.setenv(R_PARALLELLY_MAKENODEPSOCK_SETUP_STRATEGY = "sequential")
Sys.setenv(R_PARALLELLY_MAKENODEPSOCK_TRIES=10)
Sys.setenv(R_PARALLELLY_RANDOM_PORTS = "10000:39999")

source("Generate_Data_Setting_2.R")
source("internalFxns.R")
##### Set-Up ####
predC = as.data.frame(matrix(NA, nrow = 1, ncol = 9))

settings = NULL
# this is for cluster environment set up.
alpha = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
p.censoring = as.numeric(Sys.getenv("p_censoring"))
correl = as.numeric(Sys.getenv("correl"))

print(alpha)
print(paste0("correl: ", correl))
print(paste0("censoring: ", p.censoring))
#for(alpha in 1:numTimes){
# print(paste("time: ", alpha))

predList=finalPredictions = NULL #reinitialize predictions after each run.

study_length.set = 2
space.set = round(1/12, 3)
window.set = round(1/6, 3)
num_covs.set = 7
study_size.set = 500
estimand.set = 2
discard_excess.set = T
  
seln<-c(paste0("X_", 1:num_covs.set), "marker")
  
data_package = gen.data(study_length = study_length.set, 
                          space = space.set, 
                          window = window.set, 
                          num_covs = num_covs.set, 
                          study_size = study_size.set, 
                          p.censoring = p.censoring, 
                          correl = correl, 
                          estimand =estimand.set, 
                          discard_excess = discard_excess.set)
  
  data = data_package$data
  generatingLambdas = data_package$generatingLambdas
  names(data)
  
  data$sub.id = data$ID
  names(data)
  data$sub.id = data$ID
  
  
  #GEE glm
  #dat = cbind.data.frame(X1, X2, X3, X4, X5, X6, X7,  
#  exp(X2*sin(X1/X6))+ ifelse(X2>2, -X3, X3)+ 3*X1*X6 + X2^2*X4)
  # get variables into the correct form
  data$fooVar = ifelse(data$X_2>2, -data$X_3, data$X_3)
  data$X2sinX1 = data$X_2* sin(data$X_1/data$X_6)
  data$inter = data$X_1*data$X_6
  data$complicated = data$X_2*data$X_2*data$X_4
  
#   dd.gee <- cbind.data.frame(y=data$pseudoEst, data[, c("X2sinX1", "fooVar", "inter", "complicated", "avg_exacerb_impute")])
#   fit.geebin <- geese(y ~ X2sinX1 + fooVar +inter + complicated+ avg_exacerb_impute, 
#                       data=dd.gee,scale.fix=TRUE,family=gaussian, id=data$sub.id,
#                       mean.link="logit",corstr="unstructured",jack=TRUE)
#   # no way to predict out of the geese method:
#   ntest=nrow(data)
#   cov.test.gee=as.matrix(cbind(rep(1,ntest), data[, c("X2sinX1", "fooVar", "inter", "complicated", "avg_exacerb_impute")]))
#   beta=matrix(fit.geebin$beta)
#   predList[["magicImpute"]]<-1/(1+exp(-cov.test.gee %*% beta))
  
  
  #Oracle no history
  
  dd.gee <- cbind.data.frame(y=data$pseudoEst, data[, c("X2sinX1", "fooVar", "inter", "complicated")])
  fit.geebin <- geese(y ~ X2sinX1 + fooVar +inter + complicated, 
                      data=dd.gee,scale.fix=TRUE,family=gaussian, id=data$sub.id,
                      mean.link="logit",corstr="unstructured",jack=TRUE)
  # no way to predict out of the geese method:
  ntest=nrow(data)
  cov.test.gee=as.matrix(cbind(rep(1,ntest), data[, c("X2sinX1", "fooVar", "inter", "complicated")]))
  beta=matrix(fit.geebin$beta)
  predList[["magicNoHist"]]<-1/(1+exp(-cov.test.gee %*% beta))
  
  
  dd.gee <- cbind.data.frame(y=data$pseudoEst, data[, c("X2sinX1", "fooVar", "inter", "complicated", "avg_exacerb_true")])
  fit.geebin <- geese(y ~ X2sinX1 + fooVar +inter + complicated + avg_exacerb_true, 
                      data=dd.gee,scale.fix=TRUE,family=gaussian, id=data$sub.id,
                      mean.link="logit",corstr="unstructured",jack=TRUE)
  # no way to predict out of the geese method:
  ntest=nrow(data)
  cov.test.gee=as.matrix(cbind(rep(1,ntest), data[, c("X2sinX1", "fooVar", "inter", "complicated", "avg_exacerb_true")]))
  beta=matrix(fit.geebin$beta)
  predList[["magicHist"]]<-1/(1+exp(-cov.test.gee %*% beta))
  
  
  
  
  #Naive LM
  dd.naive <- cbind.data.frame(y=data$pseudoEst, data[, c(paste0("X_", 1:num_covs.set))])
  fit.geebin <- geese(y ~ . ,
                      data=dd.naive,scale.fix=TRUE,family=gaussian, id=data$sub.id,
                      mean.link="logit",corstr="unstructured",jack=TRUE)
  # no way to predict out of the geese method:
  ntest=nrow(data)
  cov.test.gee=as.matrix(cbind(rep(1,ntest), data[, c(paste0("X_", 1:num_covs.set)) ]))
  beta=matrix(fit.geebin$beta)
  predList[["naiveNoHist"]] = 1/(1+exp(-cov.test.gee %*% beta))
  
  
#     #Naive LM
#   dd.naive <- cbind.data.frame(y=data$pseudoEst, data[, c(paste0("X_", 1:num_covs.set), "avg_exacerb_impute")])
#   fit.geebin <- geese(y ~ . ,
#                       data=dd.naive,scale.fix=TRUE,family=gaussian, id=data$sub.id,
#                       mean.link="logit",corstr="unstructured",jack=TRUE)
#   # no way to predict out of the geese method:
#   ntest=nrow(data)
#   cov.test.gee=as.matrix(cbind(rep(1,ntest), data[, c(paste0("X_", 1:num_covs.set), "avg_exacerb_impute") ]))
#   beta=matrix(fit.geebin$beta)
#   predList[["naiveImpute"]] = 1/(1+exp(-cov.test.gee %*% beta))
  
  
    dd.naive <- cbind.data.frame(y=data$pseudoEst, data[, c(paste0("X_", 1:num_covs.set), "avg_exacerb_true")])
  fit.geebin <- geese(y ~ . ,
                      data=dd.naive,scale.fix=TRUE,family=gaussian, id=data$sub.id,
                      mean.link="logit",corstr="unstructured",jack=TRUE)
  # no way to predict out of the geese method:
  ntest=nrow(data)
  cov.test.gee=as.matrix(cbind(rep(1,ntest), data[, c(paste0("X_", 1:num_covs.set), "avg_exacerb_true") ]))
  beta=matrix(fit.geebin$beta)
  predList[["naiveHist"]] = 1/(1+exp(-cov.test.gee %*% beta))
  
  
#   # hrf HISTORY
#   trainDf = as.data.frame(data[,c("pseudoEst", seln, "avg_exacerb_impute")])
#   htree.rf=hrf(x=trainDf,time=data$marker,
#               id=as.numeric(data$ID),yindx="pseudoEst",ntrees=500,
#               historical=FALSE,nsamp=30, se=FALSE, control = list(nodesize = 40))
#   # this currently does not have the correct minimal number of leaves.
#   predList[["hrfImpute"]] = htree.rf$pred_oob
  
  
    trainDf = as.data.frame(data[,c("pseudoEst", seln, "avg_exacerb_true")])
  htree.rf=hrf(x=trainDf,time=data$marker,
               id=as.numeric(data$ID),yindx="pseudoEst",ntrees=500,
               historical=FALSE,nsamp=30, se=FALSE, control = list(nodesize = 40))
  # this currently does not have the correct minimal number of leaves.
  predList[["hrfHistory"]] = htree.rf$pred_oob
  
    # hrf NO HISTORY
  trainDf = as.data.frame(data[,c("pseudoEst", seln)])
  htree.rf=hrf(x=trainDf,time=data$marker,
               id=as.numeric(data$ID),yindx="pseudoEst",ntrees=500,
               historical=FALSE,nsamp=30, se=FALSE, control = list(nodesize = 40))
  # this currently does not have the correct minimal number of leaves.
  predList[["hrfNoHistory"]] = htree.rf$pred_oob
  
  
  finalPredictions = as.data.frame(predList)

rcorr.cens(currPred, Surv(data$calTime - data$marker, data$delta))


  print("calculating c...")
  #order is true, hrf, naive
my.surv.obj = Surv(data$calTime - data$marker, data$delta)
  # true:
  m1 = rcorr.cens(finalPredictions$magicHist, 
                  my.surv.obj)

m2 = rcorr.cens(finalPredictions$magicNoHist, 
                  my.surv.obj)
# here, we are concerned with mean behavior, so while one could possibly
# impute and average, we found little difference in simulation between imputations
# with m= 10, and a sinlge partial history case. code should be modified here
# to deal with the partial history case.

m3 = rcorr.cens(finalPredictions$magicImpute, 
                  my.surv.obj); m3
                    
h1 = rcorr.cens(finalPredictions$hrfHistory, 
                  my.surv.obj); h1
h2 = rcorr.cens(finalPredictions$hrfNoHistory, 
                  my.surv.obj); h2

h3 = rcorr.cens(finalPredictions$hrfImpute, 
                  my.surv.obj); h3
n1 = rcorr.cens(finalPredictions$naiveHist, 
                  my.surv.obj); n1

n2 = rcorr.cens(finalPredictions$naiveNoHist, 
                  my.surv.obj); n2
 
n3 = rcorr.cens(finalPredictions$naiveImpute, 
                  my.surv.obj); n3                   

  predC[1, ] = c(m1, m2, m3, h1, h2, h3, n1, n2, n3)
    print(predC)
  


resultDataFrame = NULL
resultDataFrame = cbind.data.frame(predC)
names(resultDataFrame) = c("magicHist", "magicNoHist","magicImpute", 
                           "hrfHist", "hrfNoHist", "hrfImpute",
                           "naiveHist", "naiveNoHist", "naiveImpute")

r_data_path=paste0("r_data/", p.censoring ,"/" , correl, "/")
if(!file.exists(r_data_path))
  dir.create(file.path(r_data_path), recursive = TRUE)
save(resultDataFrame, file=paste0(r_data_path, "/resultDataFrame_", alpha, ".Rdata"))

print("censoring:")
print(p.censoring)
print("correlation:")
print(correl)




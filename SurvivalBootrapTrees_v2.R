###Survival bootstrap v2
#rpart
rm(list=ls())
library(readxl)
library(openxlsx)
library(rpart)
library(MCMCpack)
library(FNN)
library(reshape2)
library(ggplot2)
library(pROC)
library(boot)
library(simpleboot)
library(caret)
library(resample)
library(evtree)
library(dbarts)
library(partykit)
library(survival)
library(survminer)
library(dplyr)
library(caret)
library(tictoc)
library(treeClust)
library(simsurv)
#library(ggfortify)
library(rpart)
library(randomForestSRC)
library(pec)
library(ROCR)
library(survivalROC)
library(timeROC)
library(tgp)
library(extraDistr)
#library(PCAtools)
library(dplyr)
library(coxed)
library(ipred)
library(SurvMetrics)
source("./Survival bootstrap/getpartition.R")
source("./Survival bootstrap/SurvBootTrees_cat_function.R")
source("./Survival bootstrap/SurvBootTrees_aggregating_function.R")
set.seed(123)
N<-50
n_var<-5
#nr_fold=5
n_sim<-100

# summary(simdata[[1]]$data)
# data<-simdata[[1]]$data
# model <- coxph(Surv(y, failed) ~ ., data=data)
# model$coefficients ## model-estimated coefficients
# simdata[[1]]$betas
int_brier_fin<-matrix(nrow=n_sim,ncol=3)
int_brier_rf_fin<-vector()
int_brier_cox_fin<-vector()
cat_val<-c("A","B","C")
for(sim in 1:n_sim){
  #sim<-1
  tic()
  rcat(N,c(0.5,0.5))
  X<-data.frame(X1=rcat(N,c(0.3,0.7)),X2=rcat(N,c(0.5,0.5)),X3=rcat(N,c(0.33,0.34,0.33)),X4=rcat(N,c(0.25,0.75)),X5=rcat(N,c(0.4,0.6)))
  #X<-data.frame(X1=runif(N,0,1),X2=rcat(N,c(0.5,0.5)),X3=rcat(N,c(0.33,0.34,0.33)),X4=rcat(N,c(0.25,0.75)),X5=runif(N,0,1))
 # X<-data.frame(X1=runif(N,0,1),X2=runif(N,0,1),X3=rcat(N,c(0.33,0.34,0.33)),X4=runif(N,0,1),X5=runif(N,0,1))
  
  #Here we can sample from the multivariate distribution
  #Estimate the covariance matrix of the simulated data
 # X<-data.frame(X1=runif(N,0,1),X2=runif(N,0,1),X3=runif(N,0,1),X4=runif(N,0,1),X5=runif(N,0,1))
  # simdata <- sim.survdata(N, T=100, xvars=n_var, censor=.20, num.data.frames = 1,X=X)
  #simdata<-simsurv(dist="weibull",lambdas=0.1,gammas=0.5,x=X,betas=c(X1=-1,X2=0.5,X3=2,X4=-0.5,X5=0.1))
  simdata<-simsurv(dist="weibull",lambdas=0.1,gammas=0.5,x=X,betas=c(X1=1,X2=1,X3=1,X4=1,X5=1))
  tmax<-quantile(simdata$eventtime,0.9)
  ind_max<-which(simdata$eventtime>tmax)
  simdata$eventtime[ind_max]<-tmax
  simdata$status[ind_max]<-0
  #print(table(simdata$status)/N)
  print(sim)
  
  #data<-simdata$data
  data<-X
  data$X1<-as.factor(cat_val[X$X1])
  data$X2<-as.factor(cat_val[X$X2])
  data$X3<-as.factor(cat_val[X$X3])
  data$X4<-as.factor(cat_val[X$X4])
  data$X5<-as.factor(cat_val[X$X5])
  data$time<-simdata$eventtime
  data$event<-simdata$status
  #names(data)[c(n_var+1,n_var+2)]<-c("time","event")
  rm(X,simdata)
  print(summary(coxph(Surv(time,event)~X1+X2+X3+X4+X5, data =data, x=TRUE)))
  
  
  # bayesian dirichlet with k = 0 
  #start_time <- Sys.time()
  trainind<-createDataPartition(data$event, p
                                = .75,
                                list = FALSE,
                                 times = 1)
  # trainind<-createDataPartition(data$time, p 
  #                               = .75, 
  #                               list = FALSE, 
  #                               times = 1)
  train <- data[trainind, ]
  test <- data[-trainind, ]
  # flds <- createFolds(data$time, k = nr_fold, list = TRUE, returnTrain = FALSE)
  #prior <- "uniform"
  #pred_val_bayes0_tot <- data.frame()
  
  #test_bayes0_auc <- numeric()
  #test_bayes0_mis_err <- numeric()
  
  #int_brier_boot<-matrix(nrow=nr_fold,ncol=3)
  int_brier_boot<-numeric()
  int_brier_rf<-numeric()
  int_brier_cox<-numeric()
  #for(iter in 1:nr_fold){
  #print(iter)
  #iter<-1
  # train <- data[-flds[[iter]], ]
  # test <- data[flds[[iter]], ]
  n <- nrow(train)
  m<-n
  k_all<-c(0,n*0.25/0.75,n)
  for (k_ind in 1:3){
    #k_ind<-1
    k<-k_all[k_ind]
    B <- 100
    weights <- rdirichlet(B, rep((n+k)/m,m))
    
    #oob_bayes0 <- numeric()
    #test_bayes0 <- numeric()
    #test_bayes0_all <- numeric()
    #pred_test_all <- data.frame()
    #pred_test_all <- matrix(data=0, nrow=nrow(test), ncol=length(unique(data$time)))
    #pred_val_bayes0 <- data.frame()
    #pred_all<-list()
    #mod_int <- vector(mode="list", length=B)
 
  trees_mod<-SurvBootTrees(B,k,train,weight_par=weights,sampling="empirical")
  KM_agg<-SBTaggregating(trees_mod,test,n_var,N)
    
    surv_obj<-Surv(test$time,test$event)
    int_brier_boot[k_ind]<-ipred::sbrier(surv_obj,KM_agg)
    #predict(survfit(surv_obj~X1+X2+X3+X4+X5,data=test))
    #csurv(test,RET)
    # bag_pred[[iter]]<-RET
    
  }
  # int_brier_fin<-rbind(int_brier_fin,int_brier_boot)
  # }
  #int_brier_fin[sim,]<-colMeans(int_brier_boot)
  int_brier_fin[sim,]<-int_brier_boot
  # trainind<-createDataPartition(data$event, p = .75, 
  #                               list = FALSE, 
  #                               times = 1)
  #for(iter in 1:nr_fold){
  #print(iter)
  # train <- data[-flds[[iter]], ]
  # # test <- data[flds[[iter]], ]
  # train <- data[trainind, ]
  # test <- data[-trainind, ]
  # n <- nrow(train)
  # 
  #surv_object <- Surv(time = train$surv_time, event = train$status)
  random_f<-rfsrc(Surv(time,event)~X1+X2+X3+X4+X5, data =train,
                  importance = "permute",verbose = TRUE, probability = TRUE)# splitrule = "extratrees" #### 
  dis_time = random_f$time.interest
  mat_rsf = predict(random_f, test)$survival
  cox_mod<-coxph(Surv(time,event)~X1+X2+X3+X4+X5, data =train, x=TRUE)
  mat_cox = predictSurvProb(cox_mod, test, dis_time)
  # pred<-predict(random_f,test,outcome="test",perf.type="brier")
  #pred<-surv_fit(random_f,test)
  # metrics_rsf[i] = Cindex(surv_obj, predicted = mat_rsf[, med_index])
  surv_obj<-Surv(test$time,test$event)
  #get.brier.survival(random_f,test,"rfsrc")
  int_brier_rf<-c(int_brier_rf,IBS(surv_obj, sp_matrix = mat_rsf, dis_time))
  int_brier_cox<- c(int_brier_cox,IBS(surv_obj, sp_matrix = mat_cox, dis_time))
  
  #}
  # int_brier_rf_fin[sim]<-mean(int_brier_rf)
  # int_brier_cox_fin[sim]<-mean(int_brier_cox)
  int_brier_rf_fin[sim]<-int_brier_rf
  int_brier_cox_fin[sim]<-int_brier_cox
  toc()
}



# int_brier_final<-cbind(colMeans(int_brier_fin),rbind(quantile(int_brier_fin[,1],c(0.05,0.95)),quantile(int_brier_fin[,2],c(0.05,0.95)),quantile(int_brier_fin[,3],c(0.05,0.95))))
# # int_brier_rf_final<-c(mean(int_brier_rf_fin[-which(int_brier_rf_fin>1)]),quantile(int_brier_rf_fin[-which(int_brier_rf_fin>1)],c(0.05,0.95)))
# # int_brier_cox_final<-c(mean(int_brier_cox_fin[-which(int_brier_cox_fin>1)]),quantile(int_brier_cox_fin[-which(int_brier_cox_fin>1)],c(0.05,0.95)))
int_brier_final<-cbind(colMeans(int_brier_fin),rbind(quantile(int_brier_fin[,1],c(0.05,0.95)),quantile(int_brier_fin[,2],c(0.05,0.95)),quantile(int_brier_fin[,3],c(0.05,0.95))))
int_brier_rf_final<-c(mean(int_brier_rf_fin),quantile(int_brier_rf_fin,c(0.05,0.95)))
int_brier_cox_final<-c(mean(int_brier_cox_fin),quantile(int_brier_cox_fin,c(0.05,0.95)))



avg_IBS<-c(int_brier_final[,1],int_brier_rf_final[1],int_brier_cox_final[1])
L<-c(int_brier_final[,2],int_brier_rf_final[2],int_brier_cox_final[2])
U<-c(int_brier_final[,3],int_brier_rf_final[3],int_brier_cox_final[3])

model<-c("GBEST w=0","GBEST w=0.25","GBEST w=0.5","Random forest","Cox model")
Perf_df<-data.frame(model,avg_IBS,L,U)
save(Perf_df,file="./Res_rpart/res_N50_cens10_comparison_5cat_emp.RData")
#load("./Res_rpart/Results_categorical/res_N50_cens10_comparison_5cat.RData")
plot_res<-ggplot(Perf_df,aes(x=model,y=avg_IBS))+geom_point()+geom_errorbar(aes(ymin=L,ymax=U))#+ylim(0.4,1)
ggplot(Perf_df,aes(x=model,y=avg_IBS))+geom_point()+geom_errorbar(aes(ymin=L,ymax=U))+ylim(0,1)


raw_results<-data.frame(cbind(int_brier_fin,int_brier_rf_fin,int_brier_cox_fin))
names(raw_results)<-model
save(raw_results,file="./Res_rpart/Raw_res_N50_cens10_comparison_1cat_emp.RData")
raw_results_long<-data.frame(IBS=c(raw_results$`GBEST w=0`,raw_results$`GBEST w=0.25`,raw_results$`GBEST w=0.5`,raw_results$`Cox model`,raw_results$`Random forest`),sampling=rep(model,each=nrow(raw_results)))

ggplot(raw_results_long,aes(x=sampling,y=IBS))+geom_boxplot()


library(ggplot2)

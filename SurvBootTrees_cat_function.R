

SurvBootTrees<-function(B,k,train,weight_par,sampling){
  #weight_par<-weights
  #sampling="empirical"
  models<-list()
  datasets<-list()
  index_node<-list()
  n <- nrow(train)
  m<-n
  #k<-n #Solo ora, da cancellare nella funzione
 # weights <- rdirichlet(B, rep((n+k)/m,m))
  mod_int <- vector(mode="list", length=B)
  for(i in 1:B){
    #i<-1
    test_par<-1
    while(test_par==1){
    #Genero il dataset su cui allenare l'albero
    xstar <- data.frame()
    index <- numeric()
    for(j in 1:m){ 
      soglia <- runif(1,0,1)
      if(soglia>(k/(k+n))){
        ind <- sample.int(n,size=1)
        xstar <- rbind(xstar,train[ind,])
        index <- c(index,ind)
      } else{
       # x_star <- vector()
      
        x_star<-list()
        #################### Loop that generates the new values of covariates, column by column
        for(l in 1:(ncol(train)-2)){
          if (class(train[,l])=="numeric"){
            x_star[[l]] <- runif(1,min(train[,l], na.rm=TRUE),max(train[,l],na.rm=TRUE))
          } else {
            if (sampling=="empirical"){
              tab<-table(train[,l])
              x_star[[l]] <- sample(names(tab),1,prob=tab/n) #con i pesi del dataset_ giusto usare il train 
            } else if(sampling=="uniform"){
              tab<-table(train[,l])
              x_star[[l]] <- sample(names(tab),1,prob=rep(1/length(names(tab)),length(names(tab))))
            }
            
          }
        }
        #Replace it with a multimensional samplig, giving the estimated covariance matrix
        dfx_star <- data.frame(x_star)
        colnames(dfx_star) <- paste0("X", seq(1:(ncol(dfx_star))))
        dfx_star<-dfx_star %>% mutate_if(is.character,as.factor)
        #y_star <-  knn(train[-c(which(names(train)=='y'),which(names(train)=='censored'))],cl=train$time,test=dfx_star, k = 5)[1]
       mod<-survreg(Surv(time,event)~.,data=train,dist="exponential")#"lognormal")
       time_star <- predict(mod,newdata=dfx_star,type="response")
        xstar <- rbind(xstar,data.frame(dfx_star,time=time_star,event=1))
        ####Introdurre la censura?
        }
        }
   
    tab<-data.frame(table(xstar$event))
    if(length(tab$Freq[tab$Var1==1])>0){
      test_par<-0
    }
    }
      
   
    if(sum(xstar$time==0)>0) xstar$time[xstar$time==0]<-0.5
    
    tree<-rpart(Surv(time,event) ~ .,data = xstar,  weight_par[i,],method ="exp")
    #summary(tree)
    #library(rpart.plot)
    #rpart.plot(tree)
    #tree$where
    #tree<-bagging(Surv(time,event) ~ .,data = xstar,  weights = weights[i,],nbagg=1)
    # tab <- tabulate(bindx, nbins = nobs)
    #mc$data <- m[bindx,,drop = FALSE] ### tab * weights
    #this <- irpart(Surv(xstar$time,xstar$event) ~ X1+X2+X3+X4+X5, data=xstar, control=rpart.control(xval=0), 
                   # bcontrol=list(nbagg=1, ns=dim(train)[1]))
                   # this <- irpart(Surv(xstar$time,xstar$event) ~ ., data=xstar, control=rpart.control(xval=0), 
                  # bcontrol=list(nbagg=1, ns=dim(train)[1]))
    models[[i]] <- tree
    datasets[[i]]<-xstar
    index_node[[i]]<-tree$where
  }
  return(list(models=models,datasets=datasets,index_node=index_node))
}

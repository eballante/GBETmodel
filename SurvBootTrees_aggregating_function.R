# 
# newdata<-test
# OOB=FALSE
# comb=FALSE
SBTaggregating<-function(trees_mod,newdata,n_var,N,OOB=FALSE,comb=FALSE) {
  

object<-list()
object$comb<-comb
object$OOB<-OOB
object$mtrees<-trees_mod$models
#object$y<-Surv(train$time,train$event)
#object$X<-train[,-c(n_var+1,n_var+2)]
class(object)<-"survbagg"
#newdata=test
if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
NN <- nrow(newdata)
if (!object$comb) {
  tree_mod <- object$mtrees[[1]]
  Terms <- delete.response(trees_mod$models[[1]]$terms)
  act <- (tree_mod$call)$na.action
  if (is.null(act)) act<- na.rpart
  newdata <- model.frame(Terms, newdata, na.action = act,
                         xlev=attr(tree_mod, "xlevels"))
  #newdata <- getFromNamespace("rpart.matrix", ns = "rpart")(newdata) ##
}
agglsample <- list()#creo le liste delle dimensioni giuste
aggcens <- list()
for (j in 1:N) { 
  agglsample <- c(agglsample, list(c()))
  aggcens <- c(aggcens, list(c()))
}
for (i in 1:length(object$mtrees)) {
  bdata <- Surv(trees_mod$datasets[[i]]$time,trees_mod$datasets[[i]]$event)
  #newpart <- getFromNamespace("pred.rpart", ns = "rpart")(trees_mod$models[[i]], newdata)
  newpart<-rpart.predict.leaves(trees_mod$models[[i]],newdata,type = "where")
  oldpart <- trees_mod$index_node[[i]]
  # if (OOB) {
  #   if (!is.null(object$mtrees[[i]][[1]]$bfct))
  #     stop("cannot compute out-of-bag estimate for combined models!")
  #   tindx <- (1:N)[-object$mtrees[[i]][[1]]$bindx]
  # } else {
  tindx <- 1:NN #Per ogni osservazione nei dati di test 
  # }
  for (j in tindx) {
    aggobs <- bdata[oldpart == newpart[j],1] 
    
    agglsample[[j]] <- c(agglsample[[j]], aggobs)
    aggobs <- bdata[oldpart == newpart[j],2]
    aggcens[[j]] <- c(aggcens[[j]], aggobs)
  }
}
RET <- list()
for (j in 1:NN) RET <- c(RET, list(survfit(Surv(agglsample[[j]], aggcens[[j]]) ~ 1)))


return(RET)

}
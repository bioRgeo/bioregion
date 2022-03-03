controls=function(X, weight = TRUE, name = "net", type = c("data.frame, matrix, dist")){

  # data.frame
  if(type == "date.frame"){

    if(!is.data.frame(X)){
      stop(paste0(name, " must be a two- or three-columns data.frame"))
    }

    if(dim(X)[2]!=2 & dim(X)[2]!=3){
      stop(paste0(name, " must be a two- or three-columns data.frame"))
    }

    sco=sum(is.na(X))
    if(sco>0){
      stop(paste0("NA(s) detected in ", name))
    }

    if(weight & dim(X)[2]==2){
      stop(paste0(name, " must be a three-columns data.frame if weight equal TRUE"))
    }

    if(weight & dim(X)[2]==3){
      if(class(X[,3])!="numeric" & class(X[,3])!="integer"){
        stop(paste0("The third column of ", name," must be numeric"))
      }
    }

  }

}

knbclu=function(partitions, method = "max"){

  # Number of partitions
  nb = dim(partitions)[2] - 1

  # Number of clusters per partitions
  if(nb==1){
    if(method=="length"){
      colnames(partitions)[2]=paste0("K_", length(table(partitions[,2])))
    }
    if(method=="max"){
      colnames(partitions)[2]=paste0("K_", max(partitions[,2]))
    }
  }else{
    if(method == "length"){
      nbclus=NULL
      for(i in 1:nb){
        nbclus=c(nbclus,length(table(partitions[,i+1])))
      }
    }
    if(method=="max"){
      nbclus=apply(partitions[,-1], 2, max)
    }

    ord=cbind(2:(nb+1),nbclus)
    ord=as.numeric(ord[order(ord[,2]),1])
    partitions=partitions[,c(1,ord)]
    colnames(partitions)[2:(nb+1)]=paste0("K_",nbclus)
  }

  partitions

}




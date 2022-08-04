# controls=function(X, weight = TRUE, name = "net", type = c("data.frame, matrix, dist")){
#
#   # data.frame
#   if(type == "date.frame"){
#
#     if(!is.data.frame(X)){
#       stop(paste0(name, " must be a two- or three-columns data.frame"))
#     }
#
#     if(dim(X)[2]!=2 & dim(X)[2]!=3){
#       stop(paste0(name, " must be a two- or three-columns data.frame"))
#     }
#
#     sco=sum(is.na(X))
#     if(sco>0){
#       stop(paste0("NA(s) detected in ", name))
#     }
#
#     if(weight & dim(X)[2]==2){
#       stop(paste0(name, " must be a three-columns data.frame if weight equal TRUE"))
#     }
#
#     if(weight & dim(X)[2]==3){
#       if(class(X[,3])!="numeric" & class(X[,3])!="integer"){
#         stop(paste0("The third column of ", name," must be numeric"))
#       }
#     }
#
#   }
#
# }

knbclu <- function(partitions, method = "length",
                   reorder = TRUE, rename_duplicates = TRUE) {

  # Identify the number of clusters per partition
  nb <- dim(partitions)[2] - 1


  if (method == "max") {
    nbclus <- as.numeric(apply(
      partitions[, 2:(nb + 1), drop = FALSE],
      2,
      function(x) max(x)
    ))
  } else if (method == "length") {
    nbclus <- apply(
      partitions[, 2:(nb + 1), drop = FALSE],
      2,
      function(x) length(unique(x))
    )
  }

  # Rename and reorder
  if (reorder) {
    ord <- cbind(2:(nb + 1), nbclus)
    ord <- ord[order(ord[, 2]), , drop = FALSE]
    partitions <- partitions[, c(1, ord[, 1])]
    colnames(partitions)[2:(nb + 1)] <- paste0("K_", ord[, 2])
  } else {
    colnames(partitions)[2:(nb + 1)] <- paste0("K_", nbclus)
  }

  # Rename duplicates
  if (rename_duplicates) {
    colnames(partitions)[2:(nb + 1)] <- make.unique.2(colnames(partitions)[2:(nb + 1)], sep = "_")
  }

  # Convert in character
  for (k in 1:(nb + 1)) {
    partitions[, k] <- as.character(partitions[, k])
  }

  partitions
}

make.unique.2 <- function(x, sep = ".") { # From https://stackoverflow.com/questions/7659891/r-make-unique-starting-in-1
  stats::ave(x, x, FUN = function(a) {
    if (length(a) > 1) {
      paste(a, 1:length(a), sep = sep)
    } else {
      a
    }
  })
}

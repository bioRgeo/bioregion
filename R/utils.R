reformat_hierarchy <- function(input, integerize = FALSE){
  
  input=as.character(as.vector(as.matrix(input)))
  
  print(input)
  
  nblev <- max(lengths(regmatches(input, gregexpr("\\.", input)))) + 1
  
  print(nblev)
  
  output <- tidyr::separate(data = data.frame(input),
                            col = input,
                            remove = FALSE,
                            into = paste0("lvl", 1:nblev),
                            sep = "\\.",
                            fill = "right")

  output[which(is.na(output), arr.ind = TRUE)] <- 0
  
  for(lvl in grep("lvl", colnames(output))[2:nblev]) {
    output[, lvl] <- paste(output[, lvl - 1], output[, lvl],sep = ".")
  }
  
  print(output)
  
  output[grep("lvl", colnames(output))] <- lapply(output[grep("lvl", colnames(output))],
                                                  function(x) gsub("\\.0", "", x))

  if(integerize){
    for(k in 2:(nblev+1)){
      output[,k]=as.numeric(as.factor(output[,k]))
    }
  }  

  return(output)

}

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
  
  # Change colnames 1 en ID
  colnames(partitions)[1]="ID"

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

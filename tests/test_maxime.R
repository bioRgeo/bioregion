# Remove bioRgeo
unlink("/home/maxime/Applications/R/bioRgeo", recursive=TRUE)

# bioRgeo
devtools::install_github("bioRgeo/bioRgeo", build_vignettes = TRUE)
vignette("bioRgeo", package = "bioRgeo")

# Import packages
library(sf)
library(bioRgeo)

# Working directory
setwd("/home/maxime/mmmycloud/Research/Articles/InProgress/bioRgeo")

# Import data
data(vegedf)
data(vegemat)
data(vegesp)

tab=vegedf
mat=vegemat
G=vegesp

############################################################################################
# bin
#source("bioRgeo/R/bin.R")

bin()


# check results
net=spproject(mat,metric=c("Simpson"))

coml=louvain(net[net[,3]>0.5,], weight=TRUE, q=0, lang="Cpp")
comi=infomap(net[net[,3]>0.5,], weight=TRUE, markovtime=1)
como=oslom(net[net[,3]>0.5,], r=1, reassign="simil")

coml[1:10,]
comi[1:10,]
como[1:10,]

dim(coml)
dim(comi)
dim(como)

table(coml[,2])
table(comi[,2])
table(como[,2])

sp=G[match(coml[,1],G$ID),]
sp=cbind(sp,como[,2],comi[,2],coml[,2])

plot(sp)

# louvain
#source("bioRgeo/R/louvain.R")

net=spproject(mat,metric=c("Simpson"))

com=louvain(net, weight=TRUE, lang="all")

# infomap
#source("bioRgeo/R/infomap.R")

net=spproject(mat,metric=c("Jaccard"))

com=infomap(net[net[,3]>0.5,])

# oslom
#source("bioRgeo/R/oslom.R")

net=spproject(mat,metric=c("Sorensen"))

com=oslom(net[net[,3]>0.5,], r=1, reassign="simil")

# spproject
library(Rcpp)
#library(Matrix)
sourceCpp("bioRgeo/src/abc.cpp")
source("bioRgeo/R/spproject.R")

test=spproject(mat, metric=c("Jaccard", "Bray"), formula= c("1 - (b + c) / (a + b + c)", "1 - (B + C) / (2*A + B + C)"))
test[1:10,]



test=spproject(mat,metric=c("abc","ABC"), method="loops")
test[1:10,]
test=spproject(mat,metric=c("abc","ABC"))
test[1:10,]
test=spproject(mat,metric=c("abc","ABC","Euclidean"))
test[1:10,]
test=spproject(mat,metric=c("all"))
test[1:10,]

sum(mat[1,])
sum(mat[2,])
sum(pmin(mat[1,],mat[2,]))

matp=mat
matp[matp!=0]=1
sum(matp[1,])
sum(matp[2,])
sum(pmin(matp[1,],matp[2,]))

# prodmat
test=prodmat(matp,t(matp))

# abc
#library(Rcpp)
#sourceCpp("bioRgeo/src/abc.cpp")

test=abc(mat)

# df_to_contingency
#source("bioRgeo/R/df_to_contingency.R")

mat2=df_to_contingency(tab,weight=TRUE)
mat2[1:10,1:10]
mat2=df_to_contingency(tab,weight=TRUE,squared=TRUE)
mat2[1:10,1:10]
mat2=df_to_contingency(tab,weight=FALSE)
mat2[1:10,1:10]

# contingency_to_df
#source("bioRgeo/R/contingency_to_df.R")

tab2=contingency_to_df(mat)
tab2[1:10,]
tab2=contingency_to_df(mat,weight=TRUE)
tab2[1:10,]
tab2=contingency_to_df(mat,weight=TRUE,remove_absent_objects=FALSE)
tab2[1:10,]











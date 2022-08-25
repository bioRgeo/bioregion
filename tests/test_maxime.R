# Remove bioRgeo
unlink("/home/maxime/Applications/R/bioRgeo", recursive=TRUE)

# bioRgeo
devtools::install_github("bioRgeo/bioRgeo")

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

# install_binaries
#source("bioRgeo/R/install_binaries.R")
install_binaries(binpath = NULL, infomap_version = c("2.1.0","2.6.0"))

############################################################################################

















# check results
net=similarity(mat,metric=c("Simpson"))

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







############################################################################################
# walktrap
#source("bioRgeo/R/utils.R")
#source("bioRgeo/R/netclu_walktrap.R")

net=similarity(mat,metric=c("Simpson"))

com=netclu_walktrap(net,
                  weight = TRUE,
                  #index = names(net)[3],
                  steps =4,
                  bipartite = FALSE,
                  site_col = 1,
                  species_col = 2,
                  return_node_type = "both",
                  algorithm_in_output = TRUE)
com$clusters[1:10,]

com=netclu_walktrap(tab[1:1000,],                  
                  weight = TRUE,
                  #index = names(net)[3],
                  steps = 4,
                  bipartite = TRUE,
                  site_col = 1,
                  species_col = 2,
                  return_node_type = "sites",
                  algorithm_in_output = TRUE)
com$clusters[1:10,]

# leadingeigen
#source("bioRgeo/R/utils.R")
#source("bioRgeo/R/netclu_leadingeigen.R")

net=similarity(mat,metric=c("Simpson"))

com=netclu_leadingeigen(net,
                  weight = TRUE,
                  #index = names(net)[3],
                  bipartite = FALSE,
                  site_col = 1,
                  species_col = 2,
                  return_node_type = "both",
                  algorithm_in_output = TRUE)
com$clusters[1:10,]

com=netclu_leadingeigen(tab[1:1000,],                  
                  weight = TRUE,
                  #index = names(net)[3],
                  bipartite = TRUE,
                  site_col = 1,
                  species_col = 2,
                  return_node_type = "sites",
                  algorithm_in_output = TRUE)
com$clusters[1:10,]


# labelprop
#source("bioRgeo/R/utils.R")
#source("bioRgeo/R/netclu_labelprop.R")

net=similarity(mat,metric=c("Simpson"))

com=netclu_labelprop(net,
                  weight = TRUE,
                  #index = names(net)[3],
                  bipartite = FALSE,
                  site_col = 1,
                  species_col = 2,
                  return_node_type = "both",
                  algorithm_in_output = TRUE)
com$clusters[1:10,]

com=netclu_labelprop(tab[1:1000,],                  
                  weight = TRUE,
                  #index = names(net)[3],
                  bipartite = TRUE,
                  site_col = 1,
                  species_col = 2,
                  return_node_type = "sites",
                  algorithm_in_output = TRUE)
com$clusters[1:10,]

# greedy
#source("bioRgeo/R/utils.R")
#source("bioRgeo/R/netclu_greedy.R")

net=similarity(mat,metric=c("Simpson"))

com=netclu_greedy(net,
                  weight = TRUE,
                  #index = names(net)[3],
                  bipartite = FALSE,
                  site_col = 1,
                  species_col = 2,
                  return_node_type = "both",
                  algorithm_in_output = TRUE)
com$clusters[1:10,]

com=netclu_greedy(tab[1:1000,],                  
                  weight = TRUE,
                  #index = names(net)[3],
                  bipartite = TRUE,
                  site_col = 1,
                  species_col = 2,
                  return_node_type = "sites",
                  algorithm_in_output = TRUE)
com$clusters[1:10,]

# oslom
#source("bioRgeo/R/utils.R")
#source("bioRgeo/R/netclu_oslom.R")

net=similarity(mat,metric=c("Jaccard"))

com=netclu_oslom(net[net[,3]>0.5,], 
                 weight = TRUE,
                 #index = names(net)[3],
                 reassign = "no",
                 r = 1,
                 hr = 50,
                 seed = 0,
                 t = 0.1,
                 cp = 0.5,
                 directed = FALSE,
                 bipartite = FALSE,
                 site_col = 1,
                 species_col = 2,
                 return_node_type = "both",
                 delete_temp = TRUE,
                 path_temp = "oslom_temp",
                 binpath = NULL)
com$clusters[1:10,]

com=netclu_oslom(tab[1:1000,],
                 weight = TRUE,
                 #index = names(net)[3],
                 reassign = "simil",
                 r = 1,
                 hr = 50,
                 seed = 0,
                 t = 0.1,
                 cp = 0.5,
                 directed = FALSE,
                 bipartite = TRUE,
                 site_col = 1,
                 species_col = 2,
                 return_node_type = "both",
                 delete_temp = TRUE,
                 path_temp = "oslom_temp",
                 binpath = NULL)
com$clusters[1:10,]

# louvain
#source("bioRgeo/R/utils.R")
#source("bioRgeo/R/netclu_louvain.R")

net=similarity(mat,metric=c("Jaccard"))

com=netclu_louvain(net = net[net[,3]>0.5,],
                   weight = TRUE,
                   index = names(net)[3],
                   lang = "igraph",
                   q = 0,
                   c = 0.5,
                   k = 1,
                   bipartite = FALSE,
                   site_col = 1,
                   species_col = 2,
                   return_node_type = "both",
                   delete_temp = TRUE,
                   path_temp = "louvain_temp",
                   binpath = NULL,
                   algorithm_in_output = TRUE)
com$clusters[1:10,]

com=netclu_louvain(net = tab[1:1000,], 
                   weight = TRUE,
                   #index = names(net)[3],
                   lang = "Cpp",
                   q = 0,
                   c = 0.5,
                   k = 1,
                   bipartite = TRUE,
                   site_col = 1,
                   species_col = 2,
                   return_node_type = "both",
                   delete_temp = TRUE,
                   path_temp = "louvain_temp",
                   binpath = NULL,
                   algorithm_in_output = TRUE)
com$clusters[1:10,]

# infomap
#source("bioRgeo/R/utils.R")
#source("bioRgeo/R/netclu_infomap.R")

net=similarity(mat,metric=c("Jaccard"))

com=netclu_infomap(net = net[net[,3]>0.5,],
                   weight = TRUE,
                   index = names(net)[3],
                   nbmod = 0,
                   markovtime = 1,
                   seed = 0,
                   numtrials = 1,
                   twolevel = FALSE,
                   show_hierarchy = FALSE,
                   directed = FALSE,
                   bipartite_version = FALSE,
                   bipartite = FALSE,
                   site_col = 1,
                   species_col = 2,
                   return_node_type = "both",
                   version = "2.6.0",
                   delete_temp = TRUE,
                   path_temp = "infomap_temp",
                   binpath = NULL)
com$clusters[1:10,]

com=netclu_infomap(net = tab[1:1000,],                   
                   weight = TRUE,
                   #index = names(net)[3],
                   nbmod = 0,
                   markovtime = 1,
                   seed = 0,
                   numtrials = 1,
                   twolevel = FALSE,
                   show_hierarchy = FALSE,
                   directed = FALSE,
                   bipartite_version = FALSE,
                   bipartite = TRUE,
                   site_col = 1,
                   species_col = 2,
                   return_node_type = "both",
                   version = "2.6.0",
                   delete_temp = TRUE,
                   path_temp = "infomap_temp",
                   binpath = NULL)
com$clusters[1:10,]

# beckett
#source("bioRgeo/R/utils.R")
#source("bioRgeo/R/netclu_beckett.R")

com=netclu_beckett(net = tab[1:1000,], 
                   weight = TRUE,
                   index = "Abundance",
                   site_col = 1,
                   species_col = 2,
                   return_node_type = "species",
                   forceLPA = FALSE,
                   algorithm_in_output = TRUE)
com$clusters[1:10,]
############################################################################################



############################################################################################
# similarity
#library(Rcpp)
#library(Matrix)
#sourceCpp("bioRgeo/src/abc.cpp")
#source("bioRgeo/R/similarity.R")

test=similarity(mat, metric=NULL, formula= c("1 - (b + c) / (a + b + c)", "1 - (B + C) / (2*A + B + C)"))
test[1:10,]

test=similarity(mat,metric=c("abc","ABC"), method="loops")
test[1:10,]
test=similarity(mat,metric=c("abc","ABC"))
test[1:10,]
test=similarity(mat,metric=c("abc","ABC","Euclidean"))
test[1:10,]
test=similarity(mat,metric=c("all"))
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

# net_to_mat
#source("bioRgeo/R/net_to_mat.R")

mat2=net_to_mat(tab,weight=TRUE)
mat2[1:10,1:10]
mat2=net_to_mat(tab,weight=TRUE,squared=TRUE)
mat2[1:10,1:10]
mat2=net_to_mat(tab,weight=FALSE)
mat2[1:10,1:10]

# mat_to_net
#source("bioRgeo/R/mat_to_net.R")

tab2=mat_to_net(mat)
tab2[1:10,]
tab2=mat_to_net(mat,weight=TRUE)
tab2[1:10,]
tab2=mat_to_net(mat,weight=TRUE,remove_absent_objects=FALSE)
tab2[1:10,]
############################################################################################










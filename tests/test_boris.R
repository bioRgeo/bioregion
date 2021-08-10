library(bioRgeo)


# Import data
data(VegeCBNMed)

tab=vegedf
mat=vegemat
G=vegesp


simil <- spproject(vegemat, metric = "all")

# test generic output
simil
class(simil)

# conversion
distances <- similarityToDistance(simil)
distances
class(distances)


net <- simil[, c("Site1", "Site2", "Simpson")]

coml=louvain(net[net[,3]>0.5,], weight=TRUE, q=0, lang="Cpp")
comi=infomap(net[net[,3]>0.5,], weight=TRUE, markovtime=1)
como=oslom(net[net[,3]>0.5,], r=1, reassign="simil")

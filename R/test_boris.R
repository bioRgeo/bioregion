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

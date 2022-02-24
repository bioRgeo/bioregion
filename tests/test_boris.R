library(bioRgeo)


# Import data
tab=vegedf
mat=vegemat
G=vegesp


simil <- similarity(vegemat, metric = "all")

# test generic output
simil
class(simil)

# conversion
distances <- similarity_to_distance(simil)
distances
class(distances)

# network
net <- simil[, c("Site1", "Site2", "Simpson")]

coml=louvain(net[net[,3]>0.5,], weight=TRUE, q=0, lang="Cpp")
comi=infomap(net[net[,3]>0.5,], weight=TRUE, markovtime=1)
como=oslom(net[net[,3]>0.5,], r=1, reassign="simil")

# distance-based htree
tree1 <- clustering_hierarchical(distances,
                                 n_clust = 5)
tree1
plot(tree1)
str(tree1)
tree1$clusters

# User-defined height cut
# Only one height
tree2 <- clustering_hierarchical(distances,
                                 cut_height = .05)
tree2
tree2$clusters

# Multiple heights
tree3 <- clustering_hierarchical(distances,
                                 index = "Simpson",
                                 cut_height = c(.15, .35, .4))
tree3
tree3$clusters # Mind the order of height cuts: from deep to shallow cuts

a <- cut_tree(tree3, n_clust = 5)
a <- cut_tree(tree3, n_clust = 5, find_h = FALSE)
cut_tree(tree3$final.tree, n_clust = 5)
cut_tree(tree3$final.tree, cut_height =  0.3)

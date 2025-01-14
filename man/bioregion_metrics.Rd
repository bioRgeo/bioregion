% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bioregion_metrics.R
\name{bioregion_metrics}
\alias{bioregion_metrics}
\title{Calculate contribution metrics for bioregions}
\usage{
bioregion_metrics(cluster_object, comat, map = NULL, col_bioregion = NULL)
}
\arguments{
\item{cluster_object}{A \code{bioregion.clusters} object, a \code{data.frame}, or a
list of \code{data.frame}s containing multiple partitions. At least two
partitions are required. If a list of \code{data.frame}s is provided, they must
all have the same number of rows (i.e., same items in the clustering).}

\item{comat}{A co-occurrence \code{matrix} with sites as rows and species as
columns.}

\item{map}{A spatial \verb{sf data.frame} with sites and bioregions. It is the
output of the function \code{map_bioregions}. \code{NULL} by default.}

\item{col_bioregion}{An \code{integer} specifying the column position of the
bioregion.}
}
\value{
A \code{data.frame} with 5 columns, or 6 if spatial coherence is computed.
}
\description{
This function calculates the number of sites per bioregion, as well as the
number of species these sites have, the number of endemic species, and the
proportion of endemism.
}
\details{
Endemic species are species found only in the sites belonging to one
bioregion.
}
\examples{
comat_1 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
prob = 1/1:1001), 10, 12)
rownames(comat_1) <- paste0("Site", 1:10)
colnames(comat_1) <- paste0("Species", 1:12)
comat_1 <- cbind(comat_1,
                 matrix(0, 10, 8,
                        dimnames = list(paste0("Site", 1:10),
                                        paste0("Species", 13:20))))

comat_2 <- matrix(sample(0:1000, size = 10*12, replace = TRUE,
                         prob = 1/1:1001), 10, 12)
rownames(comat_2) <- paste0("Site", 11:20)
colnames(comat_2) <- paste0("Species", 9:20)
comat_2 <- cbind(matrix(0, 10, 8,
                        dimnames = list(paste0("Site", 11:20),
                                        paste0("Species", 1:8))),
                 comat_2)

comat <- rbind(comat_1, comat_2)

dissim <- dissimilarity(comat, metric = "Simpson")
clust1 <- nhclu_kmeans(dissim, n_clust = 3, index = "Simpson")

net <- similarity(comat, metric = "Simpson")
com <- netclu_greedy(net)

bioregion_metrics(cluster_object = clust1, comat = comat) 

# Spatial coherence
vegedissim <- dissimilarity(vegemat)
hclu <- nhclu_kmeans(dissimilarity = vegedissim, n_clust = 4)
vegemap <- map_bioregions(hclu, vegesf, write_clusters = TRUE, plot = FALSE)

bioregion_metrics(cluster_object = hclu, comat = vegemat, map = vegemap,
col_bioregion = 2) 

}
\seealso{
For more details illustrated with a practical example,
see the vignette:
\url{https://biorgeo.github.io/bioregion/articles/a5_3_summary_metrics.html}.

Associated functions:
\link{site_species_metrics} \link{bioregionalization_metrics}
}
\author{
Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
Boris Leroy (\email{leroy.boris@gmail.com}) \cr
Maxime Lenormand (\email{maxime.lenormand@inrae.fr})
}

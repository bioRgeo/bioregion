% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cut_tree.R
\name{cut_tree}
\alias{cut_tree}
\title{Cut a hierarchical tree}
\usage{
cut_tree(
  tree,
  n_clust = NULL,
  cut_height = NULL,
  find_h = TRUE,
  h_max = 1,
  h_min = 0,
  dynamic_tree_cut = FALSE,
  dynamic_method = "tree",
  dynamic_minClusterSize = 5,
  dissimilarity = NULL,
  ...
)
}
\arguments{
\item{tree}{A \code{bioregion.hierar.tree} or an \code{hclust} object.}

\item{n_clust}{An \code{integer} vector or a single \code{integer} indicating the
number of clusters to be obtained from the hierarchical tree, or the output
from \code{\link[=bioregionalization_metrics]{bioregionalization_metrics()}}. This should not be used concurrently
with \code{cut_height}.}

\item{cut_height}{A \code{numeric} vector specifying the height(s) at which the
tree should be cut. This should not be used concurrently with \code{n_clust} or
\code{optim_method}.}

\item{find_h}{A \code{boolean} indicating whether the cutting height should be
determined for the requested \code{n_clust}.}

\item{h_max}{A \code{numeric} value indicating the maximum possible tree height
for determining the cutting height when \code{find_h = TRUE}.}

\item{h_min}{A \code{numeric} value specifying the minimum possible height in the
tree for determining the cutting height when \code{find_h = TRUE}.}

\item{dynamic_tree_cut}{A \code{boolean} indicating whether the dynamic tree cut
method should be used. If \code{TRUE}, \code{n_clust} and \code{cut_height} are ignored.}

\item{dynamic_method}{A \code{character} string specifying the method to be used
for dynamically cutting the tree: either \code{"tree"} (clusters searched only
within the tree) or \code{"hybrid"} (clusters searched in both the tree and the
dissimilarity matrix).}

\item{dynamic_minClusterSize}{An \code{integer} indicating the minimum cluster
size for the dynamic tree cut method (see
\link[dynamicTreeCut:cutreeDynamic]{dynamicTreeCut::cutreeDynamic()}).}

\item{dissimilarity}{Relevant only if \code{dynamic_method = "hybrid"}. Provide
the dissimilarity \code{data.frame} used to build the \code{tree}.}

\item{...}{Additional arguments passed to
\link[dynamicTreeCut:cutreeDynamic]{dynamicTreeCut::cutreeDynamic()} to
customize the dynamic tree cut method.}
}
\value{
If \code{tree} is an output from \code{\link[=hclu_hierarclust]{hclu_hierarclust()}}, the same
object is returned with updated content (i.e., \code{args} and \code{clusters}). If
\code{tree} is an \code{hclust} object, a \code{data.frame} containing the clusters is
returned.
}
\description{
This function is designed to work on a hierarchical tree and cut it
at user-selected heights. It works with outputs from either
\code{hclu_hierarclust} or \code{hclust} objects. The function allows for cutting
the tree based on the chosen number(s) of clusters or specified height(s).
Additionally, it includes a procedure to automatically determine the cutting
height for the requested number(s) of clusters.
}
\details{
The function supports two main methods for cutting the tree. First, the tree
can be cut at a uniform height (specified by \code{cut_height} or determined
automatically for the requested \code{n_clust}). Second, the dynamic tree cut
method (Langfelder et al., 2008) can be applied, which adapts to the shape
of branches in the tree, cutting at varying heights based on cluster
positions.

The dynamic tree cut method has two variants:
\itemize{
\item{The tree-based variant (\code{dynamic_method = "tree"}) uses a top-down
approach, relying solely on the tree and the order of clustered objects.}
\item{The hybrid variant (\code{dynamic_method = "hybrid"}) employs a bottom-up
approach, leveraging both the tree and the dissimilarity matrix to identify
clusters based on dissimilarity among sites. This approach is useful for
detecting outliers within clusters.}
}
}
\note{
The \code{find_h} argument is ignored if \code{dynamic_tree_cut = TRUE},
as cutting heights cannot be determined in this case.
}
\examples{
comat <- matrix(sample(0:1000, size = 500, replace = TRUE, prob = 1/1:1001),
20, 25)
rownames(comat) <- paste0("Site", 1:20)
colnames(comat) <- paste0("Species", 1:25)

simil <- similarity(comat, metric = "all")
dissimilarity <- similarity_to_dissimilarity(simil)

# User-defined number of clusters
tree1 <- hclu_hierarclust(dissimilarity,
                          n_clust = 5)
tree2 <- cut_tree(tree1, cut_height = .05)
tree3 <- cut_tree(tree1, n_clust = c(3, 5, 10))
tree4 <- cut_tree(tree1, cut_height = c(.05, .1, .15, .2, .25))
tree5 <- cut_tree(tree1, n_clust = c(3, 5, 10), find_h = FALSE)

hclust_tree <- tree2$algorithm$final.tree
clusters_2 <- cut_tree(hclust_tree, n_clust = 10)

cluster_dynamic <- cut_tree(tree1, dynamic_tree_cut = TRUE,
                            dissimilarity = dissimilarity)

}
\references{
Langfelder P, Zhang B & Horvath S (2008) Defining clusters from a
hierarchical cluster tree: the Dynamic Tree Cut package for R.
\emph{BIOINFORMATICS} 24, 719-720.
}
\seealso{
For more details illustrated with a practical example,
see the vignette:
\url{https://biorgeo.github.io/bioregion/articles/a4_1_hierarchical_clustering.html}.

Associated functions:
\link{hclu_hierarclust}
}
\author{
Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
Boris Leroy (\email{leroy.boris@gmail.com})
}

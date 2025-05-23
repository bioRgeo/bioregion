% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netclu_beckett.R
\name{netclu_beckett}
\alias{netclu_beckett}
\title{Community structure detection in weighted bipartite networks via modularity
optimization}
\usage{
netclu_beckett(
  net,
  weight = TRUE,
  cut_weight = 0,
  index = names(net)[3],
  seed = NULL,
  forceLPA = FALSE,
  site_col = 1,
  species_col = 2,
  return_node_type = "both",
  algorithm_in_output = TRUE
)
}
\arguments{
\item{net}{A \code{data.frame} representing a bipartite network with the first
two columns representing undirected links between pairs of nodes, and the
next column(s) representing the weights of the links.}

\item{weight}{A \code{boolean} indicating whether weights should be considered
if there are more than two columns (see Note).}

\item{cut_weight}{A minimal weight value. If \code{weight} is TRUE, links
with weights strictly lower than this value will not be considered
(\code{0} by default).}

\item{index}{The name or number of the column to use as weight. By default,
the third column name of \code{net} is used.}

\item{seed}{The seed for the random number generator (\code{NULL} for random
by default).}

\item{forceLPA}{A \code{boolean} indicating whether the even faster pure
LPA-algorithm of Beckett should be used. DIRT-LPA (the default) is less
likely to get trapped in a local minimum but is slightly slower. Defaults
to \code{FALSE}.}

\item{site_col}{The name or number of the column for site nodes
(i.e., primary nodes).}

\item{species_col}{The name or number of the column for species nodes
(i.e., feature nodes).}

\item{return_node_type}{A \code{character} indicating which types of nodes
(\code{"site"}, \code{"species"}, or \code{"both"}) should be returned in the output
(\code{"both"} by default).}

\item{algorithm_in_output}{A \code{boolean} indicating whether the original
output of \link[bipartite:computeModules]{computeModules} should be returned
in the output (\code{TRUE} by default, see Value).}
}
\value{
A \code{list} of class \code{bioregion.clusters} with five slots:
\enumerate{
\item{\strong{name}: A \code{character} containing the name of the algorithm.}
\item{\strong{args}: A \code{list} of input arguments as provided by the user.}
\item{\strong{inputs}: A \code{list} of characteristics of the clustering process.}
\item{\strong{algorithm}: A \code{list} of all objects associated with the
clustering procedure, such as original cluster objects (only if
\code{algorithm_in_output = TRUE}).}
\item{\strong{clusters}: A \code{data.frame} containing the clustering results.}}

If \code{algorithm_in_output = TRUE}, users can find the output of
\link[bipartite:computeModules]{computeModules} in the \code{algorithm} slot.
}
\description{
This function takes a bipartite weighted graph and computes modules by
applying Newman’s modularity measure in a bipartite weighted version.
}
\details{
This function is based on the modularity optimization algorithm provided by
Stephen Beckett (Beckett, 2016) as implemented in the
\href{https://cran.r-project.org/package=bipartite}{bipartite}
package (\link[bipartite:computeModules]{computeModules}).
}
\note{
Beckett's algorithm is designed to handle weighted bipartite networks. If
\code{weight = FALSE}, a weight of 1 will be assigned to each pair of nodes.
Ensure that the \code{site_col} and \code{species_col} arguments correctly identify
the respective columns for site nodes (primary nodes) and species nodes
(feature nodes). The type of nodes returned in the output can be selected
using the \code{return_node_type} argument: \code{"both"} to include both node types,
\code{"site"} to return only site nodes, or \code{"species"} to return only species
nodes.
}
\examples{
net <- data.frame(
  Site = c(rep("A", 2), rep("B", 3), rep("C", 2)),
  Species = c("a", "b", "a", "c", "d", "b", "d"),
  Weight = c(10, 100, 1, 20, 50, 10, 20))

com <- netclu_beckett(net)

}
\references{
Beckett SJ (2016) Improved community detection in weighted bipartite
networks. \emph{Royal Society Open Science} 3, 140536.
}
\seealso{
For more details illustrated with a practical example,
see the vignette:
\url{https://biorgeo.github.io/bioregion/articles/a4_3_network_clustering.html}.

Associated functions:
\link{netclu_infomap} \link{netclu_louvain} \link{netclu_oslom}
}
\author{
Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
Boris Leroy (\email{leroy.boris@gmail.com})
}

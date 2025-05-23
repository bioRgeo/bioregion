% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netclu_louvain.R
\name{netclu_louvain}
\alias{netclu_louvain}
\title{Louvain community finding}
\usage{
netclu_louvain(
  net,
  weight = TRUE,
  cut_weight = 0,
  index = names(net)[3],
  lang = "igraph",
  resolution = 1,
  seed = NULL,
  q = 0,
  c = 0.5,
  k = 1,
  bipartite = FALSE,
  site_col = 1,
  species_col = 2,
  return_node_type = "both",
  binpath = "tempdir",
  check_install = TRUE,
  path_temp = "louvain_temp",
  delete_temp = TRUE,
  algorithm_in_output = TRUE
)
}
\arguments{
\item{net}{The output object from \code{\link[=similarity]{similarity()}} or
\code{\link[=dissimilarity_to_similarity]{dissimilarity_to_similarity()}}.
If a \code{data.frame} is used, the first two columns represent pairs of sites
(or any pair of nodes), and the next column(s) are the similarity indices.}

\item{weight}{A \code{boolean} indicating if the weights should be considered
if there are more than two columns.}

\item{cut_weight}{A minimal weight value. If \code{weight} is TRUE, the links
between sites with a weight strictly lower than this value will not be
considered (\code{0} by default).}

\item{index}{The name or number of the column to use as weight. By default,
the third column name of \code{net} is used.}

\item{lang}{A string indicating which version of Louvain should be used
(\code{"igraph"} or \code{"cpp"}, see Details).}

\item{resolution}{A resolution parameter to adjust the modularity
(1 is chosen by default, see Details).}

\item{seed}{The random number generator seed (only when \code{lang = "igraph"},
NULL for random by default).}

\item{q}{The quality function used to compute the partition of the graph
(modularity is chosen by default, see Details).}

\item{c}{The parameter for the Owsinski-Zadrozny quality function
(between 0 and 1, 0.5 is chosen by default).}

\item{k}{The kappa_min value for the Shi-Malik quality function
(it must be > 0, 1 is chosen by default).}

\item{bipartite}{A \code{boolean} indicating if the network is bipartite
(see Details).}

\item{site_col}{The name or number for the column of site nodes
(i.e., primary nodes).}

\item{species_col}{The name or number for the column of species nodes
(i.e., feature nodes).}

\item{return_node_type}{A \code{character} indicating what types of nodes
(\code{"site"}, \code{"species"}, or \code{"both"}) should be returned in the output
(\code{"both"} by default).}

\item{binpath}{A \code{character} indicating the path to the bin folder
(see \link{install_binaries} and Details).}

\item{check_install}{A \code{boolean} indicating if the function should check that
Louvain has been properly installed (see \link{install_binaries} and Details).}

\item{path_temp}{A \code{character} indicating the path to the temporary folder
(see Details).}

\item{delete_temp}{A \code{boolean} indicating if the temporary folder should
be removed (see Details).}

\item{algorithm_in_output}{A \code{boolean} indicating if the original output
of \link[igraph:cluster_louvain]{cluster_louvain} should be returned in the
output (\code{TRUE} by default, see Value).}
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

In the \code{algorithm} slot, if \code{algorithm_in_output = TRUE}, users can
find the output of \link[igraph:cluster_louvain]{cluster_louvain} if
\code{lang = "igraph"} and the following element if \code{lang = "cpp"}:

\itemize{
\item{\code{cmd}: The command line used to run Louvain.}
\item{\code{version}: The Louvain version.}
\item{\code{web}: The Louvain's website.}
}
}
\description{
This function finds communities in a (un)weighted undirected network based
on the Louvain algorithm.
}
\details{
Louvain is a network community detection algorithm proposed in
(Blondel et al., 2008). This function offers two
implementations of the Louvain algorithm (controlled by the \code{lang} parameter):
the \href{https://cran.r-project.org/package=igraph}{igraph}
implementation (\link[igraph:cluster_louvain]{cluster_louvain}) and the C++
implementation (\url{https://sourceforge.net/projects/louvain/}, version 0.3).

The \href{https://cran.r-project.org/package=igraph}{igraph}
implementation allows adjustment of the resolution parameter of
the modularity function (\code{resolution} argument) used internally by the
algorithm. Lower values typically yield fewer, larger clusters. The original
definition of modularity is recovered when the resolution parameter
is set to 1 (by default).

The C++ implementation provides several quality functions:
\code{q = 0} for the classical Newman-Girvan criterion (Modularity),
\code{q = 1} for the Zahn-Condorcet criterion, \code{q = 2} for the Owsinski-Zadrozny
criterion (parameterized by \code{c}), \code{q = 3} for the Goldberg Density criterion,
\code{q = 4} for the A-weighted Condorcet criterion, \code{q = 5} for the Deviation to
Indetermination criterion, \code{q = 6} for the Deviation to Uniformity criterion,
\code{q = 7} for the Profile Difference criterion, \code{q = 8} for the Shi-Malik
criterion (parameterized by \code{k}), and \code{q = 9} for the Balanced Modularity
criterion.

The C++ version is based on version 0.3
(\url{https://sourceforge.net/projects/louvain/}). Binary files are required to run it,
and can be installed with \link{install_binaries}.

\strong{If you changed the default path to the \code{bin} folder
while running \link{install_binaries}, PLEASE MAKE SURE to set \code{binpath}
accordingly.}

\strong{If you did not use \link{install_binaries} to change the permissions or test
the binary files, PLEASE MAKE SURE to set \code{check_install} accordingly.}

The C++ version generates temporary folders and/or files in the \code{path_temp}
folder ("louvain_temp" with a unique timestamp located in the bin folder in
\code{binpath} by default). This temporary folder is removed by default
(\code{delete_temp = TRUE}).
}
\note{
Although this algorithm was not primarily designed to deal with bipartite
networks, it is possible to consider the bipartite network as a unipartite
network (\code{bipartite = TRUE}).

Do not forget to indicate which of the first two columns is dedicated to the
site nodes (i.e., primary nodes) and species nodes (i.e., feature nodes) using
the arguments \code{site_col} and \code{species_col}. The type of nodes returned in
the output can be chosen with the argument \code{return_node_type} equal to
\code{"both"} to keep both types of nodes, \code{"site"} to preserve only the site
nodes, and \code{"species"} to preserve only the species nodes.
}
\examples{
comat <- matrix(sample(1000, 50), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

net <- similarity(comat, metric = "Simpson")
com <- netclu_louvain(net, lang = "igraph")

}
\references{
Blondel VD, Guillaume JL, Lambiotte R & Mech ELJS (2008) Fast unfolding of
communities in large networks. \emph{J. Stat. Mech.} 10, P10008.
}
\seealso{
For more details illustrated with a practical example,
see the vignette:
\url{https://biorgeo.github.io/bioregion/articles/a4_3_network_clustering.html}.

Associated functions:
\link{netclu_infomap} \link{netclu_greedy} \link{netclu_oslom}
}
\author{
Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
Pierre Denelle (\email{pierre.denelle@gmail.com}) \cr
Boris Leroy (\email{leroy.boris@gmail.com})
}

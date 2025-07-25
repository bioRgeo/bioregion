% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bind_pairwise.R
\name{bind_pairwise}
\alias{bind_pairwise}
\title{Combine and enrich bioregion (dis)similarity object(s)}
\usage{
bind_pairwise(primary_metrics, secondary_metrics, new_metrics = NULL)
}
\arguments{
\item{primary_metrics}{A \code{bioregion.pairwise} object. This is the
main set of pairwise metrics that will be used as a base for the combination.}

\item{secondary_metrics}{A second \code{bioregion.pairwise}
object to be combined with \code{primary_metrics}. It must have the same sites
identifiers and pairwise structure. Can be set to \code{NULL} if \code{new_metrics} is
specified.}

\item{new_metrics}{A \code{character} vector or a single \code{character} string
specifying custom formula(s) based on the column names of \code{primary_metrics}
and \code{secondary_metrics} (see Details).}
}
\value{
A new \code{bioregion.pairwise} object containing the combined and/or
enriched data. It includes all original metrics from the inputs, as well as
any newly computed metrics.
}
\description{
Combine two \code{bioregion.pairwise} objects and/or compute new pairwise
metrics based on the columns of the object(s).
}
\details{
When both \code{primary_metrics} and \code{secondary_metrics} are provided and if the
pairwise structure is identical the function combine the two objects. If
\code{new_metrics} is provided, each formula is evaluated based on the column
names of \code{primary_metrics} (and \code{secondary_metrics} if provided).
}
\examples{
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("s", 1:5)
colnames(comat) <- paste0("sp", 1:10)

sim <- bind_pairwise(primary_metrics = similarity(comat, 
                                                               metric = "abc"),
                                  secondary_metrics = similarity(comat, 
                                                                 metric = "Simpson"),
                                  new_metrics = "1 - (b + c) / (a + b + c)")

}
\seealso{
For more details illustrated with a practical example,
see the vignette:
\url{https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html}.

Associated functions:
\link{dissimilarity} \link{similarity} \link{as_bioregion_pairwise}
}
\author{
Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
Boris Leroy (\email{leroy.boris@gmail.com}) \cr
Pierre Denelle (\email{pierre.denelle@gmail.com})
}

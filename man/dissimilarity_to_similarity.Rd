% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similarity_dissimilarity_conversion.R
\name{dissimilarity_to_similarity}
\alias{dissimilarity_to_similarity}
\title{Convert dissimilarity metrics to similarity metrics}
\usage{
dissimilarity_to_similarity(dissimilarity, include_formula = TRUE)
}
\arguments{
\item{dissimilarity}{the output object from \code{\link[=dissimilarity]{dissimilarity()}} or
\code{\link[=similarity_to_dissimilarity]{similarity_to_dissimilarity()}}.}

\item{include_formula}{a \code{boolean} indicating whether metrics based on custom
formula(s) should also be converted (see Details). The default is \code{TRUE}.}
}
\value{
A \code{data.frame} with the additional class
\code{bioregion.pairwise}, providing similarity metrics for each pair of
sites based on a dissimilarity object.
}
\description{
This function converts a \code{data.frame} of dissimilarity metrics
(beta diversity) between sites into similarity metrics.
}
\note{
The behavior of this function changes depending on column names. Columns
\code{Site1} and \code{Site2} are copied identically. If there are columns called
\code{a}, \code{b}, \code{c}, \code{A}, \code{B}, \code{C} they will also be copied identically. If there
are columns based on your own formula (argument \code{formula} in
\code{\link[=dissimilarity]{dissimilarity()}}) or not in the original list of dissimilarity metrics
(argument \code{metrics} in \code{\link[=dissimilarity]{dissimilarity()}}) and if the argument
\code{include_formula} is set to \code{FALSE}, they will also be copied identically.
Otherwise there are going to be converted like they other columns (default
behavior).

If a column is called \code{Euclidean}, the similarity will be calculated based
on the following formula:

Euclidean similarity = 1 / (1 - Euclidean distance)

Otherwise, all other columns will be transformed into dissimilarity with the
following formula:

similarity = 1 - dissimilarity
}
\examples{
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("s", 1:5)
colnames(comat) <- paste0("sp", 1:10)

dissimil <- dissimilarity(comat, metric = "all")
dissimil

similarity <- dissimilarity_to_similarity(dissimil)
similarity

}
\seealso{
For more details illustrated with a practical example,
see the vignette:
\url{https://biorgeo.github.io/bioregion/articles/a3_pairwise_metrics.html}.

Associated functions:
\link{similarity} \link{dissimilarity_to_similarity}
}
\author{
Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
Boris Leroy (\email{leroy.boris@gmail.com}) \cr
Pierre Denelle (\email{pierre.denelle@gmail.com})
}

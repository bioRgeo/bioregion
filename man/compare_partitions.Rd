% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_partitions.R
\name{compare_partitions}
\alias{compare_partitions}
\title{Compare cluster memberships among multiple partitions}
\usage{
compare_partitions(
  cluster_object,
  sample_items = NULL,
  indices = c("rand", "jaccard"),
  cor_frequency = FALSE,
  store_pairwise_membership = TRUE,
  store_confusion_matrix = TRUE
)
}
\arguments{
\item{cluster_object}{a \code{bioregion.clusters} object or a \code{data.frame} or a
list of \code{data.frame} containing multiple partitions. At least two partitions
are required. If a list of \code{data.frame} is provided, they should all have
the same number of rows (i.e., same items in the clustering for all
partitions).}

\item{sample_items}{\code{NULL} or a positive integer. Reduce the number of items
to be used in the comparison of partitions. Useful if the number of items is
high and pairwise comparisons cannot be computed. Suggested values
5000 or 10000 computation}

\item{indices}{\code{NULL} or \code{character}. Indices to compute for the pairwise
comparison of partitions. Current available metrics are \code{"rand"} and
\code{"jaccard"}}

\item{cor_frequency}{a boolean. If \code{TRUE}, then computes the correlation
between each partition and the total frequency of co-membership of items
across all partitions. Useful to identify which partition(s) is(are) most
representative of all the computed partitions.}

\item{store_pairwise_membership}{a boolean. If \code{TRUE}, the pairwise
membership of items is stored in the output object.}

\item{store_confusion_matrix}{a boolean. If \code{TRUE}, the confusion matrices
of pairwise partition comparisons are stored in the output object.}
}
\value{
A \code{list} with 4 to 7 elements:
\itemize{
\item \code{args}: arguments provided by the user
\item \code{inputs}: information on the input partitions, such as the number of items
being clustered
\item (facultative) \code{pairwise_membership}: only if
\code{store_pairwise_membership = TRUE}. This
element contains the pairwise memberships of all items for each
partition, in the form of a \verb{boolean matrix} where \code{TRUE} means that
two items are in the same cluster, and \code{FALSE} means that two items
are not in the same cluster
\item \code{freq_item_pw_membership}: A \verb{numeric vector}
containing the number of times each pair of items are clustered
together. It corresponds to the sum of rows of the table in
\code{pairwise_membership}
\item (facultative) \code{partition_freq_cor}:  only if \code{cor_frequency = TRUE}.
A \verb{numeric vector}
indicating the correlation between individual partitions and the total
frequency of pairwise membership across all partitions. It corresponds to
the correlation between individual columns in \code{pairwise_membership} and
\code{freq_item_pw_membership}
\item (facultative) \code{confusion_matrix}: only if \code{store_confusion_matrix = TRUE}.
A \code{list}
containing all confusion matrices between each pair of partitions.
\item \code{partition_comparison}: a \code{data.frame} containing the results of the
comparison of partitions, where the first column indicates which partitions
are compared, and the next columns correspond to the requested \code{indices}.
}
}
\description{
This function aims at computing pairwise comparisons for several
partitions, usually on outputs from \code{netclu_}, \code{hclu_} or \code{nhclu_} functions.
It also provides the confusion matrix from pairwise comparisons, so that
the user can compute additional comparison metrics.
}
\details{
\loadmathjax
This function proceeds in two main steps:
\enumerate{
\item The first step is done within each partition. It will compare all pairs of
items and document if they are clustered together (\code{TRUE}) or separately
(\code{FALSE}) in each partition. For example, if site 1 and site 2 are clustered
in the same cluster in partition 1, then the pairwise membership site1_site2
will be \code{TRUE}. The output of this first step is stored in the slot
\code{pairwise_membership} if \code{store_pairwise_membership = TRUE}.
\item The second step compares all pairs of partitions by analysing if their
pairwise memberships are similar or not. To do so, for each pair of
partitions, the function computes a confusion matrix with four elements:
}
\itemize{
\item \emph{a}: number of pairs of items grouped in partition 1 and in partition 2
\item \emph{b}: number of pairs of items grouped in partition 1 but not in partition
2
\item \emph{c}: number of pairs of items not grouped in partition 1 but grouped in
partition 2
\item \emph{d}: number of pairs of items not grouped in both partition 1 & 2
}

The confusion matrix is stored in \code{confusion_matrix} if
\code{store_confusion_matrix = TRUE}.

Based on the confusion matrices, we can compute a range of indices to
indicate the agreement among partitions. As of now, we have implemented:
\itemize{
\item \emph{Rand index}
\mjeqn{(a + d)/(a + b + c + d)}{(a + d)/(a + b + c + d)}
The Rand index measures agreement among partitions by accounting for both
the pairs of sites that are grouped, but also the pairs of sites that are
not grouped.
\item \emph{Jaccard index}
\mjeqn{(a)/(a + b + c)}{(a)/(a + b + c)}
The Jaccard index measures agreement among partitions by only accounting
for pairs of sites that are grouped - it is
}

These two metrics are complementary, because the Jaccard index will tell
if partitions are similar in their clustering structure, whereas the
Rand index will tell if partitions are similar not only in the pairs of
items clustered together, but also in terms of the pairs of sites that are
not clustered together. For example, take two partitions which
never group together the same pairs of sites. Their Jaccard index will be 0,
whereas the Rand index can be > 0 due to the sites that are not grouped
together.

Additional indices can be manually computed by the users on the basis of the
list of confusion matrices.

In some cases, users may be interested in finding which of the partitions
is most representative of all partitions. To find it out, we can
compare the pairwise membership of each partition with the total frequency
of pairwise membership across all partitions. This correlation can be
requested with \code{cor_frequency = TRUE}
}
\examples{
# A simple case with four partitions of four items
partitions <- data.frame(matrix(nr = 4, nc = 4, 
                                c(1,2,1,1,1,2,2,1,2,1,3,1,2,1,4,2),
                                byrow = TRUE))
partitions
compare_partitions(partitions)

# Find out which partitions are most representative
compare_partitions(partitions,
                   cor_frequency = TRUE)
                                

}
\seealso{
\link{partition_metrics}
}
\author{
Boris Leroy (\email{leroy.boris@gmail.com}),
Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) and
Pierre Denelle (\email{pierre.denelle@gmail.com})
}

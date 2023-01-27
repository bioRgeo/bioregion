#' Spatial distribution of Mediterranean vegetation (data.frame)
#'
#' A dataset containing the abundance of 3,697 species in 715 sites.
#'
#' @format A `data.frame` with 460,878 rows and 3 columns:
#' \describe{
#'   \item{Site}{Unique site identifier (corresponding to the field ID of vegesp).}
#'   \item{Species}{Unique species identifier.}
#'   \item{Abundance}{Species abundance}
#' }
#' @source <https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4718>
"vegedf"

#' Spatial distribution of Mediterranean vegetation (co-occurrence matrix)
#'
#' A dataset containing the abundance of each of the 3,697 species in each of
#' the 715 sites.
#'
#' @format A co-occurrence `matrix` with sites as rows and species as
#' columns. Each element of the matrix
#' represents the abundance of the species in the site.
#' @source <https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4718>
"vegemat"

#' Spatial distribution of Mediterranean vegetation (spatial grid)
#'
#' A dataset containing the geometry of the 715 sites.
#'
#' @format A
#' \describe{
#'   \item{ID}{Unique site identifier.}
#'   \item{geometry}{Geometry of the site.}
#' }
#' @source <https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4718>
"vegesp"

#' Spatial distribution of benthic invertebrates in Brittany (data.frame)
#'
#' A dataset containing the abundance of 3,697 species in 715 sites.
#'
#' @format A `data.frame` with 28,563 rows and 3 columns:
#' \describe{
#'   \item{Site}{Unique site identifier (corresponding to the field ID of vegesp).}
#'   \item{Species}{Unique species identifier.}
#'   \item{Abundance}{Species abundance}
#' }
"benthicdf"

#' Spatial distribution of benthic invertebrates in Brittany
#' (co-occurrence matrix)
#'
#' A dataset containing the abundance of each of the 838 species in each of
#' the 162 sites.
#'
#' @format A co-occurrence `matrix` with sites as rows and species as
#' columns. Each element of the matrix
#' represents the abundance of the species in the site.
"benthicmat"

#' Spatial distribution of benthic invertebrates in Brittany (spatial points)
#'
#' A dataset containing the geometry of the 162 sites.
#'
#' @format A
#' \describe{
#'   \item{ID}{Unique site identifier.}
#'   \item{geometry}{Geometry of the site.}
#' }
"benthicsp"

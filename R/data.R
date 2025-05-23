#' Spatial distribution of Mediterranean vegetation (data.frame)
#'
#' A dataset containing the abundance of 3,697 species in 715 sites.
#'
#' @format A `data.frame` with 460,878 rows and 3 columns:
#' \describe{
#'   \item{Site}{Unique site identifier (corresponding to the field ID of vegesp)}
#'   \item{Species}{Unique species identifier}
#'   \item{Abundance}{Species abundance}
#' }
#' @source \doi{10.1002/ece3.4718} 
"vegedf"

#' Spatial distribution of Mediterranean vegetation (co-occurrence matrix)
#'
#' A dataset containing the abundance of each of the 3,697 species in each of
#' the 715 sites.
#'
#' @format A co-occurrence `matrix` with sites as rows and species as
#' columns. Each element of the matrix
#' represents the abundance of the species in the site.
#' @source \doi{10.1002/ece3.4718} 
"vegemat"

#' Spatial distribution of Mediterranean vegetation (spatial grid)
#'
#' A dataset containing the geometry of the 715 sites.
#'
#' @format A
#' \describe{
#'   \item{ID}{Unique site identifier}
#'   \item{geometry}{Geometry of the site}
#' }
#' @source \doi{10.1002/ece3.4718} 
"vegesf"

#' Spatial distribution of fish in Europe (data.frame)
#'
#' A dataset containing the abundance of 195 species in 338 sites.
#'
#' @format A `data.frame` with 2,703 rows and 3 columns:
#' \describe{
#'   \item{Site}{Unique site identifier (corresponding to the field ID of fishsf)}
#'   \item{Species}{Unique species identifier}
#'   \item{Abundance}{Species abundance}
#' }
"fishdf"

#' Spatial distribution of fish in Europe (co-occurrence matrix)
#'
#' A dataset containing the abundance of each of the 195 species in each of
#' the 338 sites.
#'
#' @format A co-occurrence `matrix` with sites as rows and species as
#' columns. Each element of the matrix
#' represents the abundance of the species in the site.
"fishmat"

#' Spatial distribution of fish in Europe
#'
#' A dataset containing the geometry of the 338 sites.
#'
#' @format A
#' \describe{
#'   \item{ID}{Unique site identifier}
#'   \item{geometry}{Geometry of the site}
#' }
"fishsf"

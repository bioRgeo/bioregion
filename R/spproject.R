#' Compute similarity metrics between sites based on species composition
#'
#' This function creates a \code{data.frame} where each row provides one or several similarity
#' metric(s) between each pair of sites from a co-occurence \code{matrix} with sites as rows and species as columns.
#'
#' @param comat a co-occurence \code{matrix} with sites as rows and species as columns
#' @param metric a vector of character(s) indicating which similarity
#' metric(s) to chose (see Details)
#' @export
#' @details
#' With \code{a} the number of species shared by a pair of sites, \code{b} species only present in the first site
#' and \code{c} species only present in the second site.
#'
#'\eqn{Jaccard = 1 - (b+c)/(a + b + c)}
#'
#'\eqn{Jaccardturn = 1 - 2min(b, c)/(a + 2min(b, c))} \insertCite{Baselga2012}{bioRgeo}
#'
#'\eqn{Sorensen = 1 - (b+c)/(2a + b + c)}
#'
#'\eqn{Simpson = 1 - min(b, c)/(a + min(b, c))}
#'
#'If abundances data are available, Bray-Curtis and its turnover component can also be computed with the
#'following equation:
#'
#'\eqn{Bray = 1- (B+C)/(2A + B + C)}
#'
#'\eqn{Brayturn = 1 - min(B, C)/(A + min(B, C))} \insertCite{Baselga2013}{bioRgeo}
#'
#'with A the sum of the lesser values for common species shared by a pair of sites.
#'B and C are the total number of specimens counted at both sites minus A.
#'
#' Euclidean computes the Euclidean similarity between each pair of site following this equation:
#'
#' \eqn{Euclidean = 1 /(1 + dij)}
#'
#' Where dij is the Euclidean distance between site i and site j in terms of species composition.
#'
#' @return A \code{data.frame} providing one or several similarity
#' metric(s) between each pair of sites. The two first columns represents each pair of sites.
#' One column per similarity metric provided in \code{metric} except for the metric \emph{abc} and \emph{ABC} that
#' are stored in three columns (one for each letter).
#' @author
#' Pierre Denelle \email{pierre.denelle@gmail.com}
#' Maxime Lenormand \email{maxime.lenormand@inrae.fr}
#' Boris Leroy \email{leroy.boris@gmail.com}
#' @examples
#' comat=matrix(sample(1000,50),5,10)
#' rownames(comat)=paste0("Site",1:5)
#' colnames(comat)=paste0("Species",1:10)
#'
#' simil=spproject(comat,metric=c("abc","ABC","Simpson","Brayturn"))
#' @references
#' \insertRef{Baselga2012}{bioRgeo}
#' \insertRef{Baselga2013}{bioRgeo}
#' @export
spproject <- function(comat, metric = "Simpson"){

  # list of metrics based on abc
  lsmetricabc=c("abc","Jaccard","Jaccardturn","Sorensen","Simpson")

  # list of metrics based on ABC
  lsmetricABC=c("ABC","Bray","Brayturn")

  # list of metrics based on other features
  lsmetrico=c("Euclidean")


  # Controls
  require(Rcpp)

  if(length(intersect(c(lsmetricabc,lsmetricABC,lsmetrico),metric))!=length(metric)){
    stop("One or several similarity metric(s) chosen is not available.
     Please chose among the followings:
         abc, Jaccard, Jaccardturn, Sorensen, Simpson, ABC, Bray, Brayturn or Euclidean")
  }

  if(!is.matrix(comat)){
    stop("Co-occurence matrix should be a matrix")
  }

  sco=sum(is.na(comat))
  minco=min(comat)
  if(sco>0){
    stop("Co-occurence matrix should contains only positive real: NA(s) detected!")
  }
  if(minco<0){
    stop("Co-occurence matrix should contains only positive real: negative value detected!")
  }

  # Initialize output res (two-column data.frame containing each pair of sites to be compared)
  siteid=rownames(comat)
  res=matrix(0, nrow=length(siteid), ncol=length(siteid), dimnames=list(siteid,siteid))
  res[upper.tri(res)]=1
  res=contingency_to_df(res, weight=FALSE, remove_absent_objects = TRUE)
  colnames(res)=c("Site1","Site2")

  # abcp: compute abc for presence data
  if(length(intersect(lsmetricabc,metric))>0){

      # Compute the number of species in common "a" with matricial product a%*%t(a)
      comatp=comat
      comatp[comatp!=0]=1
      sumrow=apply(comatp,1,sum)
      abcp=prodmat(comatp,t(comatp))
      rownames(abcp)=siteid
      colnames(abcp)=siteid

      # Create a data.frame from the matrix with contingency_to_df (little trick to deal with 0s)
      abcp[abcp==0]=-1
      abcp[lower.tri(abcp, diag=TRUE)]=0
      abcp=contingency_to_df(abcp, weight=TRUE, remove_absent_objects = TRUE)
      colnames(abcp)=c("Site1","Site2","a")
      abcp[abcp[,3]==-1,3]=0

      # Compute b and c
      abcp$b=0
      abcp$c=0
      abcp[,4]=sumrow[match(abcp[,1],siteid)]-abcp[,3]
      abcp[,5]=sumrow[match(abcp[,2],siteid)]-abcp[,3]

      # Compute metrics based on abc
      if("abc" %in% metric){
        res$a=abcp$a
        res$b=abcp$b
        res$c=abcp$c
      }
      if("Jaccard" %in% metric){
        res$Jaccard = 1 - (abcp$b+abcp$c)/(abcp$a+abcp$b+abcp$c)
      }
      if("Jaccardturn" %in% metric){
        res$Jaccardturn = 1 - 2*pmin(abcp$b,abcp$c)/(abcp$a + 2*pmin(abcp$b,abcp$c))
      }
      if("Sorensen" %in% metric){
        res$Sorensen = 1 - (abcp$b+abcp$c)/(2*abcp$a+abcp$b+abcp$c)
      }
      if("Simpson" %in% metric){
        res$Simpson = 1 - pmin(abcp$b,abcp$c)/(abcp$a + pmin(abcp$b,abcp$c))
      }
  }

  # abca: compute ABC for abundance data
  if(length(intersect(lsmetricABC,metric))>0){

    # Use abc Rcpp function in src (three loops)
    abca=abc(comat)
    abca=data.frame(Site1=siteid[abca[,1]], Site2=siteid[abca[,2]], A=abca[,3], B=abca[,4]-abca[,3], C=abca[,5]-abca[,3])

    # Compute metrics based on ABC
    if("ABC" %in% metric){
      res$A=abca$A
      res$B=abca$B
      res$C=abca$C
    }
    if("Bray" %in% metric){
      res$Bray = 1 - (abca$B+abca$C)/(2*abca$A+abca$B+abca$C)
    }
    if("Brayturn" %in% metric){
      res$Brayturn = 1 - pmin(abca$B,abca$C)/(abca$A + pmin(abca$B,abca$C))
    }
  }

  # Compute Euclidean similarity between site using dist()
  if("Euclidean" %in% metric){
    eucl=as.matrix(dist(comat))
    rownames(eucl)=siteid
    colnames(eucl)=siteid
    eucl[eucl==0]=-1
    eucl[lower.tri(eucl, diag=TRUE)]=0
    eucl=contingency_to_df(eucl, weight=TRUE, remove_absent_objects = TRUE)
    colnames(eucl)=c("Site1","Site2","Euclidean")
    eucl[eucl[,3]==-1,3]=0

    res$Euclidean=1/(1+eucl$Euclidean)
  }

  # Return the output
  return(res)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install_binaries.R
\name{install_binaries}
\alias{install_binaries}
\title{Download, unzip, check permission and test the bioregion's binary files}
\usage{
install_binaries(
  binpath = "tempdir",
  download_only = FALSE,
  infomap_version = c("2.1.0", "2.6.0", "2.7.1", "2.8.0")
)
}
\arguments{
\item{binpath}{a \code{character} indicating the path to the folder that will
host the 'bin' folder containing the binary files (see Details).}

\item{download_only}{a \code{boolean} indicating whether the function should only
download the 'bin.zip' file or execute the entire process (see Details).}

\item{infomap_version}{a \code{character} vector or a single \code{character} string
indicating the Infomap version(s) to install.}
}
\value{
No return value
}
\description{
This function downloads and unzips the 'bin' folder needed to run some
functions of bioregion. It also checks if the files have the permissions to
be executed as programs. It finally tests if the binary files
are running properly.
}
\details{
By default, the binary files are installed in R's temporary
directory (\code{binpath = "tempdir"}). In this case the \code{bin} folder will be
automatically removed at the end of the R session. Alternatively, the binary
files can be installed in the bioregion's package folder
(\code{binpath = "pkgfolder"}).

Finally, a path to a folder of your choice can be specified. In this case
, and only in this case, \code{download_only} can be set to \code{TRUE}, but you must
ensure that the files have the necessary permissions to be executed as
programs.

\strong{In any case, PLEASE MAKE SURE to update the \code{binpath} and \code{check_install}
accordingly in \link{netclu_infomap}, \link{netclu_louvain} and \link{netclu_oslom}).}
}
\note{
Only the Infomap version 2.1.0, 2.6.0, 2.7.1 and 2.8.0 are available for now.
}
\author{
Maxime Lenormand (\email{maxime.lenormand@inrae.fr}) \cr
Boris Leroy (\email{leroy.boris@gmail.com}) \cr
Pierre Denelle (\email{pierre.denelle@gmail.com})
}

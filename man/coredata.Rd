\name{coredata}
\alias{coredata}
\alias{coredata.colorSpec}
\alias{as.matrix.colorSpec}

\title{Extract the Core Data of a colorSpec Object}
\description{
functions for extracting the core data contained in a \bold{colorSpec} object.
}

\usage{
\S3method{coredata}{colorSpec}( x, forcemat=FALSE )

\S3method{as.matrix}{colorSpec}( x, ... )
}

\arguments{
  \item{x}{a \bold{colorSpec} object}
  \item{forcemat}{if \code{x} has only 1 spectrum, return a matrix with 1 column instead of a vector}
  \item{...}{extra arguments ignored}  
}

\value{
\item{coredata}{
If \code{x} has \code{organization} equal to \code{'vector'} then it returns \code{x},
unless \code{forcemat} is \code{TRUE} when it returns \code{x} as a matrix with 1 column.\cr
If \code{x} has any other \code{organization} then it returns a matrix
with spectra in the columns.
All of the \bold{colorSpec} attributes are stripped except the column names,
and the row names are set to \code{as.character(wavelength(x))}.}

\item{as.matrix}{a wrapper for \code{coredata} with \code{forcemat=TRUE}}
}


\seealso{
\code{\link{organization}}
}


\keyword{colorSpec}

% File specnames.Rd

\name{specnames}
\title{specnames of a colorSpec object}
\alias{specnames}
\alias{specnames<-}
\alias{numSpectra}
\alias{specnames.colorSpec}
\alias{specnames<-.colorSpec}
\alias{numSpectra.colorSpec}

\description{
  Retrieve or set the specnames of a \bold{colorSpec} object.  Retrieve the number of spectra.
}
\usage{
\S3method{specnames}{colorSpec}(x)

\S3method{specnames}{colorSpec}(x) <- value

\S3method{numSpectra}{colorSpec}(x)

}

\arguments{
  \item{x}{a \bold{colorSpec} \R object}
  
  \item{value}{a character vector with length equal to the number of spectra in \code{x}.}
}

\details{
  If the \code{organization} of \code{x} is \code{"vector"} then \code{x} is a vector
  and \code{value} is a single string, which is stored as \code{attr(x,'specname')}.

  If the \code{organization} of \code{x} is \code{"matrix"}, then \code{x} is a matrix
  and \code{value} is stored as \code{colnames(x)}.

  If the \code{organization} of \code{x} is \code{"df.col"}, then \code{x} is a \code{data.frame}
  with \code{N+1} columns, where \code{N} is the number of spectra.
  \code{value} is stored as \code{colnames(x)[2:(N+1)]}.

  If the \code{organization} of \code{x} is \code{"df.row"}, then \code{x} is a \code{data.frame}
  and \code{value} is stored as \code{row.names(x)}.
}

\value{
\code{specnames()} returns a character vector with the names of the spectra.

\code{numSpectra(x)} is equal to \code{length(specnames(x))} but much more efficient.
}

\seealso{
  \code{\link{rownames}}, \code{\link{colnames}}
}

\keyword{colorSpec}

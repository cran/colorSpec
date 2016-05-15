% File organization.Rd


\name{organization}
\title{organization of a colorSpec object}
\alias{organization}
\alias{organization<-}
\alias{organization.colorSpec}
\alias{organization<-.colorSpec}

\description{
  Retrieve or set the organization of a \bold{colorSpec} object.
}
\usage{
\S3method{organization}{colorSpec}(x)
organization(x) <- value

}

\arguments{
  \item{x}{a \bold{colorSpec} \R object}
  
  \item{value}{a valid organization: \code{'vector'}, \code{'matrix'}, \code{'df.col'}, or \code{'df.row'}.}
}

\details{
  If the \code{organization} of \code{x} is \code{"vector"}, then \code{x} is a vector
  representing a single spectrum.  Compare this with \code{\link{ts}}.

  If the \code{organization} of \code{x} is \code{"matrix"}, then \code{x} is a matrix
  and the spectra are stored in the columns.

  If the \code{organization} of \code{x} is \code{"df.col"}, then \code{x} is a \code{data.frame}
  with \code{n+1} columns, where \code{n} is the number of spectra.
  The wavelengths are stored in column 1, and the spectra in columns 2:(n+1).
  This organization is good for printing to the console, and writing to files.

  If the \code{organization} of \code{x} is \code{"df.row"}, then \code{x} is a \code{data.frame}
  with \code{n} rows, where \code{n} is the number of spectra.
  The spectra are stored in the last column, which is a \code{"model.matrix"} with the name \code{"spectra"}.
  The other columns preceding \code{spectra} (if present) contain additional data associated with
  the spectra; see \code{\link{extradata}}.
}

\value{
  \code{organization} returns a valid organization: \code{'vector'}, \code{'matrix'}, \code{'df.col'}, or \code{'df.row'}.
}

\note{
If \code{x} has more than 1 spectrum, then \code{organization} of \code{"vector"} is invalid and ignored.

If \code{x} has \code{organization} equal to \code{df.row} and also has \code{\link{extradata}},
then changing the \code{organization} will discard the \code{extradata}.
}

\seealso{
\code{\link{colorSpec}}; \code{\link{extradata}}
}

\examples{
organization(Hoya)              # returns 'df.row'
organization(Hoya) = 'matrix'   # extradata in Hoya is silently discarded
}

\keyword{colorSpec}

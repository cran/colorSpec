% File metadata.Rd

\name{metadata}
\title{metadata of a colorSpec object}
\alias{metadata}
\alias{metadata<-}
\alias{metadata.colorSpec}
\alias{metadata<-.colorSpec}

\description{
  Retrieve or set the metadata of a \bold{colorSpec} object.
}
\usage{
\S3method{metadata}{colorSpec}(x, ...)
metadata(x) <- value

}

\arguments{
  \item{x}{a \bold{colorSpec} \R object}
  
  \item{...}{optional names of metadata to return}
  
  \item{value}{a named \code{list} that is appended to the existing list of metadata.
  Unnamed items are ignored.  If a name already exists, its value is updated.
  }
}

\details{
The metadata list is stored as \code{attr(x,'metadata')}.
Initially this list is empty.
\cr
The assignment operator looks like it replaces the internal metadata list of \code{x},
but it actually appends to it.
}

\value{
\code{metadata()} with no argument returns the complete named list of metadata.
If arguments are present, then only those \code{metadata} items are returned.

}


\seealso{
  \code{\link{colorSpec}}; \code{\link{extradata}}
}

\examples{

\dontrun{
# get list of *all* metadata
metadata(x)

# get just the file 'path'
metadata( x, 'path' )

# set the 'date'
metadata( x ) = list( date="2016-04-01" )
}
}

\keyword{colorSpec}

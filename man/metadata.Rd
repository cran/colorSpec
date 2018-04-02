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

\S3method{metadata}{colorSpec}(x, add=FALSE ) <- value
}

\arguments{
  \item{x}{a \bold{colorSpec} \R object}
  
  \item{...}{optional names of metadata to return}
  
  \item{value}{a named \code{list}.  
  If \code{add} is \code{FALSE}, \code{value} replaces any existing \code{metadata}.   
  If \code{add} is \code{TRUE},  \code{value} is appended to the existing list of metadata.
  If a name already exists, its value is updated using \code{\link{modifyList}()}.
  Unnamed items in  \code{value} are ignored.}
  
  \item{add}{if \code{add=FALSE}, any existing metadata is discarded.
  If \code{add=TRUE} then existing metadata is preserved, using \code{modifyList()}.
}
}

\details{
The metadata list is stored as \code{attr(x,'metadata')}.
After construction this list is empty.
}

\value{
\code{metadata(x)} with no additional arguments returns the complete named list of metadata.
If arguments are present, then only those \code{metadata} items are returned.

}

\note{
Do not confuse \code{extradata} and \code{metadata}.\cr
\code{metadata} is unstructured data that is attached to the entire \bold{colorSpec} object.
\code{extradata} is structured data, with a row of data for each spectrum in the object.
}


\seealso{
\code{\link{extradata}},
\code{\link{modifyList}}
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

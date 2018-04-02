% File extradata.Rd

\name{extradata}
\title{extradata of a colorSpec object}
\alias{extradata}
\alias{extradata<-}
\alias{extradata.colorSpec}
\alias{extradata<-.colorSpec}

\description{
  Retrieve or set the extradata of a \bold{colorSpec} object.
}
\usage{
\S3method{extradata}{colorSpec}(x)

\S3method{extradata}{colorSpec}(x,add=FALSE) <- value
}

\arguments{
  \item{x}{a \bold{colorSpec}  object with M spectra}
  
  \item{value}{a \code{data.frame} with M rows.
  It is OK for \code{value} to have 0 columns,
  and \code{value} can also be \code{NULL}; see \code{add}.
}
  
  \item{add}{ If \code{add} is \code{FALSE},
  then any existing extradata is discarded and replaced by \code{value},
  except when \code{value} is \code{NULL} when \code{x} is left with no \code{extradata}.\cr
  If \code{add} is \code{TRUE}, then \code{value} is appended to the existing extradata,
  except when \code{value} is \code{NULL} when \code{x} is left unchanged.}
}

\details{
  If the \code{organization} of \code{x} is not \code{'df.row'}, 
  then \code{extradata} cannot be stored in \code{x} and the assignment is ignored, with a warning.
  First change the \code{\link{organization}} to \code{'df.row'},
  and \emph{then} assign the \code{extradata}.
  
  If the \code{organization} of \code{x} is \code{'df.row'}, but \code{value} does not have
  the right number of rows, the assignment is ignored, with a warning.\cr
  
}

\value{
  \code{extradata(x)} returns a \code{data.frame} with M rows,
  where M is the number of spectra in \code{x}.  
  The \code{rownames} are set to the \code{specnames} of \code{x}.
  If there is no extra data then the number of columns in this \code{data.frame} is 0.
}

\note{
Do not confuse \code{extradata} and \code{metadata}.\cr
\code{metadata} is unstructured data that is attached to the entire \bold{colorSpec} object.
\code{extradata} is structured data, with a row of data for each spectrum in the object.
}

\seealso{
  \code{\link{metadata}}, 
  \code{\link{organization}}
}

\keyword{colorSpec}

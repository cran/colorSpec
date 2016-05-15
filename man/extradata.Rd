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
extradata(x) <- value
}

\arguments{
  \item{x}{a \bold{colorSpec} \R object}
  
  \item{value}{a \code{data.frame} with \code{m} rows, where \code{m} is the number of spectra in \code{x}.}
}

\details{
  If the \code{organization} of \code{x} is not \code{'df.row'}, 
  then \code{extradata} cannot be stored in \code{x} and the assignment is ignored, with a warning.
  First change the \code{\link{organization}} to \code{'df.row'}, and then assign the \code{extradata}.
  
  If the \code{organization} of \code{x} is \code{"df.row"}, but \code{value} does not have
  the right number of rows, the assignment is also ignored.
}

\value{
  \code{extradata} returns a \code{data.frame} with \code{m} rows,
  where \code{m} is the number of spectra in \code{x}.  
  The \code{rownames} are set to the \code{specnames} of \code{x}.
  If there is no extra data then the number of columns in this \code{data.frame} is 0.
}

\note{
If the assigment of \code{value} is successful, any existing \code{extradata} in \code{x} is discarded.
}

\seealso{
  \code{\link{colorSpec}}, 
  \code{\link{organization}}
}

\keyword{colorSpec}

% File bind.Rd

\name{bind}
\title{Combine colorSpec Objects}
\alias{bind}
\alias{bind.colorSpec}

\description{
  Take a sequence of \bold{colorSpec} objects and combine their spectra
}
\usage{
\S3method{bind}{colorSpec}( \ldots )

}

\arguments{
  \item{...}{\bold{colorSpec} objects with the same \code{\link{wavelength}} and \code{\link{quantity}}, and with distinct \code{\link{specnames}} }
    
}
    
\value{
  \code{bind} returns a \bold{colorSpec} object, or \code{NULL} in case of ERROR.  If the bind is
  successful, the number of spectra in the output object is the sum of the number of spectra
  in the input objects.
           
}    


\details{
  The \code{\link{organization}} of the returned object is the most complex of those in the input list,
  where the order of complexity is:
  
  \code{'matrix'}  <  \code{'df.col'}  <  \code{'df.row'}
  
  If the selected \code{organization} is \code{'df.row'} the \code{extradata} is combined using \code{\link{merge}}.
  
  The \code{\link{metadata}} of the returned object is copied from the first object in the input list.
  
}

\seealso{
  \code{\link{wavelength}},
  \code{\link{quantity}},
  \code{\link{specnames}},    
  \code{\link{organization}},  
  \code{\link{extradata}},    
  \code{\link{metadata}},
  \code{\link[base]{merge}}  
}


\keyword{colorSpec}

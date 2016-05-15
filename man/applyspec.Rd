% File applyspec.Rd

\name{applyspec}
\title{apply a function to each spectrum in a colorSpec object}
\alias{applyspec.colorSpec}
\alias{applyspec}

\description{
apply a function to each spectrum in a colorSpec object
}

\usage{
\S3method{applyspec}{colorSpec}( x, FUN, \ldots )
}

\arguments{
  \item{x}{a \bold{colorSpec} object with N wavelengths}
  
  \item{FUN}{a function that takes an N-vector as argument and returns an N-vector}  
  
  \item{\ldots}{additional arguments passed to \code{FUN}}

}
    
\value{
a \bold{colorSpec} object with the same dimensions,
\code{\link{wavelength}}, \code{\link{quantity}}, and \code{\link{organization}}.
If \code{FUN} does not return an N-vector, it is an ERROR and
\code{applyspec} returns \code{NULL}.
}


\details{
\code{applyspec} simply calls \code{\link[base]{apply}} with the right \code{MARGIN}.
}


\seealso{
  \code{\link{quantity}},
  \code{\link{wavelength}},
  \code{\link{linearize}},  
  \code{\link{organization}}
}

\examples{
#  convert absorbance to transmittance
path = system.file( "extdata/stains/Hematoxylin.txt", package='colorSpec' )
x = readSpectra( path )
x = applyspec( x, function(y) {10^(-y)} ) # this is what linearize(x) does
}

\keyword{colorSpec}

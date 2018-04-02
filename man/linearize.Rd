% File linearize.Rd

\name{linearize}
\title{linearize a colorSpec object - to make it ready for colorimetric calculations}
\alias{linearize}
\alias{linearize.colorSpec}

\description{
  linearize spectra and return modified object
}
\usage{
\S3method{linearize}{colorSpec}( x )
}

\arguments{
  \item{x}{a \bold{colorSpec} object}
}
    
\value{
  \code{linearize} returns a \bold{colorSpec} object with linear quantities.
}

\details{
If the \code{\link{quantity}(x)} is not \code{'absorbance'}
then \code{x} is returned unchanged.

If the \code{\link{quantity}(x)} is \code{'absorbance'}
then absorbance is converted to transmittance using
\deqn{ transmittance = 10^{-absorbance} }
Surprisingly, there does not seem to be a similar logarithmic version of reflectance.
Plots with log(responsivity) is somewhat common, but does not seem to have a separate name.
I have not seen log(radiometric power).
}

\seealso{
\code{\link{quantity}}
}


\keyword{colorSpec}

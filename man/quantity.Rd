% File quantity.Rd

\name{quantity}
\title{quantity of a colorSpec object}
\alias{quantity}
\alias{quantity<-}
\alias{type}
\alias{quantity.colorSpec}
\alias{quantity<-.colorSpec}
\alias{type.colorSpec}


\description{
  Retrieve or set the quantity of a \bold{colorSpec} object.
}
\usage{
\S3method{quantity}{colorSpec}(x)
quantity(x) <- value

\S3method{type}{colorSpec}(x)

}

\arguments{
\item{x}{a \bold{colorSpec} \R object}

\item{value}{a valid \code{quantity} string; see Details.}
}

\details{
There are 4 valid \code{types}, which are further divided into 14 valid \code{quantities}.
All of these are strings:

For \code{type='light'}\cr
\code{quantity = } \code{'power'} (radiometric), or \code{'photons/sec'} (actinometric)

For \code{type='responsivity.light'}\cr
\code{quantity = } \code{'power->electrical'}, \code{'power->neural'}, \code{'power->action'},\cr 
\code{'photons->electrical'}, \code{'photons->neural'}, or \code{'photons->action'} \cr

For \code{type='material'}\cr
\code{quantity = } \code{'reflectance'}, \code{'transmittance'}, or \code{'absorbance'} \cr

For \code{type='responsivity.material'}\cr
\code{quantity = } \code{'material->electrical'}, \code{'material->neural'}, or \code{'material->action'} \cr

For more explanation and examples of all these see 
\code{colorSpec-guide.pdf} in the doc directory.  
Go to the main index page and look for 
\bold{User guides, package vignettes and other documentation}.

}


\value{
  \code{quantity} returns the quantity of \code{x}
  
  \code{type} returns the type of \code{x}, which depends only on the \code{quantity}.
}

\note{
Although \code{quantity} is intended to represent a physical quantity,
the units are left arbitrary in most cases.
Exceptions are \code{reflectance}, \code{transmittance}, and \code{absorbance}
which are dimensionless.

Changing the \code{quantity} does not change the underlying numbers.
For example, changing \code{photons/sec} to \code{power} does not do numerical conversion.
For this conversion, see \code{\link{radiometric}}.

Similarly, see \code{\link{linearize}} for conversion from \code{absorbance} to \code{transmittance}.

}

\seealso{
  \code{\link{colorSpec}}; \code{\link{radiometric}}; \code{\link{linearize}}
}

\keyword{colorSpec}


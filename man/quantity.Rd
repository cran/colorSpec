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

\S3method{quantity}{colorSpec}(x) <- value

\S3method{type}{colorSpec}(x)

}

\arguments{
\item{x}{a \bold{colorSpec} \R object}

\item{value}{a valid \code{quantity} string; see Details.}
}

\details{
There are 4 valid \code{types}, which are further divided into 14 valid quantities.
All of these are strings:

For \code{type='light'}\cr
\code{quantity = } \code{'energy'} (radiometric), or \code{'photons'} (actinometric)

For \code{type='responsivity.light'}\cr
\code{quantity = } \code{'energy->electrical'}, \code{'energy->neural'}, \code{'energy->action'},\cr 
\code{'photons->electrical'}, \code{'photons->neural'}, or \code{'photons->action'} \cr

For \code{type='material'}\cr
\code{quantity = } \code{'reflectance'}, \code{'transmittance'}, or \code{'absorbance'} \cr

For \code{type='responsivity.material'}\cr
\code{quantity = } \code{'material->electrical'}, \code{'material->neural'}, or \code{'material->action'} \cr

Also see the \bold{colorSpec User Guide}.
}


\value{
  \code{quantity()} returns the quantity of \code{x}
  
  \code{type()} returns the type of \code{x}, which depends only on the \code{quantity}.
}

\note{
The \bold{colorSpec} \code{quantity} is more general than the physical SI quantity;
for example \code{quantity='energy'} really includes 10 distinct SI quantities
and maybe more.
The unit is left arbitrary in most cases.
Exceptions are \code{reflectance}, \code{transmittance}, and \code{absorbance}
which are dimensionless.

Changing the \code{quantity} should only be done if one knows what one is doing.
It does not change the underlying numbers.
For example, changing \code{photons} to \code{energy} does not do numerical conversion.
For this specific conversion, see \code{\link{radiometric}()}.

Similarly, see \code{\link{linearize}()} for conversion
from \code{absorbance} to \code{transmittance}.

}

\seealso{
\code{\link{colorSpec}},
\code{\link{radiometric}},
\code{\link{linearize}}
}

\keyword{colorSpec}


% File mean.Rd

\name{mean}
\alias{mean.colorSpec}
\alias{mean}

\title{calculate mean of multiple spectra}

\description{
compute mean of all spectra in a \bold{colorSpec} object
}

\usage{
\S3method{mean}{colorSpec}( x, ... )

}


\arguments{
\item{x}{a colorSpec object}

\item{...}{further arguments ignored}

}

\value{
a \bold{colorSpec} object with single spectrum = average of all spectra in \bold{colorSpec}.
}

\details{
This function might be useful when capturing many spectra on a spectrometer
and averaging them to reduce noise.

}

\keyword{colorSpec}

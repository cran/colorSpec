% File powerSpectra.Rd

\name{lightResponsivitySpectra}

\title{compute standard light responsivity spectra}
\alias{lightResponsivitySpectra}
\alias{erythemalSpectrum}




\usage{
erythemalSpectrum( wavelength=250:400 )

}

\arguments{
  \item{wavelength}{a vector of wavelengths, in nm}
}

\description{
Some action spectra standards are defined by simple equations;
the erythemal spectrum for human sunburn is one of them.

}

    
\value{
  For \code{erythemalSpectrum()} \cr
  A \bold{colorSpec} object with \code{\link{quantity}} equal to \code{'power->action'}.
  The responsivity is 0 for \eqn{\lambda} > 400 nm, so this putting this spectrum in the category
  of human vision is a bit of a stretch.
    
}

\details{
This erythemal spectrum is defined in 4 pieces: 
\eqn{\lambda<=298 }, 
\eqn{ 298<=\lambda<=328 },
\eqn{ 328<=\lambda<=400 }, and 
\eqn{ 400<\lambda}.  All the units are nm.
The spectrum is used in the definition of the international standard \bold{UV Index}.
}

\references{
McKinlay, A.F., and B.L. Diffey.
A reference action spectrum for ultraviolet induced erythema in human skin. 
CIE Res. Note, 6(1), 17-22. (1987)

}

\source{
http://en.wikipedia.org/wiki/Ultraviolet_index
}

\seealso{
\link{daylight},  \code{\link{quantity}},
\link{materialSpectra}, \link{lightSpectra}
}

\keyword{light}

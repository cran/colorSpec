% File wavelength.Rd


\name{wavelength}
\title{wavelength vector of a colorSpec object}
\alias{wavelength}
\alias{wavelength<-}
\alias{numWavelengths}
\alias{is.regular}
\alias{step.wl}
\alias{wavelength.colorSpec}
\alias{wavelength<-.colorSpec}
\alias{numWavelengths.colorSpec}
\alias{is.regular.colorSpec}
\alias{step.wl.colorSpec}

\description{
  Retrieve or set the wavelengths of a \bold{colorSpec} object.  
  Retrieve the number of wavelengths, and whether the wavelength sequence is regular.
}
\usage{
\S3method{wavelength}{colorSpec}(x)
wavelength(x) <- value

\S3method{numWavelengths}{colorSpec}(x)

\S3method{is.regular}{colorSpec}(x)
step.wl(x)

}

\arguments{
  \item{x}{a \bold{colorSpec} \R object}
  
  \item{value}{a numeric vector with length equal to the number of wavelengths in \code{x}.  
                The wavelengths must be increasing. The unit must be nanometers.}
}

\details{
If the \code{organization} of \code{x} is \code{'df.col'}, then \code{x} is a \code{data.frame}
and the wavelength vector is stored in the first column of \code{x}.\cr
Otherwise, the wavelength vector is stored as \code{attr(x,'wavelength')}.

}

\value{
  \code{wavelength} returns a numeric vector with the wavelength of the spectra, in nanometers.
  
  \code{numWavelengths(x)} is the same as \code{length(wavelength(x))} but much more efficient.
  
  \code{is.regular} returns TRUE or FALSE, depending on whether the step between consecutive
  wavelengths is a constant. A truncation error of 1.e-6 nm is tolerated here.  
  For example, the X-Rite ColorMunki spectrometer in hi-res mode has a step of 3.33333nm, 
  and it is considered regular.
  
  \code{step.wl} returns the \emph{mean} step in nm, whether the wavelengths are regular or not.
}

\seealso{
\code{\link{colorSpec}} 
}

\keyword{colorSpec}


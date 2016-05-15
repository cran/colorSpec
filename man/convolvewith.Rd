% File convolvewith.Rd

\name{convolvewith}
\title{Convolve each spectrum in a colorSpec object with a kernel}
\alias{convolvewith}
\alias{convolvewith.colorSpec}


\description{
This function convolves each spectrum in a colorSpec object with a kernel of odd length.
Its primary purpose is to correct raw spectrometer data
(with positive instrumental bandwidth) to have bandwidth=0.
Two popular correction kernels for this, with lengths 3 and 5,
are built-in options, see \bold{Details}.

}

\usage{
\S3method{convolvewith}{colorSpec}( x, coeff )
}

\arguments{
\item{x}{a \bold{colorSpec} object with N wavelengths}

\item{coeff}{
a convolution kernel of odd length. The center entry of this vector is taken
as index 0 in the convolution.\cr
\code{coeff} can also be the string \code{'SS3'} which means to apply the Stearns&Stearns
bandwidth correction kernel \code{coeff=c(-1,14,-1)/12}, see \bold{Details}.\cr
\code{coeff} can also be the string \code{'SS5'} which means to apply the Stearns&Stearns 
bandwidth correction kernel \code{coeff=c(1,-12,120,-12,1)/98}, see \bold{Details}.  
}
}
    
\value{
a \bold{colorSpec} object with the same dimensions,
\code{\link{wavelength}}, \code{\link{quantity}}, and \code{\link{organization}}.
If \code{coeff} is invalid it is an ERROR and
\code{convolvewith} returns \code{NULL}.
}


\details{
The built-in kernels, \code{'SS3'} and  \code{'SS5'}, were derived
by \cite{Stearns & Stearns} under specific hypotheses on the
spectrometer profile, bandpass, and pitch; see \bold{References}.\cr
Missing values at both ends are filled by copying from the nearest valid value.
}

\references{
Stearns, E.I., Stearns R.E.
An example of a method for correction radiance data for bandpass error.
\emph{Color Research and Application}.
13-4. 257-259. 1988.

Schanda, Janos.
CIE Colorimetry, in
\emph{Colorimetry: Understanding the CIE System}. 
Wiley Interscience. 2007.
p. 124.

Oleari, Claudio, Gabriele Simone.
\emph{Standard Colorimetry: Definitions, Algorithms and Software}.
John Wiley. 2016.
p. 309.
}


\seealso{
  \code{\link{quantity}},
  \code{\link{wavelength}},
  \code{\link{linearize}},  
  \code{\link{applyspec}},    
  \code{\link{organization}}
}


\keyword{colorSpec}

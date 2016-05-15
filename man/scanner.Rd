% File scanner.Rd

\name{scanner}
\alias{scanner.ACES}
\docType{data}

\title{
standard RGB scanners
}

\description{

\code{scanner.ACES} is an RGB responder to material; 
an ACES/SMPTE standard for scanning RGB film.
The 3 spectra are defined from 368 to 728 nm, at 2nm intervals.
}

\format{
A \bold{colorSpec} object with \code{\link{quantity}} equal to \code{'material->electrical'} and 3 spectra:
\code{r}, \code{g}, and \code{b}.
}

\details{
The responsivities have been scaled (by \code{\link{calibrate}}) so the response to the \emph{perfect transmitting filter} (PTF) is RGB=(1,1,1).

}


\references{
Technical Bulletin
TB-2014-005.
Informative Notes on
SMPTE ST 2065-2 - Academy Printing Density (APD).
Spectral Responsivities, Reference Measurement Device and Spectral Calculation.

SMPTE ST 2065-3 Academy Density Exchange Encoding (ADX).
Encoding Academy Printing Density (APD) Values.

The Academy of Motion Picture Arts and Sciences.
Science and Technology Council.
Academy Color Encoding System (ACES) Project Committee.
Version 1.0 December 19, 2014.
Annex A		Spectral Responsivities.
}


\seealso{
\code{\link{quantity}}, \code{\link{calibrate}}

}

\examples{
#   compute response of ACES scanner to the Hoya filters
product( Hoya, scanner.ACES, wave='auto' )
}

\keyword{datasets}

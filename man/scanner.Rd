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

##                  R            G          B
## R-60   0.902447043 2.022522e-05 0.00000000
## G-533  0.038450857 4.900983e-01 0.05431134
## B-440  0.008466317 1.686241e-02 0.42863320
## LB-120 0.184408941 3.264111e-01 0.53492533
}

\keyword{datasets}

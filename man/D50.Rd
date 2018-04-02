% File D50.Rd
\encoding{UTF-8}

\name{D50}
\alias{D50}
\alias{D50.5nm}
\docType{data}

\title{
Standard Illuminant D50 (1964)
}
\description{
\tabular{ll}{
\code{D50.5nm} \tab standard Illuminant D50, from 300 to 830 nm at 5 nm intervals.
}}

\format{
A \bold{colorSpec} object organized as a vector, with 107 data points 
and \code{\link{specnames}} equal to \code{'D50'}.
}

\details{
This spectrum is not read from a table from a CIE publication.
It is computed using the function \code{\link{daylightSpectra}()} by following
the special CIE recipe given in the \bold{References}.
The temperature is set to (14388/14380) * 5000 = 5002.781 Kelvin.
The coefficients of the daylight components \eqn{S_0, S_1}, and \eqn{S_2} are
rounded to 3 decimal places.
This linear combination is computed at 10nm intervals and then linearly interpolated to 5nm intervals.
The result is normalized to value 1 at 560nm (instead of the usual 100),
and finally rounded to 5 decimal places.  See \bold{Examples}.
}


\references{
Günther Wyszecki and W.S. Stiles.
\bold{Color Science : Concepts and Methods, Quantitative Data and Formulae}.
Second Edition.
Wiley-Interscience. 1982.
Table I(3.3.4)  pp. 754-758


CIE 15: Technical Report: Colorimetry, 3rd edition.
CIE 15:2004.
Table T.1, pp 30-32, and Note 5 on page 69.

Schanda, Janos.
CIE Colorimetry, in
\emph{Colorimetry: Understanding the CIE System}. 
Wiley Interscience. 2007.
p. 42.
}

\seealso{
\code{\link{ABC}} ,
\code{\link{D65}} ,
\code{\link{daylightSpectra}}
}

\examples{
#   the CIE recipe for computing D50.5nm
correction  = 14388 / 14380     # note 5, page 69 in CIE 15:2004
D50.10nm    = daylightSpectra( correction*5000, wavelength=seq(300,830,by=10), roundMs=TRUE )
D50.5nm     = resample( D50.10nm, seq(300,830,by=5), method='linear' )
D50.5nm     = round( D50.5nm, 5 )

summary( D50.5nm )
white.point = product( D50.5nm, xyz1931.1nm, wave='auto' )
}

\keyword{datasets}
\keyword{light}

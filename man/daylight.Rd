\encoding{UTF-8}

\name{daylight}
\alias{daylight1964}
\alias{daylight2013}
\alias{daylight}
\docType{data}
\title{
Standard Daylight Components
}
\description{
\tabular{ll}{
\code{daylight1964} \tab spectral components \eqn{S_0, S_1, S_2}; from 300 to 830 nm at 5 nm intervals \cr
\code{daylight2013} \tab smoothed spectral components  \eqn{S_0, S_1, S_2}; from 300 to 830 nm at 1 nm intervals \cr
}
}

\format{
Each is a \bold{colorSpec} object organized as a matrix with 3 columns
  \tabular{ll}{
    {\code{S0}} \tab {component 0, the mean power spectrum} \cr
    {\code{S1}} \tab {component 1, the 1st characteristic spectrum} \cr
    {\code{S2}} \tab {component 2, the 2nd characteristic spectrum} \cr    
  }
}

\source{
http://www.cie.co.at/publ/abst/datatables15_2004/CIE_sel_colorimetric_tables.xls

http://vision.vein.hu/~schanda/CIE\%20TC1-74/
}

\references{
Günther Wyszecki and W.S. Stiles.
\bold{Color Science : Concepts and Methods, Quantitative Data and Formulae}.
Second Edition.
Wiley-Interscience. 1982.
Table V(3.3.4)  p. 762.

Smoothing spectral power distribution of daylights.
Zsolt Kosztyan and Janos Schanda.
Color Research & Application.
Volume 38, Issue 5, pages 316-321, October 2013.

CIE 15: Technical Report: Colorimetry, 3rd edition.
CIE 15:2004.
Table T.2, pp 33-34

JUDD, D.B., MACADAM, D.L. and WYSZECKI, G.,
with the collaboration of BUDDE, H.W, CONDIT, H.R, HENDERSON, S.T, and SIMONDS, J.L. 
Spectral distribution of typical daylight as a function of correlated color temperature. 
J Opt. Soc. Am. 54, 1031-1040, 1964.

Zsolt Kosztyan and Janos Schanda.
Smoothing spectral power distribution of daylights.
Color Research & Application.
Volume 38, Issue 5, pages 316-321, October 2013.
}


\seealso{
\code{\link{D65}},
\code{\link{D50}}, 
\code{\link{daylightSpectra}} 

}

\examples{
summary( daylight1964 )
day1964 = daylightSpectra( c(5000,6500), comp=daylight1964 )
day2013 = daylightSpectra( c(5000,6500), comp=daylight2013 )

plot( day1964, col='black' )
plot( day2013, col='black', add=TRUE )

}
\keyword{datasets}
\keyword{light}

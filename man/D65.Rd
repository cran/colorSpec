% File D65.Rd
\encoding{UTF-8}

\name{D65}
\alias{D65}
\alias{D65.1nm}
\alias{D65.5nm}
\docType{data}
\title{
Standard Illuminant D65 (1964)
}
\description{
\tabular{ll}{
\code{D65.1nm} \tab standard Illuminant D65, 300 to 830 nm at 1 nm intervals \cr
\code{D65.5nm} \tab standard Illuminant D65, 380 to 780 nm at 5 nm intervals \cr
}
}


\format{
Each is a \bold{colorSpec} object organized as a vector, 
with \code{\link{specnames}} equal to \code{'D65'}.
}

\details{
Both of these have been divided by 100, 
to make the values at 560nm equal to 1 instead of 100.
}

\source{
\url{http://www.cvrl.org}

}

\references{
Günther Wyszecki and W.S. Stiles.
\bold{Color Science : Concepts and Methods, Quantitative Data and Formulae}.
Second Edition.
Wiley-Interscience. 1982.
Table I(3.3.4)  pp. 754-758


ASTM E 308-01.
Standard Practice for Computing the Colors of Objects by Using the CIE System.
Table 3. pages 3-4.

}


\seealso{
\code{\link{ABC}}, \code{\link{D50}}, \code{\link{daylightSpectra}} , \code{\link{daylight}} 
}

\examples{
summary( D65.1nm )
white.point = product( D65.1nm, xyz1931.1nm, wave='auto' )
}
\keyword{datasets}
\keyword{light}

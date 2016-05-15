\encoding{UTF-8}

\name{ABC}
\alias{A.1nm}
\alias{B.5nm}
\alias{C.5nm}
\alias{ABC}
\docType{data}

\title{
Standard Illuminants A, B, and C (1931)
}
\description{
\tabular{ll}{
\code{A.1nm} \tab standard Illuminant A, 360 to 780 nm at 1 nm intervals \cr
\code{B.5nm} \tab standard Illuminant B, 320 to 780 nm at 5 nm intervals \cr
\code{C.5nm} \tab standard Illuminant C, 320 to 780 nm at 5 nm intervals
}
}

\format{
Each is a \bold{colorSpec} object organized as a vector, with \code{\link{quantity}} equal to \code{"power".}

}

\details{
All of these have been divided by 100, to make the values at 560nm near 1 instead of 100.
}

\source{
http://www.cvrl.org

}

\references{
Günther Wyszecki and W. S. Stiles.
\bold{Color Science : Concepts and Methods, Quantitative Data and Formulae}.
Second Edition.
Wiley-Interscience. 1982.

\tabular{ll}{
A \tab Table I(3.3.4)  pp. 754-758. \cr
B \tab Table II(3.3.4)  pp. 759. \cr
C \tab Table II(3.3.4)  pp. 759. \cr
}
}

\seealso{
\code{\link{D50}}
\code{\link{D65}} 
}

\examples{
summary(xyz1931.1nm)
white.point = product( D65.1nm, xyz1931.1nm, wave='auto' )
}
\keyword{datasets}
\keyword{light}

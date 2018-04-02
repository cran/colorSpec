\encoding{UTF-8}

\name{xyz1931}
\alias{xyz1931.1nm}
\alias{xyz1931.5nm}
\alias{xyz1931}
\docType{data}

\title{
CIE Color Matching Functions - 2-degree (1931)
}
\description{
\tabular{ll}{
\code{xyz1931.1nm} \tab the 1931 2° functions from 360 to 830 nm, at 1nm intervals \cr
\code{xyz1931.5nm} \tab the 1931 2° functions from 380 to 780 nm, at 5nm intervals \cr
}
}

\format{
  Each is a \bold{colorSpec} object organized as a matrix with 3 variables.
  \tabular{ll}{
    {\code{x}} \tab {the x-bar responsivity function} \cr
    {\code{y}} \tab {the y-bar responsivity function} \cr
    {\code{z}} \tab {the z-bar responsivity function} \cr
  }
}

\references{
Günther Wyszecki and W.S. Stiles.
\bold{Color Science : Concepts and Methods, Quantitative Data and Formulae}.
Second Edition.
Wiley-Interscience. 1982.
Table I(3.3.1).
pp. 723-735.

ASTM E 308 - 01.
Standard Practice for Computing the Colors of Objects by Using the CIE System.
Table 1 
}

\source{
\url{http://www.cvrl.org}
}

\seealso{
\code{\link{xyz1964}} 
}

\examples{
summary(xyz1931.1nm)
white.point = product( D65.1nm, xyz1931.1nm, wave='auto' )
}

\keyword{datasets}
\keyword{eyes}

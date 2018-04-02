\name{Fluorescents}
\alias{Fs.5nm}
\alias{Fluorescents}
\docType{data}

\title{
Standard series F Illuminants F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, and F12
}
\description{
\code{Fs.5nm} contains 12 CIE Fluorescent Illuminants, from 380 to 780 nm, at 5nm intervals.
}

\format{
\code{Fs.5nm} is a \bold{colorSpec} object with 12 spectra.
It is organized as a data frame  with \code{\link{quantity}} equal to \code{"energy".}
}

\note{
The series F illuminants do not seem to be normalized in a consistent way.
}

\source{
\url{http://www.rit-mcsl.org/UsefulData/Fluorescents.htm}
}

\seealso{
\code{\link{ABC}},  \code{\link{D50}},   \code{\link{D65}}  
}

\examples{
#   plot only F4
plot( subset(Fs.5nm,"F4") )
}

\keyword{datasets}

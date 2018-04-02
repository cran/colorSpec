% File F96T12.Rd

\name{F96T12}
\alias{F96T12}
\docType{data}

\title{
Photon Irradiance of F96T12 Fluorescent Bulb
}
\description{
\code{F96T12}\cr
Sylvania F96T12 CW/VHO 215-Watt fluorescent bulb photon irradiance, 
measured with a LI-COR LI-1800 spectroradiometer,
from 300 to 900 nm at 1 nm intervals.

}

\format{
A \bold{colorSpec} object organized as a vector, with 601 data points 
and \code{\link{specnames}} equal to \code{'F96T12'}.
}

\details{
The unit is (\eqn{\mu}mole of photons)\eqn{ *sec^{-1}*m^{-2}*nm^{-1}}.
}


\source{
Pedro J. Aphalo.
\url{http://www.mv.helsinki.fi/aphalo/photobio/lamps.html}
}

\seealso{
\code{\link{ABC}} ,
\code{\link{D65}} ,
\code{\link{daylightSpectra}}
}

\examples{
sum( F96T12 ) 
# [1] 320.1132  photon irradiance, (micromoles of photons)*m^{-2}

sum( radiometric(F96T12) )
# [1] 68.91819  irradiance, watts*m^{-2}
}

\keyword{datasets}
\keyword{light}

% File radiometric.Rd

\name{radiometric}
\title{force a colorSpec object to be radiometric}
\alias{radiometric}
\alias{radiometric.colorSpec}
\alias{actinometric}

\description{
convert a \bold{colorSpec} object to have quantity that is radiometric (power of photons) - to make it ready for colorimetric calculations
}
\usage{
\S3method{radiometric}{colorSpec}( x )
}

\arguments{
  \item{x}{a \bold{colorSpec} object}

}
    
\value{
  \code{radiometric} returns a \bold{colorSpec} object with \code{\link{quantity}} that is
  radiometric (power of photons) and not actinometric (number of photons).
}


\details{
If the \code{\link{quantity}} starts with \code{'power'}
then \code{x} is returned unchanged.
And if the \code{\link{type}} of \code{x} is \code{'material'} or \code{'responsivity.material'}
then \code{x} is returned unchanged.


Otherwise, the \code{\link{quantity}} starts with \code{'photons'} 
so the \code{quantity}  of \code{x} is actinometric (proportional to number of photons/sec).

If the \code{\link{type}} is \code{'light'} then
the most common actinometric unit of photon flux is (\eqn{\mu}mole of photons)/sec.
The conversion equation is:
\deqn{ p = x * (10^{-6} * N_A * h * c / \lambda) }
where \eqn{p} is proportional to the power of photons,
\eqn{x} is the photon flux,
\eqn{N_A} is Avogadro's constant,
\eqn{h} is Planck's constant, \eqn{c} is the speed of light,
and \eqn{\lambda} is the wavelength in meters.
The output power unit is watts.\cr
If the unit of \code{x} is not (\eqn{\mu}mole of photons)/sec, 
then the output should be scaled appropriately.
For example, if the photon flux is exaphotons/sec,
then divide the output by 0.6022.

If the \code{\link{type}} is \code{'responsivity.light'},
then the most common actinometric unit of responsivity to light is quantum efficiency (QE).
The conversion equation is:
\deqn{ r = x * \lambda * e / (h * c) }
where \eqn{r} is the responsivity,
\eqn{x} is the quantum efficiency,
and \eqn{e} is the charge of an electron (in C).
The output responsivity unit is amps/watt (A/W) or coulombs/joule (C/J).\cr
If the unit of \code{x} is not quantum efficiency, 
then the output should be scaled appropriately.
}

\source{
Wikipedia.
\bold{Photon counting}.
http://en.wikipedia.org/wiki/Photon_counting
}

\seealso{
  \code{\link{quantity}}, \code{\link{type}}, \code{\link{F96T12}}
}

\examples{
sum( F96T12 ) 
# [1] 320.1132  photon irradiance, (micromoles of photons)*m^{-2}

sum( radiometric(F96T12) )
# [1] 68.91819  irradiance, watts*m^{-2}
}

\keyword{light}

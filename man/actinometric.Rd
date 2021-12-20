% File radiometric.Rd

\name{actinometric}
\title{convert a colorSpec object to be actinometric}
\alias{actinometric.colorSpec}
\alias{actinometric}
\alias{is.actinometric.colorSpec}
\alias{is.actinometric}

\description{
Convert a radiometric \bold{colorSpec} object to have
quantity that is actinometric (number of photons).
Test an object for whether it is actinometric.
}
\usage{
\S3method{actinometric}{colorSpec}( x, multiplier=1, warn=FALSE )

\S3method{is.actinometric}{colorSpec}( x )
}

\arguments{
\item{x}{a \bold{colorSpec} object}
\item{multiplier}{a scalar which is multiplied by the output, and intended for unit conversion}
\item{warn}{if \code{TRUE} and a conversion actually takes place, the a \code{WARN} message is issued.
This makes the user aware of the conversion, so units can be verified.  This can be useful when \code{actinometric()} is called from another \bold{colorSpec} function.}
}
    
\value{
\code{actinometric()} returns a \bold{colorSpec} object
with \code{\link{quantity}} that is
actinometric (photon-based) and not radiometric (energy-based).
If \code{type(x)} is a material type
(\code{'material'} or \code{'responsivity.material'})
then \code{x} is returned unchanged.

If \code{quantity(x)} starts with \code{'photons'},
then \code{is.actinometric()} returns \code{TRUE}, and otherwise \code{FALSE}.
}


\details{
If the \code{\link{quantity}} of \code{x} does not start with \code{'energy'}
then the quantity is not radiometric and so \code{x} is returned unchanged.
Otherwise \code{x} is radiometric (energy-based), and must be converted.

If \code{\link{type}(x)} is \code{'light'} then
the most common radiometric energy unit is joule.\cr
The conversion equation is:
\deqn{ Q = E * \lambda * 10^6 / (N_A * h * c)}
wher \eqn{Q} is the photon count,
\eqn{E} is the energy of the photons,
\eqn{N_A} is Avogadro's constant,
\eqn{h} is Planck's constant, \eqn{c} is the speed of light,
and \eqn{\lambda} is the wavelength.
The output unit of photon count is
(\eqn{\mu}mole of photons) = (\eqn{6.02214 * 10^{17}} photons).
If a different unit for \code{Q} is desired, 
then the output should be scaled appropriately.
For example, if the desired unit of photon count is exaphotons,
then set \code{multiplier=0.602214}.

If the \code{\link{quantity}(x)} is \code{'energy->electrical'},
then the most common radiometric unit of responsivity to light
is coulombs/joule (C/J) or  amps/watt (A/W).
The conversion equation is:
\deqn{ QE = R_e * ((h * c)/e) / \lambda }
where \eqn{QE} is the quantum efficiency,
\eqn{R_e} is the energy-based responsivity,
and \eqn{e} is the charge of an electron (in C).\cr
If the unit of \code{x} is not C/J, 
then \code{multiplier} should be set appropriately.

If the \code{\link{quantity}(x)} is 
\code{'energy->neural'} or  \code{'energy->action'},
the most common radiometric unit of energy is joule (J).

The conversion equation is:
\deqn{ R_p = R_e  * 10^{-6} * ( N_A * h * c) / \lambda }
where \eqn{R_p} is the photon-based responsivity,
and \eqn{R_e} is the energy-based responsivity,
The output unit of photon count is
(\eqn{\mu}mole of photons) = (\eqn{6.02214 * 10^{17}} photons).
This essentially the reciprocal of the first conversion equation.

The argument \code{multiplier} is applied to the right side of all the above
conversion equations.
}


\note{
To log the executed conversion equation,
execute \code{cs.options(loglevel='INFO')}.
}

\source{
Wikipedia.
\bold{Photon counting}.
\url{https://en.wikipedia.org/wiki/Photon_counting}
}

\seealso{
\code{\link{quantity}}, 
\code{\link{type}},
\code{\link{cs.options}},
\code{\link{radiometric}}
}

\examples{
colSums( solar.irradiance ) # step size is 1nm, from 280 to 1000 nm. organized as a matrix
# AirMass.0  GlobalTilt AirMass.1.5 
#  944.5458    740.3220    649.7749  # irradiance, watts*m^{-2}


colSums( actinometric(solar.irradiance) )
# AirMass.0  GlobalTilt AirMass.1.5 
#  4886.920    3947.761    3522.149  # photon irradiance, (umoles of photons)*sec^{-1}*m^{-2}

colSums( actinometric(solar.irradiance,multiplier=0.602214) )
# AirMass.0  GlobalTilt AirMass.1.5 
#  2942.972    2377.397    2121.088  # photon irradiance, exaphotons*sec^{-1}*m^{-2}

}

\keyword{light}

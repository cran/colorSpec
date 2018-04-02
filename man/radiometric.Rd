% File radiometric.Rd

\name{radiometric}
\title{convert a colorSpec object from actinometric to radiometric}
\alias{radiometric}
\alias{radiometric.colorSpec}
\alias{is.radiometric}
\alias{is.radiometric.colorSpec}

\description{
Convert a \bold{colorSpec} object to have quantity that is radiometric (energy of photons) - to prepare it for colorimetric calculations.
Test an object for whether it is radiometric.
}
\usage{
\S3method{radiometric}{colorSpec}( x, multiplier=1, warn=FALSE )

\S3method{is.radiometric}{colorSpec}( x )
}

\arguments{
  \item{x}{a \bold{colorSpec} object}
  \item{multiplier}{a scalar which is multiplied by the output, and intended for unit conversion}
\item{warn}{if \code{TRUE} and a conversion actually takes place, the a \code{WARN} message is issued.
This makes the user aware of the conversion, so units can be verified.  This can be useful when \code{radiometric()} is called from another \bold{colorSpec} function.}
}

    
\value{
\code{radiometric()} returns a \bold{colorSpec} object with
\code{\link{quantity}} that is
radiometric (energy-based) and not actinometric (photon-based).
If \code{type(x)} is a material type
(\code{'material'} or \code{'responsivity.material'})
then \code{x} is returned unchanged.
 
If \code{quantity(x)} starts with \code{'energy'},
then \code{is.radiometric()} returns \code{TRUE}, and otherwise \code{FALSE}.
}


\details{
If the \code{\link{quantity}} of \code{x} does not start with \code{'photons'}
then the quantity is not actinometric
and so \code{x} is returned unchanged.
Otherwise \code{x} is actinometric (photon-based).

If \code{\link{type}(x)} is \code{'light'} then
the most common actinometric unit of photon count is
(\eqn{\mu}mole of photons) = (\eqn{6.02214 x 10^{17}} photons).
The conversion equation is:
\deqn{ E = Q * 10^{-6} * N_A * h * c / \lambda }
where \eqn{E} is the energy of the photons,
\eqn{Q} is the photon count,
\eqn{N_A} is Avogadro's constant,
\eqn{h} is Planck's constant, \eqn{c} is the speed of light,
and \eqn{\lambda} is the wavelength in meters.
The output energy unit is joule.\cr
If the unit of \code{Q} is not (\eqn{\mu}mole of photons), 
then the output should be scaled appropriately.
For example, if the unit of photon count is exaphotons,
then set \code{multiplier=1/0.602214}.

If the \code{\link{quantity}(x)} is \code{'photons->electrical'},
then the most common actinometric unit of responsivity to light is quantum efficiency (QE).
The conversion equation is:
\deqn{ R_e = QE * \lambda * e / (h * c) }
where \eqn{R_e} is the energy-based responsivity,
\eqn{QE} is the quantum efficiency,
and \eqn{e} is the charge of an electron (in C).
The output responsivity unit is coulombs/joule (C/J) or  amps/watt (A/W).\cr
If the unit of \code{x} is not quantum efficiency, 
then \code{multiplier} should be set appropriately.

If the \code{\link{quantity}(x)} is 
\code{'photons->neural'} or  \code{'photons->action'},
the most common actinometric unit of photon count is
(\eqn{\mu}mole of photons) = (\eqn{6.02214 x 10^{17}} photons).
The conversion equation is:
\deqn{ R_e = R_p * \lambda * 10^6 / ( N_A * h * c) }
where \eqn{R_e} is the energy-based responsivity,
\eqn{R_p} is the photon-based responsivity.
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
\url{http://en.wikipedia.org/wiki/Photon_counting}
}

\seealso{
\code{\link{quantity}}, 
\code{\link{type}}, 
\code{\link{F96T12}},
\code{\link{cs.options}},
\code{\link{actinometric}}
}


\examples{
sum( F96T12 )    # the step size is 1nm, from 300 to 900nm
# [1] 320.1132  photon irradiance, (micromoles of photons)*m^{-2}*sec^{-1}

sum( radiometric(F96T12) )
# [1] 68.91819  irradiance, watts*m^{-2}
}

\keyword{light}

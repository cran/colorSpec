\name{solar.irradiance}

\alias{solar.irradiance}
\alias{atmosphere2003}

\docType{data}
\title{
Standard Solar Irradiance - Extraterrestrial and Terrestrial
}

\description{
\describe{
\item{\code{solar.irradiance}}{Three power spectra; from 280 to 1000 nm at 1 nm intervals. The unit is \eqn{W*m^{-2}*nm^{-1}}.}

\item{\code{atmosphere2003}}{a transmittance spectrum = the quotient of 2 spectra from \code{solar.irradiance}}
}
}

\format{
\code{solar.irradiance} is a
\bold{colorSpec} object with \code{quantity} equal to \code{'energy'} and with 3 spectra:

\describe{
\item{\code{AirMass.0}}{Extraterrestrial Radiation (solar spectrum at top of atmosphere) at mean Earth-Sun distance}
\item{\code{GlobalTilt}}{spectral radiation from solar disk plus sky diffuse
and diffuse reflected from ground on south facing surface tilted 37 deg from horizontal}
\item{\code{AirMass.1.5}}{the sum of Direct and Circumsolar irradiance, when the optical path is 1.5 times that of the sun at zenith, see \bold{Details}}
}
 
\code{atmosphere2003} is a \bold{colorSpec} object with \code{quantity} equal to \code{'transmittance'} and with 1 spectrum:
\describe{
\item{\code{AirMass.1.5}}{the quotient \code{AirMass.1.5} / \code{AirMass.0} from \code{solar.irradiance} }
}}

\details{
\bold{Direct} is Direct Normal Irradiance Nearly parallel (0.5 deg divergent cone) radiation
on surface with surface normal tracking (pointing to) the sun,
excluding scattered sky and reflected ground radiation.

\bold{Circumsolar} is Spectral irradiance within +/- 2.5 degree (5 degree diameter) field of view
centered on the 0.5 deg diameter solar disk, but excluding the radiation from the disk.
}

\note{
The reference spectra in ASTM G173-03 are designed for
Photovoltaic Performance Evaluation.
 
The original wavelength sequence in ASTM G173-03 is irregular.
The interval is 0.5 nanometer from 280 to 400 nm,
1 nm from 400 to 1700 nm,
an intermediate wavelength at 1702 nm,
and 5 nm from 1705 to 4000 nm.
To create the object \code{solar.irradiance} with a regular step size,
the original was resampled from 280 to 1000 nm at 1nm intervals.
}

\source{
Reference Solar Spectral Irradiance: ASTM G-173.
\url{http://rredc.nrel.gov/solar/spectra/am1.5/astmg173/astmg173.html}
}

\references{
ASTM G173-03 Reference Spectra Derived from SMARTS v. 2.9.2.\cr
Standard Tables for Reference Solar Spectral Irradiances: Direct Normal and Hemispherical on 37-deg Tilted Surface (2003)
}


\seealso{
\code{\link{D65}},
\code{\link{D50}}, 
\code{\link{daylightSpectra}},
\code{\link{resample}},
vignette \href{../doc/blueflame.pdf}{\bold{Blue Flame and Green Comet}}

}

\keyword{datasets}


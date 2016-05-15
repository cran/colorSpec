\name{solar.irradiance}

\alias{solar.irradiance}
\alias{atmosphere2003}

\docType{data}
\title{
Standard Solar Irradiance - Extraterrestrial and Terrestrial
}
\description{

\code{solar.irradiance}\cr
Three power spectra; from 280 to 1000 nm at 1 nm intervals. Units are \eqn{W*m^{-2}*nm^{-1}}. \cr

\code{atmosphere2003}\cr
a transmittance spectrum = a ratio of 2 spectra from \code{solar.irradiance} 
}

\format{
\code{solar.irradiance} is a
\bold{colorSpec} object with \code{quantity} equal to \code{'power'} and with 3 spectra:

\describe{
\item{\code{AirMass.0}}{Extraterrestrial Radiation (solar spectrum at top of atmosphere) at mean Earth-Sun distance}
\item{\code{GlobalTilt}}{spectral radiation from solar disk plus sky diffuse and diffuse reflected
                   from ground on south facing surface tilted 37 deg from horizontal}
\item{\code{AirMass.1.5}}{the sum of Direct and Circumsolar irradiance, see \bold{Details}}
}
 
\code{atmosphere2003} is a \bold{colorSpec} object with \code{quantity} equal to \code{'transmittance'} and with 1 spectrum:
\describe{
\item{\code{AirMass.1.5}}{the ratio of \code{AirMass.1.5}/\code{AirMass.0} from \code{solar.irradiance} }
}}

\details{
\bold{Direct} is Direct Normal Irradiance Nearly parallel (0.5 deg divergent cone) radiation
on surface with surface normal tracking (pointing to) the sun,
excluding scattered sky and reflected ground radiation.

\bold{Circumsolar} is Spectral irradiance within +/- 2.5 degree (5 degree diameter) field of view
centered on the 0.5 deg diameter solar disk, but excluding the radiation from the disk.
}

\source{
http://rredc.nrel.gov/solar/spectra/am1.5/astmg173/astmg173.html
}

\references{
ASTM G173-03 Reference Spectra Derived from SMARTS v. 2.9.2.\cr
Standard Tables for Reference Solar Spectral Irradiances: Direct Normal and Hemispherical on 37-deg Tilted Surface (2003)
}


\seealso{
\code{\link{D65}},
\code{\link{D50}}, 
\code{\link{daylightSpectra}},
vignette blueflame

}

\keyword{datasets}


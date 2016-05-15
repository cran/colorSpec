% File photometric.Rd

\name{photometric}
\title{convert illuminant spectra to photometric units}
\alias{photometric}
\alias{photometric.colorSpec}

\description{
Convert radiometric units of power to photometric units, 
using the standard photometric weighting curve.  
Actinometric units (number of photons) are converted to radiometric units (power of photons) on-the-fly.
}

\usage{
\S3method{photometric}{colorSpec}( x )
}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type} equal to \code{'light'}}

}
    
\value{
  \code{photometric} returns an Mx1 matrix, where M is the number of spectra in \code{x}.
  The \code{rownames} are \code{specnames(x)}, and the \code{colnames} are \code{'Y'}.\cr
  In case of ERROR it returns \code{NULL}.\cr
}

\note{
To get the right output quantity and units, 
the user must know the input quantity and units.
For example, pre-conversion from microwatts to watts might be necessary.\cr
The \emph{scoptic} weighting curve (for very low light levels) is not used.
}


\details{
The function computes the dot product of \code{x} with the 
photopic responsivity spectrum \code{y} from \code{xyz1931.1nm}.
This product is an Mx1 matrix, where M is the number of spectra in \code{x}.
It is then multiplied by the CIE-standard coefficient of 683 lumen/watt and returned.

The 5 radiometric inputs and photometric outputs are:
\tabular{lll}{
\bold{radiant flux power} [\eqn{watt}] \tab  --->  \tab  \bold{luminous flux} [\eqn{lumen}] \cr
\bold{irradiance}  [\eqn{watt*m^{-2}}] \tab  --->  \tab  \bold{illuminance} [\eqn{lumen*m^{-2} = lux}] \cr
\bold{radiant exitance}  [\eqn{watt*m^{-2}}] \tab  --->  \tab  \bold{luminous exitance} [\eqn{lumen*m^{-2} = lux}] \cr
\bold{radiant intensity}  [\eqn{watt*sr^{-1}}] \tab  --->  \tab  \bold{luminous intensity} [\eqn{lumen*sr^{-1} = candela}] \cr
\bold{radiance}  [\eqn{watt*sr^-1*m^{-2}}] \tab  --->  \tab  \bold{luminance} [\eqn{candela*m^{-2} =  nit}] \cr

}
}

\references{
Poynton, Charles.
Digital Video and HD - Algorithms and Interfaces.
Morgan Kaufmann.
Second Edition. 2012.
Appendix B, pp. 573-580.
}


\seealso{
\code{\link{quantity}}, \code{\link{type}}, \code{\link{radiometric}}
}


\examples{
photometric( solar.irradiance )  # units are watt*m^{-2}

#                     Y   # units are lux
# AirMass.0   133100.41
# GlobalTilt  109494.88
# AirMass.1.5  97142.25
}

\keyword{light}

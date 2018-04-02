% File photometric.Rd

\name{photometric}
\title{convert illuminant spectra to photometric units}
\alias{photometric}
\alias{photometric.colorSpec}

\description{
Convert radiometric units of power or energy to photometric units, 
using 4 standard photometric weighting curves.  
Actinometric units (number of photons) are converted to
radiometric units (energy of photons) on-the-fly.
}

\usage{
\S3method{photometric}{colorSpec}( x, photopic=683, scotopic=1700, multiplier=1 )
}

\arguments{
\item{x}{a \bold{colorSpec} object with \code{type} equal to \code{'light'}, and with M spectra}
\item{photopic}{the conversion factor for photopic vision, in lumen/watt.  
The CIE standard is 683, and another common value is 683.002.}
\item{scotopic}{the conversion factor for scotopic vision, in lumen/watt.  
The CIE standard is 1700, and another common value is 1700.06.}
\item{multiplier}{an conversion factor intended for conversion of units,
and applied to both photopic and scotopic vision.
For example if the input unit of \code{x} is \eqn{watt*sr^{-1}},
and the desired output unit is \eqn{candlepower}, then set
\code{multiplier=1/0.981}.}
}

\value{
\code{photometric()} returns an Mx4 matrix, where M is the number of spectra in \code{x}.
The \code{rownames} are \code{specnames(x)},
and the \code{colnames} are \code{specnames(luminsivity.1nm)}.\cr
In case of ERROR it returns \code{NULL}.\cr
}

\note{
To get the right output quantity and units, 
the user must know the input quantity and units.
If the units are different than those in the above list,
then set \code{multiplier} appropriately.\cr
It is up to the user to determine whether
\emph{photopic} or \emph{scotopic} vision (or neither) is appropriate.
The intermediate \emph{mesopic} vision is currently a subject of research by the CIE,
and might be added to this function in the future.
}


\details{
The function computes the product of \code{x} with \code{luminsivity.1nm}.
This product is an Mx4 matrix, where M is the number of spectra in \code{x}.
There are 3 columns for photopic vision, and 1 column for scotopic vision.
These columns are multiplied by the appropriate conversion factors
and the resulting Mx4 matrix is returned.

The 5 power-based input quantities and corresponding photometric outputs are:
\tabular{lll}{
\bold{radiant power} [\eqn{watt}] \tab  --->  \tab  \bold{luminous flux} [\eqn{lumen}] \cr
\bold{irradiance}  [\eqn{watt*m^{-2}}] \tab  --->  \tab  \bold{illuminance} [\eqn{lumen*m^{-2} = lux}] \cr
\bold{radiant exitance}  [\eqn{watt*m^{-2}}] \tab  --->  \tab  \bold{luminous exitance} [\eqn{lumen*m^{-2} = lux}] \cr
\bold{radiant intensity}  [\eqn{watt*sr^{-1}}] \tab  --->  \tab  \bold{luminous intensity} [\eqn{lumen*sr^{-1} = candela}] \cr
\bold{radiance}  [\eqn{watt*sr^-1*m^{-2}}] \tab  --->  \tab  \bold{luminance} [\eqn{candela*m^{-2} =  nit}] \cr
}

The 2 \emph{common} energy-based input quantities and corresponding photometric outputs are:
\tabular{lll}{
\bold{radiant energy} [\eqn{joule}] \tab  --->  \tab  \bold{luminous energy} [\eqn{talbot = lumen-second}] \cr
\bold{radiant exposure}  [\eqn{joule*m^{-2}}] \tab  --->  \tab  \bold{luminous exposure} [\eqn{talbot*m^{-2} = lux-second}] \cr
}
and there are 3 more obtained by integrating over time.
For example "time-integrated radiance" --->  "time integrated luminance".
But I have not been able to find names for these 3.
The \emph{talbot} is the unofficial name for a lumen-second.
}

\references{
Poynton, Charles.
\bold{Digital Video and HD - Algorithms and Interfaces}.
Morgan Kaufmann.
Second Edition. 2012.
Appendix B, pp. 573-580.
}


\seealso{
\code{\link{quantity}},
\code{\link{type}},
\code{\link{luminsivity.1nm}},
\code{\link{radiometric}}
}


\examples{
photometric( solar.irradiance )  # unit is watt*m^{-2}

#             photopic1924 scotopic1951 photopic1978 photopic2008  # units are lux
# AirMass.0      133100.41     313883.2    133843.65     140740.3
# GlobalTilt     109494.88     250051.5    110030.31     115650.0
# AirMass.1.5     97142.25     215837.1     97571.57     102513.7
}

\keyword{light}

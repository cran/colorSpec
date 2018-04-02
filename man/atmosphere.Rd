% File materialSpectra.Rd

\name{atmosphere}
\alias{atmosTransmittance}

\title{atmospheric transmittance along a horizontal path}

\description{
Calculate transmittance along a horizontal optical path in the atmosphere,
as a function of length (distance) and the molecular and aerosol properties.
Because the path is horizontal, the atmospheric properties are
assumed to be constant on the path.
Only molecular scattering is considered.
There is no modeling of molecular absorption;
for visible wavelengths this is reasonable.
}

\usage{
atmosTransmittance( distance, wavelength=380:720, 
                    molecules=list(N=2.547305e25,n0=1.000293),
                    aerosols=list(metrange=25000,alpha=0.8,beta=0.0001) )
}

\arguments{
\item{distance}{the length of the optical path, in meters.
It can also be a numeric vector of lengths.}

\item{wavelength}{a vector of wavelengths, in nm, for the transmittance calculations}

\item{molecules}{a list of molecular properties, see \bold{Details}.
If this is \code{NULL}, then the molecular transmittance is identically 1.}
        
\item{aerosols}{a list of aerosol properties, see \bold{Details}.
If this is \code{NULL}, then the aerosol transmittance is identically 1.}
}

\details{
The list \code{molecules} has 2 parameters that describe the molecules in the atmosphere.
\code{N} is the molecular density of the atmosphere at sea level,
in \eqn{molecules/meter^3}.
The given default is the density at sea level.
\code{n0} is the refractive index of pure molecular air (with no aerosols).
For the molecular attenuation,
the standard model for Rayleigh scattering is used,
and there is no modeling of molecular absorption.

The list \code{aerosols} has 3 parameters that describe the aerosols in the atmosphere.
The standard Angstrom aerosol attenuation model is:
\deqn{attenuation(\lambda) = \beta * (\lambda/\lambda_0)^{-\alpha}}
\eqn{\alpha} is the Angstrom exponent, and is dimensionless.
\eqn{attenuation} and \eqn{\beta} have unit \eqn{m^{-1}}.
And \eqn{\lambda_0}=550nm.

\code{metrange} is the \emph{Meteorological Range} of the atmosphere in meters,
as defined by \cite{Koschmieder}.
This is the distance at which the transmittance=0.02 at \eqn{\lambda_0}.
If \code{metrange} is not \code{NULL} (the default is 25000)
then both \eqn{\alpha} and \eqn{\beta} are calculated to achieve
this desired \code{metrange}, and the supplied \eqn{\alpha} and \eqn{\beta}
are ignored.
\eqn{\alpha} is calculated from \code{metrange} using the \cite{Kruse} model,
see \bold{Note}.
\eqn{\beta} is calculated so that the product of
molecular and aerosol transmittance yields the desired \code{metrange}.
In fact:
\deqn{\beta = -\mu_0 - log(0.02) / V_r}
where \eqn{\mu_0} is the molecular attenuation at \eqn{\lambda_0},
and \eqn{V_r} is the meteorological range.
For a log message with the calculated values,
execute \code{cs.options(loglevel='INFO')} before calling \code{atmosTransmittance()}.
}
    
\value{
\code{atmosTransmittance()} returns a
\bold{colorSpec} object with \code{\link{quantity}} equal to \code{'transmittance'}.
There is a spectrum in the object for each value in the vector \code{distance}.
The \code{specnames} are set to \code{sprintf("dist=\%gm",distance)}.\cr
The final transmittance is the product of the molecular transmittance
and the aerosol transmittance.
If both \code{molecules} and \code{aerosols} are \code{NULL},
then the final transmittance is identically 1;
the atmosphere has become a vacuum.

}

\note{
The Kruse model for \eqn{\alpha} as a function of \eqn{V_r}
is defined in 3 pieces.
For \eqn{0 \le V_r < 6000}, \eqn{\alpha = 0.585 * (V_r/1000)^{1/3}}.
For \eqn{6000 \le V_r < 50000}, \eqn{\alpha = 1.3}.
And for \eqn{V_r \ge} 50000, \eqn{\alpha = 1.6}.
So \eqn{\alpha} is increasing, but not strictly, and not continuously.
\eqn{V_r} is in meters.
See \emph{Kruse} and \emph{Kaushal}.

The built-in object \code{\link{atmosphere2003}} is transmittance along
an optical path that is \bold{NOT} horizontal,
and extends to outer space.
This is \emph{much} more complicated to calculate.
}

\references{
Angstrom, Anders.
On the atmospheric transmission of sun radiation and on dust in the air.
\emph{Geogr. Ann.},
no. 2. 1929.

Kaushal, H. and Jain, V.K. and Kar, S.
\bold{Free Space Optical Communication}.
Springer. 2017.

Koschmieder, Harald.
Theorie der horizontalen Sichtweite.
\emph{Beitrage zur Physik der Atmosphare}. 1924.
\bold{12}. pages 33-53.
  
P. W. Kruse, L. D. McGlauchlin, and R. B. McQuistan.
\bold{Elements of Infrared Technology: Generation, Transmission, and Detection}.
J. Wiley & Sons, New York, 1962.
}

\seealso{
\code{\link{solar.irradiance}},
\code{\link{specnames}}
}

\examples{
trans = atmosTransmittance( c(5,10,15,20,25)*1000 ) # 5 distances with atmospheric defaults

# verify that transmittance[550]=0.02 at dist=25000
plot( trans, legend='bottomright', log='y' )

# repeat, but this time assign alpha and beta explicitly
trans = atmosTransmittance( c(5,10,15,20,25)*1000, aero=list(alpha=1,beta=0.0001) )
}

\keyword{atmosphere}

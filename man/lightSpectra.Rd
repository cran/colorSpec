% File lightSpectra.Rd
\encoding{UTF-8}

\name{LightSpectra}

\title{compute standard light spectra}
\alias{lightSpectra}
\alias{daylightSpectra}
\alias{planckSpectra}
\alias{illuminantE}

\description{
Two families of standard illuminants that depend on temperature are the
Planckian spectra (black-body spectra), and daylight spectra.
For the daylight spectra, a smoothed version is available.
Illuminant E, a third and trivial spectrum, is also available.
}
\usage{
planckSpectra( temperature, wavelength=300:830, normalize=TRUE, c2=1.4388e7 )

daylightSpectra( temperature, wavelength=NULL, 
                    components=colorSpec::daylight1964, roundMs=FALSE )

illuminantE( energy=1, wavelength=380:780 )
}

\arguments{
\item{temperature}{a vector of temperatures, in Kelvin}

\item{wavelength}{a vector of wavelengths.
For \code{planckSpectra()} and  \code{illuminantE()} this is required. 
For \code{daylightSpectra()} this is optional. 
The default \code{wavelength=NULL} means to use the wavelengths in \code{components}, 
and otherwise \code{components} is resampled at the given \code{wavelength} vector. }

\item{normalize}{a logical value.  
If \code{TRUE} the Planck spectra are normalized to have value 1 at 560nm.
If \code{FALSE} then the quantity returned is radiant exitance with unit 
\eqn{W * m^{-2} * nm^{-1}}. }
       
\item{c2}{the value of \eqn{h*c/k} in Planck's law.
\eqn{h} is the Planck constant; \eqn{c} is the speed of light; and \eqn{k} is the Boltzmann constant.
The default value of 1.4388e7 \eqn{nm*K} was recommended by the CIE in 2005;
in 1986 the CIE recommended \code{c2}=1.438e7.
If \code{c2='calc'} then \code{c2} is calculated directly from the 3 physical constants.
}

\item{components}{a \bold{colorSpec} object with the daylight components \eqn{S_0, S_1}, and \eqn{S_2}.
The default is \code{\link{daylight1964}} and
a smoothed version \code{\link{daylight2013}} is also available.}

\item{roundMs}{a logical value.  The original CIE method for the daylight spectra requires
rounding intermediate coefficients M1 and M2 to 3 decimal places.
This rounding is necessary to reproduce the tabulated values in Table T.1 of the CIE
publication in \bold{References}.}

\item{energy}{a vector of energy levels}
}
    
\value{
For \code{planckSpectra()} and \code{daylightSpectra()} :\cr
A \bold{colorSpec} object with \code{\link{quantity}} equal to \code{'energy'},
and \code{\link{organization}} equal to \code{'matrix'} or \code{'vector'}.
The \code{specname}s are \code{PNNNN} or \code{DNNNN}
for \code{planckSpectra()} and \code{daylightSpectra()} respectively.\cr
The number of spectra in the object is the number of temperatures = \code{length(temperature)}.\cr
  
For \code{illuminantE()} :\cr
A \bold{colorSpec} object with \code{\link{quantity}} equal to \code{'energy'}. \cr
The number of spectra in the object is the number of energy levels = \code{length(energy)}.
}


\details{
For \code{daylightSpectra()} the valid range of temperatures is 4000 to 25000 K.
For a temperature outside this range the spectrum is set to all \code{NA}s.

The equations for \code{daylightSpectra()} and \code{planckSpectra()} are complex 
and can be found in the \bold{References}.

\code{IlluminantE()} is trivial - all constant energy.
}

\references{
Günther Wyszecki and W.S. Stiles.
\bold{Color Science : Concepts and Methods, Quantitative Data and Formulae}.
Second Edition.
Wiley-Interscience. 1982.
page 146.

CIE 15: Technical Report: Colorimetry, 3rd edition.
CIE 15:2004.
Table T.1, pp 30-32, and Note 5 on page 69.

Schanda, Janos.
CIE Colorimetry, in
\emph{Colorimetry: Understanding the CIE System}. 
Wiley Interscience. 2007.
p. 42.
}

\seealso{
\link{daylight}, \code{\link{resample}}, \code{\link{organization}}, \code{\link{quantity}},
\link{materialSpectra}
}

\keyword{light}
